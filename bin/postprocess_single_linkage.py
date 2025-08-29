#!/usr/bin/env python

"""
Post-process single-linkage clustering outputs 
for a single (compound) cluster and similarity threshold (TH) series.

Inputs per cluster:
- clusters_dir:       directory containing calc_distm_out files named like `mx_cluster_<ClusterID>_out_{TH}`
                      where TH in {0.030, 0.025, 0.020, 0.015, 0.010, 0.005}
- cluster_membership: tab-delimited file (no header) with ALL clusters across the run:
                      ClusterID, MemberType (Query|Ref), MemberID
                      The script subsets to the current ClusterID parsed from filenames
- best_hits_tsv:      TSV with columns: query, target, evalue, bits, alnlen, pident, qcov, tcov
- centroid2sh:        `centroid2sh_mappings.txt` file from the SH database

Output:
- Single TSV file with header and the following columns:
  query\tbest_ref\tstatus\tsh_code\textra\tthreshold
  (rows are the same as produced by `analyse_usearch_output_sh.py` of the original pipeline,
   with an added `threshold` column)

Notes:
- We detect if query clusters with its best reference at threshold TH (status "present")
- If not present, we check the query's cluster:
    only queries -> "new cluster"
    single query -> "singleton"
- If a cluster contains references but not the best-hit reference, 
    mark "present" if any ref in the cluster maps to the same SH as the best hit;
    otherwise keep strict behavior.
 - centroid2sh is a prebuilt mapping (per threshold code 1..6) from reference centroid IDs to SH codes, 
     it is used to translate the best-hit reference into the correct SH at each TH
     (the same centroid can belong to different SH codes at different cutoffs, 
      so this mapping must be threshold-specific)

Usage:
  postprocess_single_linkage.py \
    --clusters-dir       calc_distm_out \
    --cluster-membership cluster_membership.txt \
    --best-hits          best_hits.tsv \
    --centroid2sh        centroid2sh_mappings.txt \
    --output             pp_matches/matches.tsv
"""

import argparse
import csv
import os
from pathlib import Path
import re


# Threshold names
THRESHOLDS = ["0.030", "0.025", "0.020", "0.015", "0.010", "0.005"]


def load_cluster_members(cluster_membership_fp: Path, cluster_id: str):
    query_ids = set()
    ref_ids = set()
    if not cluster_membership_fp.exists():
        return query_ids, ref_ids
    with open(cluster_membership_fp) as f:
        reader = csv.reader(f, delimiter="\t")
        for row in reader:
            if len(row) < 3:
                continue
            cid, mtype, mid = row[0], row[1], row[2]
            if cid != cluster_id:
                continue
            if mtype.lower().startswith("q"):
                query_ids.add(mid)
            elif mtype.lower().startswith("r"):
                ref_ids.add(mid)
    return query_ids, ref_ids


def load_clusters(cluster_file: Path):
    # cluster_file has rows: cluster_id\tseq_id
    clusters = {}
    if not cluster_file.exists():
        return clusters
    with open(cluster_file) as f:
        reader = csv.reader(f, delimiter="\t")
        for row in reader:
            if not row:
                continue
            cid, sid = row[0], row[1]
            clusters.setdefault(cid, []).append(sid)
    return clusters


def invert_membership(clusters):
    m = {}
    for cid, members in clusters.items():
        for s in members:
            m[s] = cid
    return m


def main():
    ap = argparse.ArgumentParser()
    ap.add_argument("--clusters-dir",       required=True, help="Directory with calc_distm_out/*.txt per threshold")
    ap.add_argument("--cluster-membership", required=True, help="TSV (no header): ClusterID, MemberType(Query|Ref), MemberID")
    ap.add_argument("--best-hits",          required=True, help="Output of MMseqs convertalis with best hits (TSV without header: query,target,evalue,bits,alnlen,pident,qcov,tcov)")
    ap.add_argument("--centroid2sh",        required=True, help="centroid2sh_mappings.txt for mapping ref -> SH per threshold code 1..6")
    ap.add_argument("--output",             required=True, help="Output TSV path (single file)")
    args = ap.parse_args()

    out_path = Path(args.output)
    out_path.parent.mkdir(parents=True, exist_ok=True)

    # Determine cluster ID from any threshold filename in the directory
    any_file = None
    for p in sorted(Path(args.clusters_dir).glob("*_out_*")):
        any_file = p
        break
    if any_file is None:
        raise RuntimeError(f"No clustering output files found in {args.clusters_dir}")
    m = re.search(r"mx_cluster_(.+?)_out_", any_file.stem)
    if not m:
        raise RuntimeError(f"Cannot parse cluster ID from filename: {any_file.name}")
    cluster_id = m.group(1)

    # Load query/ref members for this cluster
    query_ids, ref_ids = load_cluster_members(Path(args.cluster_membership), cluster_id)

    # best hits (select best per query by evalue asc, bits desc, pident desc)
    def parse_best_hits(tsv_path: str):
        with open(tsv_path) as f:
            reader = csv.reader(f, delimiter="\t")
            first = next(reader, None)
            if first is None:
                return {}, {}
            # Detect header
            columns = [c.strip().lower() for c in first]
            has_header = "query" in columns and "target" in columns
            if has_header:
                col_idx = {name: i for i, name in enumerate(columns)}
            else:
                # Assume default MMseqs convertalis order
                col_idx = {"query": 0, "target": 1, "evalue": 2, "bits": 3, "alnlen": 4, "pident": 5, "qcov": 6, "tcov": 7}
            best = {}
            if not has_header and first:
                # Treat first as data row
                rows_iter = [first]
            else:
                rows_iter = []
            # Chain remaining rows
            for row in rows_iter + list(reader):
                try:
                    qid = row[col_idx["query"]]
                    tid = row[col_idx["target"]]
                    evalue = float(row[col_idx.get("evalue", 2)])
                    bits = float(row[col_idx.get("bits", 3)])
                    pident = float(row[col_idx.get("pident", 5)])
                except Exception:
                    continue
                prev = best.get(qid)
                key = (evalue, -bits, -pident)
                if prev is None or key < prev[0]:
                    best[qid] = (key, tid, pident)
            q2best_ref = {q: v[1] for q, v in best.items()}
            q2pident = {q: v[2] for q, v in best.items()}
            return q2best_ref, q2pident

    q2best, q2pident = parse_best_hits(args.best_hits)

    # centroid2sh mappings: columns threshold_code(1..6) keyed later, for quick lookup per ref id
    # File format in old pipeline: fields[2] is ref id, fields[0] is SH code, fields[1] is threshold code
    # We'll build: map_by_th[th_code][ref_id] = SH_code
    map_by_th = {"1": {}, "2": {}, "3": {}, "4": {}, "5": {}, "6": {}}
    with open(args.centroid2sh) as f:
        r = csv.reader(f, delimiter="\t")
        for row in r:
            if len(row) < 3:
                continue
            sh_code, th_code, ref_id = row[0], row[1], row[2]
            if th_code in map_by_th:
                map_by_th[th_code][ref_id] = sh_code

    # threshold -> coded mapping like in `parse_matches_sh.pl`
    # legacy keys (03, 025, ...) and new decimal keys (0.030, 0.025, ...)
    th_to_code = {
        "03":  "1",  "0.030": "1",
        "02":  "2",  "0.020": "2",
        "01":  "3",  "0.010": "3",
        "025": "4",  "0.025": "4",
        "015": "5",  "0.015": "5",
        "005": "6",  "0.005": "6",
    }

    # Write results (all thresholds are aggregated into a single file)
    with open(out_path, "w") as o:
        w = csv.writer(o, delimiter="\t", lineterminator="\n")
        # Header
        w.writerow(["query", "best_ref", "status", "sh_code", "extra", "threshold"])

        for th in THRESHOLDS:
            th_code = th_to_code[th]
            # goclust produces files named like mx_*.txt_out_TH with lines: clusterId\tseqId
            files = sorted(Path(args.clusters_dir).glob(f"*_out_{th}"))
            clusters = {}
            for f in files:
                part = load_clusters(f)
                for cid, members in part.items():
                    # ensure cluster IDs are unique per cluster file by prefixing
                    key = f"{f.stem}_{cid}"
                    clusters[key] = members
            membership = invert_membership(clusters)

            for q in sorted(query_ids):
                best_ref = q2best.get(q, "")
                # default values
                status = "missing"
                sh_code = ""
                extra = ""

                cid = membership.get(q)
                if cid is None:
                    w.writerow([q, best_ref, status, sh_code, extra, th])
                    continue
                members = clusters[cid]

                members_set = set(members)
                ref_in_cluster = members_set & ref_ids
                has_best_ref = best_ref in ref_in_cluster if best_ref else False
                best_sh = map_by_th[th_code].get(best_ref, "") if best_ref else ""

                if has_best_ref:
                    status = "present"
                    # Map best_ref to SH at this threshold code
                    sh_code = best_sh
                    extra = ""
                else:
                    # If no refs, check singleton vs new cluster
                    if not ref_in_cluster:
                        only_q = all(m in query_ids for m in members)
                        if only_q:
                            if len(members) == 1:
                                status = "singleton"
                                # Keep sh_code as best_sh (context from best hit)
                                sh_code = best_sh
                            else:
                                status = "new cluster"
                                # join members with space as in legacy
                                extra = " ".join(members)
                                # Provide SH context from best hit if available
                                sh_code = best_sh
                        else:
                            # mixed but missing best ref; consider present if any ref maps to same SH as best_ref
                            if best_ref:
                                best_sh = map_by_th[th_code].get(best_ref, None)
                                if best_sh is not None:
                                    for r in ref_in_cluster:
                                        if map_by_th[th_code].get(r, None) == best_sh:
                                            status = "present"
                                            sh_code = best_sh
                                            break
                            # otherwise leave as missing
                    else:
                        # There are references but not best_ref; mark present if any ref present
                        status = "present"
                        # pick SH of any ref deterministically (sorted)
                        pick = sorted(ref_in_cluster)[0]
                        sh_code = map_by_th[th_code].get(pick, "")

                w.writerow([q, best_ref, status, sh_code, extra, th])


if __name__ == "__main__":
    main()

