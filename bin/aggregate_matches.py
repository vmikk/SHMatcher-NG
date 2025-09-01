#!/usr/bin/env python

import argparse
import csv
from pathlib import Path


THRESHOLDS = ["03", "025", "02", "015", "01", "005"]
THRESHOLDS_DEC = ["0.030", "0.025", "0.020", "0.015", "0.010", "0.005"]
TH_ALIAS = {"03": "0.030", "025": "0.025", "02": "0.020", "015": "0.015", "01": "0.010", "005": "0.005"}


def read_matches_aggregated(files):
    # Input = a single `matches.tsv` with header [query, best_ref, status, sh_code, extra, threshold]
    per_th = {th: {} for th in THRESHOLDS}
    all_queries = set()
    for f in files:
        with open(f) as h:
            r = csv.reader(h, delimiter="\t")
            header = next(h, None)
            if header is None:
                continue
            for row in r:
                # Skip empty/whitespace-only or malformed rows
                if not row or len(row) < 6:
                    continue
                q = row[0]
                all_queries.add(q)
                # Normalize threshold to legacy key if decimal
                th = row[5]
                for k, v in TH_ALIAS.items():
                    if th == v:
                        th = k
                        break
                per_th[th][q] = row
    return per_th, sorted(all_queries)


def load_map_file(fp, key_idx, val_idx):
    m = {}
    with open(fp) as f:
        r = csv.reader(f, delimiter="\t")
        for row in r:
            if len(row) <= max(key_idx, val_idx):
                continue
            m[row[key_idx]] = row[val_idx]
    return m


def main():
    ap = argparse.ArgumentParser(description="Aggregate per-cluster matches into matches_out_all.csv compatible with return_common_taxonomy.py")
    ap.add_argument("--matches-dirs",   nargs="+",     help="Directories or files to search for matches (all .tsv files will be read)")
    ap.add_argument("--best-hits",      required=True, help="Output of MMseqs convertalis (TSV without herader: query,target,evalue,bits,alnlen,pident,qcov,tcov)")
    ap.add_argument("--shs-file",       required=True, help="shs_out.txt (SH_code \t taxonomy)")
    ap.add_argument("--sh2compound",    required=True, help="sh2compound_mapping.txt (SH_code \t UCL_code)")
    ap.add_argument("--compounds-file", required=True, help="compounds_out.txt (.. \t UCL_code \t taxonomy)")
    ap.add_argument("--out",            required=True, help="Output CSV path: matches_out_all.csv")
    args = ap.parse_args()

    # Collect aggregated matches files
    agg_files = []
    for root in args.matches_dirs:
        p = Path(root)
        if p.is_dir():
            # Read all .tsv files under the specified directories
            agg_files.extend(sorted(str(x) for x in p.rglob("*.tsv") if x.is_file()))
        elif p.is_file():
            agg_files.append(str(p))

    per_th, all_queries = read_matches_aggregated(agg_files)

    # Best hits mapping (parse MMseqs convertalis; pick best per query)
    def parse_best_hits(tsv_path: str):
        with open(tsv_path) as f:
            reader = csv.reader(f, delimiter="\t")
            first = next(reader, None)
            if first is None:
                return {}, {}
            columns = [c.strip().lower() for c in first]
            has_header = "query" in columns and "target" in columns
            if has_header:
                col_idx = {name: i for i, name in enumerate(columns)}
            else:
                col_idx = {"query": 0, "target": 1, "evalue": 2, "bits": 3, "alnlen": 4, "pident": 5, "qcov": 6, "tcov": 7}
            best = {}
            if not has_header and first:
                rows_iter = [first]
            else:
                rows_iter = []
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

    # Helper to extract SH code (e.g., SH0989441.10FU) from best-hit target id
    def extract_sh_code_from_target(target_id: str):
        if not target_id:
            return ""
        parts = target_id.split("_", 1)
        if len(parts) == 2 and parts[1].startswith("SH"):
            return parts[1]
        return ""

    # Taxonomies
    sh2tax = load_map_file(args.shs_file, 0, 1)
    sh2ucl = load_map_file(args.sh2compound, 0, 1)
    # compounds_out.txt: columns [..., UCL, TAX]
    ucl2tax = {}
    with open(args.compounds_file) as f:
        r = csv.reader(f, delimiter="\t")
        for row in r:
            if len(row) >= 3:
                ucl2tax[row[1]] = row[2]

    with open(args.out, "w") as o:
        # Header mirrors merge_matches.py final header
        header = [
            "seq_id_tmp",
            "seq_accno",
            "status (3.0)", "SH code (3.0)", "SH/compound taxonomy (3.0)",
            "status (2.5)", "SH code (2.5)", "SH/compound taxonomy (2.5)",
            "status (2.0)", "SH code (2.0)", "SH/compound taxonomy (2.0)",
            "status (1.5)", "SH code (1.5)", "SH/compound taxonomy (1.5)",
            "status (1.0)", "SH code (1.0)", "SH/compound taxonomy (1.0)",
            "status (0.5)", "SH code (0.5)", "SH/compound taxonomy (0.5)",
            "compound_cl_code (0.5)", "Compound taxonomy (0.5)",
            "Matched sequence", "Similarity percentage"
        ]
        o.write("\t".join(header) + "\n")

        for q in all_queries:
            seq_id_tmp = q  # keep as-is; no renaming
            seq_accno = q

            row_out = [seq_id_tmp, seq_accno]

            compound_ucl = ""
            compound_tax = ""

            for th in ["03", "025", "02", "015", "01", "005"]:
                m = per_th[th].get(q)
                status = ""
                sh_code = ""
                tax = ""
                if m:
                    status_raw = m[2]
                    best_ref = m[1]
                    sh_maybe = m[3]
                    if status_raw == "present":
                        status = "present_in"
                        # Determine SH code: prefer explicit SH in input; otherwise derive from best hit
                        if sh_maybe and sh_maybe.startswith("SH"):
                            sh_code = sh_maybe
                        else:
                            sh_code = extract_sh_code_from_target(best_ref)
                        tax = sh2tax.get(sh_code, "")
                        # Derive compound from SH if available
                        if not compound_ucl and sh_code and sh_code in sh2ucl:
                            compound_ucl = sh2ucl[sh_code]
                            compound_tax = ucl2tax.get(compound_ucl, "")
                    elif status_raw == "new cluster":
                        status = "new_sh_in"
                        if sh_maybe:
                            # Input may contain either SH or UCL; convert to UCL and taxonomy
                            if sh_maybe.startswith("UCL"):
                                ucl = sh_maybe
                            else:
                                ucl = sh2ucl.get(sh_maybe, "")
                            tax = ucl2tax.get(ucl, "")
                            if not compound_ucl and ucl:
                                compound_ucl = ucl
                                compound_tax = tax
                    elif status_raw == "singleton":
                        status = "new_singleton_in"
                        if sh_maybe:
                            if sh_maybe.startswith("UCL"):
                                ucl = sh_maybe
                            else:
                                ucl = sh2ucl.get(sh_maybe, "")
                            tax = ucl2tax.get(ucl, "")
                            if not compound_ucl and ucl:
                                compound_ucl = ucl
                                compound_tax = tax
                row_out.extend([status, sh_code, tax])

            # compound columns (from 0.5 threshold context)
            row_out.extend([compound_ucl, compound_tax])

            matched_seq = q2best.get(q, "")
            pident = q2pident.get(q, "")
            row_out.extend([matched_seq, pident])

            o.write("\t".join(map(str, row_out)) + "\n")


if __name__ == "__main__":
    main()

