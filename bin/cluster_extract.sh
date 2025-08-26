#!/usr/bin/env bash
set -euo pipefail

## TO DO - restrict the number of CPU threads used by mmseqs to 1

usage() {
  cat <<EOF
Usage: $0 \
  --rep            <repKey> \
  --q-db-dir       <Q_DB_DIR> \
  --ref-db-dir     <REF_DB_DIR> \
  --all-hits-dir   <ALL_HITS_DIR> \
  --members-tsv    <MEMBERS_NUMERIC_TSV> \
  --out-with-ref   <OUT_WITH_REF_DIR> \
  --out-query-only <OUT_QUERY_ONLY_DIR>
EOF
}

repKey=""                     # Cluster representative key (e.g., "31")
Q_DB_DIR=""                   # Query database directory
REF_DB_DIR=""                 # Reference database directory
ALL_HITS_DIR=""               # Global search results directory
MEMBERS_NUMERIC_TSV=""        # Cluster membership TSV
OUT_WITH_REF_DIR=""           # Output directory for FASTA files with reference sequences
OUT_QUERY_ONLY_DIR=""         # Output directory for FASTA files with only query sequences

while [[ $# -gt 0 ]]; do
  case "$1" in
    --rep)
      repKey="$2"; shift 2 ;;
    --q-db-dir)
      Q_DB_DIR="$2"; shift 2 ;;
    --ref-db-dir)
      REF_DB_DIR="$2"; shift 2 ;;
    --all-hits-dir)
      ALL_HITS_DIR="$2"; shift 2 ;;
    --members-tsv)
      MEMBERS_NUMERIC_TSV="$2"; shift 2 ;;
    --out-with-ref)
      OUT_WITH_REF_DIR="$2"; shift 2 ;;
    --out-query-only)
      OUT_QUERY_ONLY_DIR="$2"; shift 2 ;;
    -h|--help)
      usage; exit 0 ;;
    *)
      echo "Unknown option: $1" >&2; usage; exit 1 ;;
  esac
done

if [[ -z "${repKey}" || -z "${Q_DB_DIR}" || -z "${REF_DB_DIR}" || -z "${ALL_HITS_DIR}" || -z "${MEMBERS_NUMERIC_TSV}" || -z "${OUT_WITH_REF_DIR}" || -z "${OUT_QUERY_ONLY_DIR}" ]]; then
  echo "Missing required arguments" >&2
  usage
  exit 1
fi

# Resolve base paths
Q_DB_BASE=$(find "${Q_DB_DIR}" -name "*.lookup" | parallel -j1 "echo {/}" | sed 's/\.lookup//')
REF_DB_BASE=$(find "${REF_DB_DIR}" -name "*.lookup" | parallel -j1 "echo {/}" | sed 's/\.lookup//')
ALL_HITS_BASE=$(find "${ALL_HITS_DIR}" -name "*.dbtype" | parallel -j1 "echo {/}" | sed 's/\.dbtype//')

## Extract member keys for this cluster
awk -v r="${repKey}" '$1==r{print $2}' "${MEMBERS_NUMERIC_TSV}" \
  > "cluster_${repKey}.keys"

## Create sub-DB for cluster queries
mmseqs createsubdb \
  "cluster_${repKey}.keys" \
  "${Q_DB_DIR}"/"${Q_DB_BASE}" \
  "cluster_${repKey}_q" \
  --id-mode 0

## Subset global hits to these queries (subset by entry keys)
mmseqs createsubdb \
  "cluster_${repKey}.keys" \
  "${ALL_HITS_DIR}"/"${ALL_HITS_BASE}" \
  "cluster_${repKey}_hits" \
  --id-mode 0

## Export TSV to get target accessions (second column are targets)
mmseqs createtsv \
  "${Q_DB_DIR}"/"${Q_DB_BASE}" \
  "${REF_DB_DIR}"/"${REF_DB_BASE}" \
  "cluster_${repKey}_hits" \
  "cluster_${repKey}_hits.tsv"

awk '{print $2}' "cluster_${repKey}_hits.tsv" \
  | runiq - | sort -h \
  > "cluster_${repKey}_targets.acc" || true

## Map target accessions to numeric keys in SH_db
awk 'NR==FNR{t[$2]=$1; next} ($1 in t){print t[$1]}' \
  "${REF_DB_DIR}"/"${REF_DB_BASE}".lookup \
  "cluster_${repKey}_targets.acc" \
  > "cluster_${repKey}_t.keys" || true

## Create sub-DB for targets and export FASTAs
if [[ -s "cluster_${repKey}_t.keys" ]]; then
  mmseqs createsubdb \
    "cluster_${repKey}_t.keys" \
    "${REF_DB_DIR}"/"${REF_DB_BASE}" \
    "cluster_${repKey}_t"
  mmseqs convert2fasta \
    "cluster_${repKey}_t" \
    "cluster_${repKey}_t.fasta"
else
  : > "cluster_${repKey}_t.fasta"
fi

mmseqs convert2fasta \
  "cluster_${repKey}_q" \
  "cluster_${repKey}_q.fasta"

if [[ -s "cluster_${repKey}_t.fasta" ]]; then
  cat "cluster_${repKey}_q.fasta" "cluster_${repKey}_t.fasta" > "${OUT_WITH_REF_DIR%/}/cluster_${repKey}.fasta"
else
  cp "cluster_${repKey}_q.fasta" "${OUT_QUERY_ONLY_DIR%/}/cluster_${repKey}.fasta"
fi

## Cleanup intermediates for this cluster (keep final outputs)
rm -f "cluster_${repKey}_targets.acc" "cluster_${repKey}.keys"
rm -f "cluster_${repKey}_hits"*
rm -f "cluster_${repKey}_q"*
rm -f "cluster_${repKey}_t"*
