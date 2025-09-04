#!/usr/bin/env bash
set -euo pipefail

## TODO:
# - replace custom FASTA writing tricks with `exon`?
#   https://github.com/wheretrue/exon-duckdb
#   Current problem is that DuckDB's partitioned COPY doesn't guarantee ordering within partitions
#   which lead to errors in the FASTA format (e.g., sequence comes before header)
#   --> for now added explicit ordering hierarchy:


usage(){
  echo "Usage: $(basename "$0") [options]" 
  echo
  echo "Required:"
  echo "  -q, --query_seqs FILE         Input query sequences (FASTA)"
  echo "  -c, --cluster_members FILE    Query-cluster membership TSV (Cluster_ID, Query_ID)"
  echo "  -H, --top_hits FILE           Top N hits per query (TSV with header)"
  echo "  -d, --db_seqs FILE            Database of SH sequences (FASTA)"
  echo
  echo "Optional:"
  echo "  -T, --threads INT             Number of CPU threads"
  echo "  -M, --memory SIZE             Memory limit for DuckDB, e.g. 100G"
  echo "  -h, --help                    Show this help and exit"
  exit 1
}


## Initialize variables
QUERY_SEQS=""
CLUSTER_MEMBERS=""
TOP_HITS=""
DB_SEQS=""
THREADS=""
MEMORY=""

## Parse command-line options
LONG_TO_SHORT_ARGS=()
while [[ $# -gt 0 ]]; do
  case "$1" in
    --query_seqs=*)        LONG_TO_SHORT_ARGS+=("-q" "${1#*=}"); shift ;;
    --query_seqs)          shift; LONG_TO_SHORT_ARGS+=("-q" "${1:-}"); shift || true ;;
    --cluster_members=*)   LONG_TO_SHORT_ARGS+=("-c" "${1#*=}"); shift ;;
    --cluster_members)     shift; LONG_TO_SHORT_ARGS+=("-c" "${1:-}"); shift || true ;;
    --top_hits=*)          LONG_TO_SHORT_ARGS+=("-H" "${1#*=}"); shift ;;
    --top_hits)            shift; LONG_TO_SHORT_ARGS+=("-H" "${1:-}"); shift || true ;;
    --db_seqs=*)           LONG_TO_SHORT_ARGS+=("-d" "${1#*=}"); shift ;;
    --db_seqs)             shift; LONG_TO_SHORT_ARGS+=("-d" "${1:-}"); shift || true ;;
    --threads=*)           LONG_TO_SHORT_ARGS+=("-T" "${1#*=}"); shift ;;
    --threads)             shift; LONG_TO_SHORT_ARGS+=("-T" "${1:-}"); shift || true ;;
    --memory=*)            LONG_TO_SHORT_ARGS+=("-M" "${1#*=}"); shift ;;
    --memory)              shift; LONG_TO_SHORT_ARGS+=("-M" "${1:-}"); shift || true ;;
    --help|-h)             usage ;;
    --)                    shift; break ;;
    -*)                    LONG_TO_SHORT_ARGS+=("$1"); shift ;;
    *)                     LONG_TO_SHORT_ARGS+=("$1"); shift ;;
  esac
done
set -- "${LONG_TO_SHORT_ARGS[@]}"

OPTIND=1
while getopts ":q:c:H:d:T:M:h" opt; do
  case "$opt" in
    q) QUERY_SEQS="$OPTARG" ;;
    c) CLUSTER_MEMBERS="$OPTARG" ;;
    H) TOP_HITS="$OPTARG" ;;
    d) DB_SEQS="$OPTARG" ;;
    T) THREADS="$OPTARG" ;;
    M) MEMORY="$OPTARG" ;;
    h) usage ;;
    :) echo "Error: Option -$OPTARG requires an argument." >&2; usage ;;
    \?) echo "Error: Invalid option: -$OPTARG" >&2; usage ;;
  esac
done
shift $((OPTIND - 1))

if [[ $# -gt 0 ]]; then
  echo "Error: Unexpected positional arguments: $*" >&2
  usage
fi

## Validate input parameters
if [[ -z "$QUERY_SEQS" || -z "$CLUSTER_MEMBERS" || -z "$TOP_HITS" || -z "$DB_SEQS" ]]; then
    echo -e "Error: Missing required parameters!\n"
    usage
fi  

## Threads should be a positive integer
if [[ -n "$THREADS" && "$THREADS" -le 0 ]]; then
    echo -e "Error: Threads must be a positive integer!\n"
    usage
fi 

## View user-supplied parameters
echo -e "\nInput parameters:"
echo "Query sequences: $QUERY_SEQS"
echo "Cluster members: $CLUSTER_MEMBERS"
echo "Top hits: $TOP_HITS"
echo "Database sequences: $DB_SEQS"
if [[ -n "$THREADS" ]]; then
    echo "Threads: $THREADS"
fi
if [[ -n "$MEMORY" ]]; then
    echo "Memory: $MEMORY"
fi


## Convert data to parquet
echo -e "\n\nConverting data to parquet\n"

seqs_to_parquet (){
  # $1 = output file (parquet)

  duckdb :memory: \
    "CREATE TEMPORARY TABLE tbl AS SELECT * FROM read_csv('/dev/stdin',
      auto_detect = false,
      header = false,
      delim = '\t',
      columns = {
        'SeqID': 'VARCHAR',
        'Seq': 'VARCHAR'
      });
     COPY tbl TO '$1'
     (FORMAT PARQUET, COMPRESSION ZSTD, COMPRESSION_LEVEL 3);"
}


echo "..Query sequences"
seqkit fx2tab "$QUERY_SEQS" | sed 's/\t$//' | seqs_to_parquet "query_seqs.parquet"

echo "..Database sequences"
seqkit fx2tab "$DB_SEQS" | sed 's/\t$//' | seqs_to_parquet "db_seqs.parquet"

echo "..Cluster membership"
cat "$CLUSTER_MEMBERS" | duckdb :memory: \
    "CREATE TEMPORARY TABLE tbl AS SELECT * FROM read_csv('/dev/stdin',
      auto_detect = false,
      header = false,
      delim = '\t',
      columns = {
        'Cluster_ID': 'VARCHAR',
        'Query_ID': 'VARCHAR'
      });
     COPY tbl TO 'query_cluster_members.parquet'
     (FORMAT PARQUET, COMPRESSION ZSTD, COMPRESSION_LEVEL 3);"

echo "..Top hits"
cat "$TOP_HITS" | duckdb :memory: \
"COPY (
    SELECT qseqid, sseqid
    FROM read_csv_auto('/dev/stdin', delim='\t', header=true)
) TO 'top_hits.parquet' (
    FORMAT PARQUET, COMPRESSION ZSTD, COMPRESSION_LEVEL 3
);"



echo -e "\n\nMerge query and target sequences\n"

echo "..Preparing DuckDB command"

DUCKDB_COMMAND=""

## Add configuration settings (if provided)
if [[ -n "$THREADS" ]]; then
    DUCKDB_COMMAND+="
SET threads TO ${THREADS};
"
fi

if [[ -n "$MEMORY" ]]; then
    DUCKDB_COMMAND+="
SET memory_limit = '${MEMORY}';
"
fi

DUCKDB_COMMAND+="
COPY (
  WITH
  cm AS (
    SELECT Cluster_ID, Query_ID
    FROM 'query_cluster_members.parquet'
  ),
  th AS (
    SELECT qseqid, sseqid
    FROM 'top_hits.parquet'
  ),

  -- unique query members per cluster
  query_members AS (
    SELECT DISTINCT
      cm.Cluster_ID        AS ClusterID,
      'Query'              AS MemberType,
      cm.Query_ID::VARCHAR AS MemberID
    FROM cm
  ),

  -- unique sseqid hits per cluster (for all members)
  sseqid_hits AS (
    SELECT DISTINCT
      cm.Cluster_ID      AS ClusterID,
      'Ref'              AS MemberType,
      th.sseqid::VARCHAR AS MemberID
    FROM cm
    JOIN th
      ON th.qseqid = cm.Query_ID
    WHERE th.sseqid IS NOT NULL
  )

  SELECT DISTINCT ClusterID, MemberType, MemberID
  FROM (
    SELECT * FROM query_members
    UNION 
    SELECT * FROM sseqid_hits
  )
  ORDER BY ClusterID, MemberType, MemberID
) TO 'cluster_membership.txt'
  (FORMAT 'csv', DELIMITER '\t', HEADER TRUE);
"

## Execute the command
echo "..Executing DuckDB command"
duckdb -c "${DUCKDB_COMMAND}"




echo -e "\n\nEnriching cluster_membership with sequences\n"

echo "..Preparing Python script for cluster processing"

# Create Python script for processing clusters with proper FASTA ordering
cat > process_clusters.py << 'EOF'
#!/usr/bin/env python3
"""
Process cluster membership and generate FASTA files
Input files:
  cluster_membership.txt  (TSV with columns: ClusterID, MemberType, MemberID)
  db_seqs.parquet         (columns: SeqID, Seq)
  query_seqs.parquet      (columns: SeqID, Seq)
Output:
  FASTA files split by cluster and reference type (parquet partitions):
  - out_with_ref/clust_xxxx.fasta    (clusters containing Ref sequences)
  - out_query_only/clust_xxxx.fasta  (clusters with only Query sequences)
"""

import duckdb
import os
import sys
from pathlib import Path
from concurrent.futures import ThreadPoolExecutor, as_completed
from typing import Dict, List, Tuple, Optional
import logging

# Configure logging
logging.basicConfig(level=logging.INFO, format='%(asctime)s - %(levelname)s - %(message)s')
logger = logging.getLogger(__name__)

def setup_duckdb_connection(threads: Optional[int] = None, memory: Optional[str] = None) -> duckdb.DuckDBPyConnection:
    """Create and configure DuckDB connection."""
    con = duckdb.connect(':memory:')
    
    if threads:
        con.execute(f"SET threads TO {threads}")
        logger.info(f"Set DuckDB threads to {threads}")
    
    if memory:
        con.execute(f"SET memory_limit = '{memory}'")
        logger.info(f"Set DuckDB memory limit to {memory}")
    
    return con

def load_and_process_data(con: duckdb.DuckDBPyConnection) -> Dict:
    """Load cluster membership and sequence data, returning processed clusters."""
    
    logger.info("Loading and processing cluster membership and sequence data")
    
    # Join the core data
    query = """
    WITH
    -- Read cluster membership data
    cm AS (
        SELECT ClusterID, MemberType, MemberID
        FROM read_csv('cluster_membership.txt', 
            columns={'ClusterID': 'VARCHAR', 'MemberType': 'VARCHAR', 'MemberID': 'VARCHAR'},
            delim='\t', header=true
        )
    ),
    
    -- Combine sequence data from both sources with priority for query sequences
    all_sequences AS (
        SELECT SeqID, Seq, 0 AS priority FROM 'db_seqs.parquet'
        UNION ALL
        SELECT SeqID, Seq, 1 AS priority FROM 'query_seqs.parquet'
    ),
    
    -- Resolve duplicates by SeqID (prefer query sequences with higher priority)
    seq_map AS (
        SELECT SeqID, arg_max(Seq, priority) AS Seq
        FROM all_sequences
        GROUP BY SeqID
    )
    
    -- Join cluster membership with sequences
    SELECT cm.ClusterID, cm.MemberType, cm.MemberID, sm.Seq
    FROM cm
    INNER JOIN seq_map sm ON sm.SeqID = cm.MemberID
    ORDER BY cm.ClusterID, cm.MemberID
    """
    
    # Execute query and get results as DataFrame
    result_df = con.execute(query).fetchdf()
    logger.info(f"Loaded {len(result_df)} sequence records")
    
    # Group by cluster and determine cluster properties
    cluster_data = {}
    for _, row in result_df.iterrows():
        cluster_id = row['ClusterID']
        if cluster_id not in cluster_data:
            cluster_data[cluster_id] = {
                'members': [],
                'has_ref': False
            }
        
        # Check if this cluster has reference sequences
        if row['MemberType'] == 'Ref':
            cluster_data[cluster_id]['has_ref'] = True
            
        cluster_data[cluster_id]['members'].append({
            'member_id': row['MemberID'],
            'member_type': row['MemberType'],
            'sequence': row['Seq']
        })
    
    # Generate cluster names with zero-padding
    cluster_count = len(cluster_data)
    pad_width = len(str(cluster_count))
    
    # Create final clusters dict with proper naming
    clusters = {}
    for i, (cluster_id, data) in enumerate(sorted(cluster_data.items()), 1):
        cluster_name = f"clust_{i:0{pad_width}d}"
        output_dir = "out_with_ref" if data['has_ref'] else "out_query_only"
        
        clusters[cluster_name] = {
            'cluster_id': cluster_id,
            'has_ref': data['has_ref'],
            'output_dir': output_dir,
            'members': data['members']
        }
    
    logger.info(f"Organized data into {len(clusters)} clusters")
    return clusters

def write_fasta_file(cluster_name: str, cluster_data: Dict, base_dir: str) -> Tuple[str, int]:
    """Write FASTA file for a single cluster with guaranteed ordering."""
    
    output_dir = Path(base_dir) / cluster_data['output_dir']
    output_dir.mkdir(parents=True, exist_ok=True)
    
    fasta_file = output_dir / f"{cluster_name}.fasta"
    
    # Sort members by MemberID for consistent ordering
    members = sorted(cluster_data['members'], key=lambda x: x['member_id'])
    
    sequence_count = 0
    with open(fasta_file, 'w') as f:
        for member in members:
            # Write header line
            f.write(f">{member['member_id']}\n")
            # Write sequence line
            f.write(f"{member['sequence']}\n")
            sequence_count += 1
    
    return str(fasta_file), sequence_count

def process_clusters_parallel(clusters: Dict, threads: int = 1) -> Dict[str, int]:
    """Process all clusters in parallel using ThreadPoolExecutor."""
    
    logger.info(f"Processing {len(clusters)} clusters using {threads} threads")
    
    # Create output directories
    Path("out_with_ref").mkdir(exist_ok=True)
    Path("out_query_only").mkdir(exist_ok=True)
    
    results = {'files_created': 0, 'sequences_written': 0, 'with_ref': 0, 'query_only': 0}
    
    with ThreadPoolExecutor(max_workers=threads) as executor:
        # Submit all cluster processing tasks
        future_to_cluster = {
            executor.submit(write_fasta_file, cluster_name, cluster_data, "."): cluster_name
            for cluster_name, cluster_data in clusters.items()
        }
        
        # Collect results
        for future in as_completed(future_to_cluster):
            cluster_name = future_to_cluster[future]
            try:
                fasta_file, seq_count = future.result()
                results['files_created'] += 1
                results['sequences_written'] += seq_count
                
                if 'out_with_ref' in fasta_file:
                    results['with_ref'] += 1
                else:
                    results['query_only'] += 1
                    
                logger.debug(f"Created {fasta_file} with {seq_count} sequences")
                
            except Exception as e:
                logger.error(f"Error processing cluster {cluster_name}: {e}")
                raise
    
    return results

def main():
    """Main function for cluster processing."""
    
    # Parse command line arguments (threads, memory)
    threads = None
    memory = None
    
    if len(sys.argv) > 1 and sys.argv[1].strip():
        try:
            threads = int(sys.argv[1])
        except ValueError:
            logger.warning(f"Invalid threads value '{sys.argv[1]}', using default")
            
    if len(sys.argv) > 2 and sys.argv[2].strip():
        memory = sys.argv[2]
    
    try:
        # Setup DuckDB connection
        con = setup_duckdb_connection(threads, memory)
        
        # Load and process data
        clusters = load_and_process_data(con)
        
        # Process clusters in parallel
        processing_threads = threads if threads else 1
        results = process_clusters_parallel(clusters, processing_threads)
        
        # Report results
        logger.info("="*50)
        logger.info("FASTA generation complete")
        logger.info(f"Files created: {results['files_created']}")
        logger.info(f"Sequences written: {results['sequences_written']}")
        logger.info(f"Clusters with references: {results['with_ref']}")
        logger.info(f"Query-only clusters: {results['query_only']}")
        logger.info("="*50)
        
        # Close connection
        con.close()
        
    except Exception as e:
        logger.error(f"Error in main processing: {e}")
        sys.exit(1)

if __name__ == "__main__":
    main()
EOF

## Execute the Python script
echo "..Executing Python cluster processing script"
python process_clusters.py "${THREADS:-}" "${MEMORY:-}"

## Clean up temporary Python script
rm -f process_clusters.py


echo -e "\n\nValidating sequence counts"

## Input seqs
seq_count_in=$(rg -c "^>" "$QUERY_SEQS" 2>/dev/null || echo "0")

## Output seqs (query + ref)
seq_count_out1=$(find out_with_ref -name "*.fasta" | parallel -j1 "cat {}" | seqkit stat -T | awk 'NR==2{print $4}')
seq_count_out2=$(find out_query_only -name "*.fasta" | parallel -j1 "cat {}" | seqkit stat -T | awk 'NR==2{print $4}')

## Number of refs
seq_count_refs=$(awk '$2 == "Ref" { count++ } END { print count }' cluster_membership.txt)

## Number of valid query seqs in output
seq_count_out=$((seq_count_out1 + seq_count_out2 - seq_count_refs))

echo "  Sequences in:  $seq_count_in"
echo "  Sequences out: $seq_count_out"

if [ "$seq_count_in" -ne "$seq_count_out" ]; then
  echo "ERROR: Number of sequences in and out do not match"
  exit 1
fi
