#!/usr/bin/env bash
set -euo pipefail

## Convert data to parquet
echo -e "Converting data to parquet\n"

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


echo -e "..Query sequences\n"
seqkit fx2tab test_inp.fasta | sed 's/\t$//' | seqs_to_parquet "query_seqs.parquet"

echo -e "..Database sequences\n"
seqkit fx2tab sanger_refs_sh_full.fasta | sed 's/\t$//' | seqs_to_parquet "db_seqs.parquet"

echo -e "..Cluster membership\n"
cat query_cluster_membership.tsv | duckdb :memory: \
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

echo -e "..Top hits\n"
cat top_hits.tsv | duckdb :memory: \
"COPY (
    SELECT qseqid, sseqid
    FROM read_csv_auto('/dev/stdin', delim='\t', header=true)
) TO 'top_hits.parquet' (
    FORMAT PARQUET, COMPRESSION ZSTD, COMPRESSION_LEVEL 3
);"



echo -e "\n\nMerge query and target sequences\n"

echo -e "..Preparing DuckDB command\n"

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
      cm.Cluster_ID    AS ClusterID,
      'Query'::VARCHAR AS MemberType,
      cm.Query_ID      AS MemberID
    FROM cm
  ),

  -- unique sseqid hits per cluster (for all members)
  sseqid_members AS (
    SELECT DISTINCT
      cm.Cluster_ID  AS ClusterID,
      'Ref'::VARCHAR AS MemberType,
      th.sseqid      AS MemberID
    FROM cm
    JOIN th
      ON th.qseqid = cm.Query_ID
    WHERE th.sseqid IS NOT NULL
  )

  SELECT DISTINCT ClusterID, MemberType, MemberID
  FROM (
    SELECT * FROM query_members
    UNION ALL
    SELECT * FROM sseqid_members
  )
  ORDER BY ClusterID, MemberType, MemberID
) TO 'cluster_membership.txt'
  (FORMAT 'csv', DELIMITER '\t', HEADER TRUE);
"

## Execute the command
echo -e "..Executing DuckDB command\n"
duckdb -c "${DUCKDB_COMMAND}"
                             

