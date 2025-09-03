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

echo -e "\n\nEnriching cluster_membership with sequences\n"

echo -e "..Preparing DuckDB command\n"

SEQS_COMMAND+="
-- Input files:
--   cluster_membership.txt  (TSV with columns: ClusterID, MemberType, MemberID)
--   db_seqs.parquet         (columns: SeqID, Seq)
--   query_seqs.parquet      (columns: SeqID, Seq)
WITH
cm AS (
  SELECT
    CAST(ClusterID  AS VARCHAR) AS ClusterID,
    CAST(MemberType AS VARCHAR) AS MemberType,
    CAST(MemberID   AS VARCHAR) AS MemberID
  FROM read_csv_auto('cluster_membership.txt', delim='\t', header=true)
),

-- Load sequences dictionaries (+ tag with a priority so we can break ties)
db AS (
  SELECT CAST(SeqID AS VARCHAR) AS SeqID,
         CAST(Seq   AS VARCHAR) AS Seq,
         0 AS src_order
  FROM 'db_seqs.parquet'
),
qry AS (
  SELECT CAST(SeqID AS VARCHAR) AS SeqID,
         CAST(Seq   AS VARCHAR) AS Seq,
         1 AS src_order
  FROM 'query_seqs.parquet'
),
all_seqs AS (
  SELECT * FROM db
  UNION ALL
  SELECT * FROM qry
),

-- Resolve duplicates by SeqID (prefer rows with higher src_order = query_seqs)
seq_map AS (
  SELECT
    SeqID,
    arg_max(Seq, src_order) AS Seq
  FROM all_seqs
  GROUP BY SeqID
)

-- Add Seq to the cluster membership table
SELECT
  cm.ClusterID,
  cm.MemberType,
  cm.MemberID,
  sm.Seq
FROM cm
LEFT JOIN seq_map sm
  ON sm.SeqID = cm.MemberID
ORDER BY cm.ClusterID, cm.MemberType, cm.MemberID;
"

## Execute the command
echo -e "..Executing DuckDB command\n"
duckdb -c "${SEQS_COMMAND}"
                             
