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
     COPY tbl TO 'cluster_membership.parquet'
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


## Start the SQL command
echo -e "\nPreparing DuckDB command"

DUCKDB_COMMAND=""

DUCKDB_COMMAND+="
WITH
cm AS (
  SELECT Cluster_ID, Query_ID
  FROM 'cluster_membership.parquet'
),
th AS (
  SELECT qseqid, sseqid
  FROM 'top_hits.parquet'
),

-- unique member IDs per cluster
members AS (
  SELECT
    Cluster_ID,
    list(DISTINCT Query_ID) AS member_query_ids
  FROM cm
  GROUP BY Cluster_ID
),

-- union of all sseqids (matches) for every member in the cluster
hits AS (
  SELECT
    cm.Cluster_ID,
    -- drop NULLs using FILTER; DISTINCT removes duplicates
    list(DISTINCT th.sseqid) FILTER (WHERE th.sseqid IS NOT NULL) AS cluster_sseqids
  FROM cm
  LEFT JOIN th
    ON th.qseqid = cm.Query_ID            -- cast if types differ, e.g. th.qseqid::TEXT = cm.Query_ID::TEXT
  GROUP BY cm.Cluster_ID
)

SELECT
  m.Cluster_ID,
  -- sort the lists to make output deterministic
  list_sort(m.member_query_ids)   AS member_query_ids,
  list_sort(h.cluster_sseqids)    AS cluster_sseqids
FROM members m
LEFT JOIN hits h USING (Cluster_ID)
ORDER BY m.Cluster_ID;
"

## Execute the command
echo -e "\nExecuting DuckDB command"

duckdb -c "${DUCKDB_COMMAND}"
                             

