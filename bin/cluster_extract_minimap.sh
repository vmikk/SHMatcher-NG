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
echo -e "..Executing DuckDB command\n"
duckdb -c "${DUCKDB_COMMAND}"





echo -e "\n\nEnriching cluster_membership with sequences\n"

echo -e "..Preparing DuckDB command\n"

SEQS_COMMAND+="
-- Input files:
--   cluster_membership.txt  (TSV with columns: ClusterID, MemberType, MemberID)
--   db_seqs.parquet         (columns: SeqID, Seq)
--   query_seqs.parquet      (columns: SeqID, Seq)
-- Output:
--   FASTA files split by cluster and reference type (parquet partitions):
--   - out_with_ref/cluster_name=clust_xxxx/data_0.fasta    (clusters containing Ref sequences)
--   - out_query_only/cluster_name=clust_xxxx/data_0.fasta  (clusters with only Query sequences)

COPY (
  WITH
  -- Read cluster membership data
  cm AS (
    SELECT
      ClusterID,
      MemberType,
      MemberID
    FROM read_csv('cluster_membership.txt', 
      columns={'ClusterID': 'VARCHAR', 'MemberType': 'VARCHAR', 'MemberID': 'VARCHAR'},
      delim='\t', 
      header=true
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
    SELECT
      SeqID,
      arg_max(Seq, priority) AS Seq
    FROM all_sequences
    GROUP BY SeqID
  ),

  -- Join cluster membership with sequences
  members_with_seqs AS (
    SELECT 
      cm.ClusterID,
      cm.MemberType,
      cm.MemberID,
      sm.Seq
    FROM cm
    INNER JOIN seq_map sm ON sm.SeqID = cm.MemberID
  ),

  -- Determine cluster type (has_ref: true if cluster contains Ref sequences)
  cluster_types AS (
    SELECT 
      ClusterID,
      bool_or(MemberType = 'Ref') AS has_ref
    FROM members_with_seqs
    GROUP BY ClusterID
  ),

  -- Create dense numeric cluster indices with zero-padding
  cluster_numbering AS (
    SELECT
      ClusterID,
      has_ref,
      ROW_NUMBER() OVER (ORDER BY ClusterID) AS cluster_no
    FROM cluster_types
  ),
  padding_params AS (
    SELECT LENGTH(CAST(MAX(cluster_no) AS VARCHAR)) AS pad_width
    FROM cluster_numbering
  ),

  -- Generate cluster directory names and output paths
  enriched_members AS (
    SELECT
      m.ClusterID,
      m.MemberType,
      m.MemberID,
      m.Seq,
      cn.has_ref,
      'clust_' || lpad(CAST(cn.cluster_no AS VARCHAR), CAST(pp.pad_width AS INTEGER), '0') AS cluster_name,
      CASE 
        WHEN cn.has_ref THEN 'out_with_ref'
        ELSE 'out_query_only'
      END AS output_dir
    FROM members_with_seqs m
    JOIN cluster_numbering cn USING (ClusterID)
    CROSS JOIN padding_params pp
  ),

  -- Generate FASTA format lines for clusters WITH references
  fasta_lines_with_ref AS (
    SELECT 
      cluster_name,
      MemberID,
      0 AS line_order,
      '>' || MemberID AS line_content
    FROM enriched_members
    WHERE output_dir = 'out_with_ref'
    
    UNION ALL
    
    SELECT 
      cluster_name,
      MemberID,
      1 AS line_order,
      Seq AS line_content
    FROM enriched_members
    WHERE output_dir = 'out_with_ref'
  ),
  
  -- Generate FASTA format lines for clusters QUERY-ONLY
  fasta_lines_query_only AS (
    SELECT 
      cluster_name,
      MemberID,
      0 AS line_order,
      '>' || MemberID AS line_content
    FROM enriched_members
    WHERE output_dir = 'out_query_only'
    
    UNION ALL
    
    SELECT 
      cluster_name,
      MemberID,
      1 AS line_order,
      Seq AS line_content
    FROM enriched_members
    WHERE output_dir = 'out_query_only'
  )

  -- First export: clusters with references
  SELECT cluster_name, line_content
  FROM fasta_lines_with_ref
  ORDER BY cluster_name, MemberID, line_order
)
TO 'out_with_ref' (
  FORMAT CSV,
  HEADER false,
  PARTITION_BY cluster_name,
  FILE_EXTENSION 'fasta'
);

-- Second export: query-only clusters
COPY (
  WITH
  -- Re-use the same CTEs as above (need to repeat them for second COPY)
  cm AS (
    SELECT ClusterID, MemberType, MemberID
    FROM read_csv('cluster_membership.txt', 
      columns={'ClusterID': 'VARCHAR', 'MemberType': 'VARCHAR', 'MemberID': 'VARCHAR'},
      delim='\t', header=true
    )
  ),
  all_sequences AS (
    SELECT SeqID, Seq, 0 AS priority FROM 'db_seqs.parquet'
    UNION ALL
    SELECT SeqID, Seq, 1 AS priority FROM 'query_seqs.parquet'
  ),
  seq_map AS (
    SELECT SeqID, arg_max(Seq, priority) AS Seq
    FROM all_sequences GROUP BY SeqID
  ),
  members_with_seqs AS (
    SELECT cm.ClusterID, cm.MemberType, cm.MemberID, sm.Seq
    FROM cm INNER JOIN seq_map sm ON sm.SeqID = cm.MemberID
  ),
  cluster_types AS (
    SELECT ClusterID, bool_or(MemberType = 'Ref') AS has_ref
    FROM members_with_seqs GROUP BY ClusterID
  ),
  cluster_numbering AS (
    SELECT ClusterID, has_ref, ROW_NUMBER() OVER (ORDER BY ClusterID) AS cluster_no
    FROM cluster_types
  ),
  padding_params AS (
    SELECT LENGTH(CAST(MAX(cluster_no) AS VARCHAR)) AS pad_width
    FROM cluster_numbering
  ),
  enriched_members AS (
    SELECT m.ClusterID, m.MemberType, m.MemberID, m.Seq, cn.has_ref,
           'clust_' || lpad(CAST(cn.cluster_no AS VARCHAR), CAST(pp.pad_width AS INTEGER), '0') AS cluster_name,
           CASE WHEN cn.has_ref THEN 'out_with_ref' ELSE 'out_query_only' END AS output_dir
    FROM members_with_seqs m
    JOIN cluster_numbering cn USING (ClusterID)
    CROSS JOIN padding_params pp
  ),
  fasta_lines_query_only AS (
    SELECT cluster_name, MemberID, 0 AS line_order, '>' || MemberID AS line_content
    FROM enriched_members WHERE output_dir = 'out_query_only'
    UNION ALL
    SELECT cluster_name, MemberID, 1 AS line_order, Seq AS line_content
    FROM enriched_members WHERE output_dir = 'out_query_only'
  )
  
  SELECT cluster_name, line_content
  FROM fasta_lines_query_only
  ORDER BY cluster_name, MemberID, line_order
)
TO 'out_query_only' (
  FORMAT CSV,
  HEADER false,
  PARTITION_BY cluster_name,
  FILE_EXTENSION 'fasta'
);
"

## Execute the command
echo -e "..Executing DuckDB command\n"
duckdb -c "${SEQS_COMMAND}"
                             
