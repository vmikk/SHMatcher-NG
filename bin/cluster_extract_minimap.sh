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


