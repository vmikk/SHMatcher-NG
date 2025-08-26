
/*
params.input = "input.fasta"
params.preclust_id   = 0.9
params.preclust_cov  = 0.7
params.preclust_kmer = 80
params.ref_db_dir    = "ref_db"
*/


// Fast pre-clustering of the dataset (to split into chunks prior processing)
process precluster {

    input:
      path input

    output:
      path "q_db/*",     emit: qdb
      path "query_cluster_membership.tsv", emit: membership

    script:
    """
    echo -e "..DB creation\\n"

    mkdir -p q_db

    mmseqs createdb \
      --dbtype 2 \
      --createdb-mode 0 \
      --shuffle 0 \
      ${input} \
      q_db/q_db

    echo -e "..Lin-Clustering\\n"

    mmseqs linclust \
      q_db/q_db \
      linclusters_db \
      tmplc \
      --min-seq-id ${params.preclust_id} \
      --cluster-mode 0 \
      --similarity-type 2 \
      --cov-mode 0 \
      -c ${params.preclust_cov} \
      --kmer-per-seq ${params.preclust_kmer} \
      --split-memory-limit 100G \
      --remove-tmp-files 1 \
      --threads ${task.cpus}

    echo -e "..Generating TSV-formatted output of clustering\\n"

    mmseqs createtsv \
      q_db/q_db \
      q_db/q_db \
      linclusters_db \
      query_cluster_membership.tsv \
      --threads ${task.cpus}

    ## Clean up
    rm -r tmplc
    # rm linclusters_db*

    """
}


// Create compund clusters by mapping the database to input sequences
process create_compound_clusters_vsearch {

    cpus = 8

    input:
      path input

    output:
      path "closedref.80.map.uc", emit: uc

    script:
    """

    vsearch \
    --usearch_global ${input} \
    --db "$udb_data_dir/sanger_refs_sh_full.udb" \
    --strand plus \
    --id 0.8 \
    --threads ${task.cpus} \
    --iddef 2 \
    --gapopen 1I/0E \
    --gapext  1I/0E \
    --uc "closedref.80.map.uc" \
    --maxaccepts 3 \
    --maxrejects 0

    """
}



workflow {

  ch_inp = channel.fromPath(params.input)

  precluster(ch_inp)


    
}
