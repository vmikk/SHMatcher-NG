
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




// Global search of all queries against reference DB
// (nucleotide-nucleotide, both strands)
process mmseqs_search {

    input:
      path(q_db, stageAs: "q_db/*")  // database of query sequences (`q_db/q_db`)
      path(refs, stageAs: 'SH_db/*') // reference database (`SH_db/...`)

    output:
      path "all_hits/*", emit: hits

    script:
    """
    echo -e "..MMseqs global search\\n"

    ## Get the base name of the reference database
    REFDB=\$(find SH_db -name "*.lookup" | sed 's/\\.lookup// ; s/SH_db\\///')

    mkdir -p all_hits tmp_search

    mmseqs search \
      q_db/q_db \
      SH_db/\${REFDB} \
      all_hits/all_hits \
      tmp_search \
      --threads ${task.cpus} \
      --alignment-mode 3 \
      --min-seq-id 0.9 \
      -s 5.7 \
      --cov-mode 0 -c 0.7 \
      --max-accept 50 \
      --search-type 3 \
      --strand 2

    ## Clean up
    rm -r tmp_search
    """
}

    """
}



workflow {

  // Input FASTA file
  ch_inp = channel.fromPath(params.input)
  
  // Reference database (MMseqs-formatted)
  ch_ref = channel.fromPath("${params.ref_db_dir}/*", checkIfExists: true).collect()

  // Lin-cluster query sequences
  precluster(ch_inp)

  // Global search (queries vs SH database)
  mmseqs_search(
    precluster.out.qdb,
    ch_ref)
    
}
