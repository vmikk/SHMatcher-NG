
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


// Prepare numeric keys for cluster membership and cluster representative list
process prepare_cluster_keys {

    input:
      path(q_db, stageAs: "q_db/*")  // database of query sequences (`q_db/q_db`)
      path membership                // query_cluster_membership.tsv

    output:
      path "cluster_members_numeric.tsv", emit: members_numeric
      path "cluster_reps.keys",           emit: reps_keys

    script:
    """
    echo -e "Preparing numeric keys for cluster membership\\n"

    awk 'NR==FNR{a[\$2]=\$1; next} {print a[\$1]"\\t"a[\$2]}' \
      q_db/q_db.lookup \
      "${membership}" \
      > cluster_members_numeric.tsv

    cut -f1 cluster_members_numeric.tsv \
      | sort -uh > cluster_reps.keys
    """
}


// Global search of all queries against reference DB
// (nucleotide-nucleotide, both strands)
process mmseqs_search {

    input:
      path(q_db, stageAs: "q_db/*")  // database of query sequences (`q_db/q_db`)
      path(refs, stageAs: 'SH_db/*') // reference database (`SH_db/...`)

    output:
      path "all_hits/*",    emit: hits
      path "best_hits.tsv", emit: best_hits

    script:
    """
    echo -e "..MMseqs global search\\n"

    ## Get the base name of the reference database
    REFDB=\$(find SH_db -name "*.lookup" | sed 's/\\.lookup// ; s/SH_db\\///')

    mkdir -p all_hits tmp_search

    ## Run the global search
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

    ### Select the best hit for each query
    echo -e "..Selecting the best hit for each query\\n"

    ## Keep only the top hit per query in the result DB
    mkdir -p best_hits

    mmseqs filterdb \
      all_hits/all_hits \
      best_hits/best_hits \
      --extract-lines 1 \
      --threads ${task.cpus}

    ## Convert to tabular format
    mmseqs convertalis \
      q_db/q_db \
      SH_db/\${REFDB} \
      best_hits/best_hits \
      best_hits.tsv \
      --format-output "query,target,evalue,bits,alnlen,pident,qcov,tcov" \
      --threads ${task.cpus}

    ## Clean up
    rm -r tmp_search

    """
}



// Prepare FASTA files for each (compound) cluster
process cluster_extract {

    input:
      path(q_db, stageAs: "q_db/*")     // database of query sequences (`q_db/q_db`)
      path(refs, stageAs: 'SH_db/*')    // reference database (`SH_db/...`)
      path(hits, stageAs: 'all_hits/*') // all hits (`all_hits/all_hits`)
      path members_numeric
      path reps_keys

    output:
      path "out_with_ref/*.fasta",   optional: true, emit: clusters_with_ref
      path "out_query_only/*.fasta", optional: true, emit: clusters_query_only

    script:
    """
    echo -e "Preparing FASTA files for each (compound) cluster\\n"

    mkdir -p out_with_ref out_query_only

    parallel --jobs ${task.cpus} \
      -a "${reps_keys}" \
      "cluster_extract.sh \
        --rep            {} \
        --q-db-dir       q_db \
        --ref-db-dir     SH_db \
        --all-hits-dir   all_hits \
        --members-tsv    cluster_members_numeric.tsv \
        --out-with-ref   out_with_ref \
        --out-query-only out_query_only"
    """
}

// Generate a distance matrix
process calc_distmx {

    input:
      path inp   // FASTA file

    output:
      path "mx_*.txt", emit: mx

    script:
    inp_base = inp.getBaseName()
    """
    echo -e "Calculating distance matrix for a cluster or sequences\\n"

    usearch \
      -calc_distmx ${inp} \
      -tabbedout mx_${inp_base}.txt \
      -maxdist   0.03 \
      -termdist  0.3 \
      -gapopen   "1.0I/0.0E" \
      -gapext    "1.0I/0.0E" \
      -threads   ${task.cpus}
    """
}

// Agglomerative clustering of a distance matrix
process cluster_aggd {

    input:
      path mx              // distance matrix in tabbed pairs format

    output:
      path "calc_distm_out/*", emit: clusters

    script:
    mxbase = mx.getBaseName()
    """
    echo -e "Clustering distance matrix\\n"

    # usearch -cluster_aggd ${mx} -clusterout clusters.txt -id 0.97 -linkage min

    mkdir -p calc_distm_out

    parallel -j1 \
      "goclust \
        --input  ${mx} \
        --output calc_distm_out/${mxbase}_out_{} \
        --method single \
        --cutoff {}" \
      ::: 0.03 0.025 0.02 0.015 0.01 0.005

    """
}





workflow {

  // Input FASTA file
  ch_inp = channel.fromPath(params.input)
  
  // Reference database (MMseqs-formatted)
  ch_ref = channel.fromPath("${params.ref_db_dir}/*", checkIfExists: true).collect()

  // Lin-cluster query sequences
  precluster(ch_inp)

  // Cluster membership preparation
  prepare_cluster_keys(
    precluster.out.qdb,
    precluster.out.membership)

  // Global search (queries vs SH database)
  mmseqs_search(
    precluster.out.qdb,
    ch_ref)

  // Prepare FASTA files for each (compound) cluster
  cluster_extract(
    precluster.out.qdb,
    ch_ref,
    mmseqs_search.out.hits,
    prepare_cluster_keys.out.members_numeric,
    prepare_cluster_keys.out.reps_keys)

  // Compound clusters
  ch_cls = cluster_extract.out.clusters_with_ref.flatten()

  // Generate a distance matrix (for each cluster of sequences)
  calc_distmx(ch_cls)
  
  // Agglomerative clustering (using a series of thresholds)
  cluster_aggd(calc_distmx.out.mx)

}
