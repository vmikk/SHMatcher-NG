
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
      -c ${params.preclust_cov} \
      --cov-mode ${params.cov_mode} \
      -k ${params.kmer_length} \
      --kmer-per-seq ${params.preclust_kmers} \
      --kmer-per-seq-scale ${params.preclust_kmerscale} \
      --spaced-kmer-mode ${params.spaced_kmer} \
      --mask 0 \
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
    exact_kmer = params.exact_kmer ? "--exact-kmer-matching" : ""
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
      --search-type 3 \
      --min-seq-id ${params.search_id} \
      -c ${params.search_cov} \
      --cov-mode ${params.cov_mode} \
      --mask 0 \
      --max-accept ${params.max_accept} \
      --max-seqs   ${params.max_seqs} \
      --strand 2 \
      ${exact_kmer} \
      --spaced-kmer-mode ${params.spaced_kmer} \
      -k ${params.kmer_length}

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
      path "cluster_membership.txt", emit: ids

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
    
    ## Pool per-cluster member IDs  (table: ClusterID, MemberType (Query/Ref), MemberID
    find out_with_ref   -name "*.ids.txt" | parallel -j1 "cat {}" >  cluster_membership.txt
    find out_query_only -name "*.ids.txt" | parallel -j1 "cat {}" >> cluster_membership.txt

    ## Clean up
    find out_with_ref   -name "*.ids.txt" | parallel -j1 -X "rm {}"
    find out_query_only -name "*.ids.txt" | parallel -j1 -X "rm {}"

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
     tuple( 
      path(mx),                 // distance matrix in tabbed pairs format
      path(cluster_membership), // cluster membership TSV (for queries and refs)
      path(best_hits),          // best hits TSV (for all queries)
      path(db_centroid2sh) )    // centroid2sh mapping (SH database part)

    output:
      path "calc_distm_out/*", emit: clusters
      path "*_matches.tsv",    emit: matches

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
      ::: 0.030 0.025 0.020 0.015 0.010 0.005

    echo -e "Processing single-linkage clustering outputs\\n"

    postprocess_single_linkage.py \
      --clusters-dir       calc_distm_out \
      --cluster-membership ${cluster_membership} \
      --best-hits          ${best_hits} \
      --centroid2sh        ${db_centroid2sh} \
      --output             ${mxbase}_matches.tsv

    """
}



// Aggregate matches
//   This process should reproduce the logic 
//   of `analyse_usearch_output_sh.py` + `merge_matches.py` 
//   without nested clustering
process aggregate_matches {

    publishDir 'Results_aggregated', mode: 'copy', overwrite: true

    input:
      path(matches, stageAs: "matches1/*")  // cluster matches (aggregated at different thresholds)
      // path(...., stageAs: "matches2/*")     // ...)
      path(shs_out)                         // SH database file `shs_out.txt`
      path(sh2compound)                     // SH database file `sh2compound_mapping.txt`
      path(compounds_out)                   // SH database file `compounds_out.txt`
      path(best_hits)                       // best hits TSV (for all queries)

    output:
      path "matches_out_all.csv", emit: all

    script:
    """
    echo -e "Aggregating matches\\n"

    aggregate_matches.py \
      --matches-dirs   matches1 \
      --best-hits      ${best_hits} \
      --shs-file       ${shs_out} \
      --sh2compound    ${sh2compound} \
      --compounds-file ${compounds_out} \
      --out            matches_out_all.csv

    """
}


    """
}





workflow {

  // SH database files
  // ch_shd = Channel.value(params.shdata)   // `sh_matching/data` path

  db_centroid2sh    = Channel.fromPath( params.shdata + '/centroid2sh_mappings.txt').first()
  db_compound2seq   = Channel.fromPath( params.shdata + '/compound2seq_mapping.txt').first()
  db_plutof_tax     = Channel.fromPath( params.shdata + '/plutof_taxonomy.txt').first()
  db_compound       = Channel.fromPath( params.shdata + '/compounds_out.txt').first() 
  db_sh2compound    = Channel.fromPath( params.shdata + '/sh2compound_mapping.txt').first()
  db_shs            = Channel.fromPath( params.shdata + '/shs_out.txt').first()
  db_sanger_sh      = Channel.fromPath( params.shdata + '/sanger_refs_sh.fasta').first()
  db_sanger_sh_full = Channel.fromPath( params.shdata + '/sanger_refs_sh_full.fasta').first()
  // TODO - check which of these are needed


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
  


  // We need to reuse shared channels in `cluster_aggd`
  // -> make a single tuple
  ch_meta = cluster_extract.out.ids
           .combine(mmseqs_search.out.best_hits)
           .combine(db_centroid2sh)
           .map { ids, best_hits, centroid2sh -> tuple(ids, best_hits, centroid2sh) }

  ch_mx = calc_distmx.out.mx
    .combine(ch_meta)
    .map { mx, ids, best_hits, centroid2sh -> tuple(
      mx,               // distance matrix for a cluster
      ids,              // cluster membership (for all queries and refs)
      best_hits,        // best hits for all queries
      centroid2sh       // centroid2sh mapping (SH database part)
    ) }

  // Agglomerative clustering (using a series of thresholds)
  cluster_aggd(ch_mx)


}
