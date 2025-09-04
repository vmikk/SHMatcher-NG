## SH-matching pipeline

## Installation

```bash
# docker version                            # Docker is needed
# java -version                             # Java >= 17 required
# curl -s https://get.nextflow.io | bash    # install Nextflow
nextflow pull vmikk/SHMatcher-NG            # pull the pipeline
```

## Prepare the database

Download the SH database
```bash
wget https://s3.hpc.ut.ee/plutof-public/original/fc513ab4-beba-48f7-9cba-71968fe108b5.zip
mv fc513ab4-beba-48f7-9cba-71968fe108b5.zip sh_matching_data_0_5.zip
unzip sh_matching_data_0_5.zip -d ./sh_matching/
rm sh_matching_data_0_5.zip
```

Index the database for Minimap2-based workflow

```bash
minimap2 \
  -k 15 \
  -d sanger_refs_sh_full.mmi \
  -t 3 \
  sh_matching/data/sanger_refs_sh_full.fasta
```

Index the database for MMseqs2-based workflow

```bash
mkdir -p SH_db

mmseqs createdb \
  sh_matching/data/sanger_refs_sh_full.fasta \
  SH_db/SH_db

mmseqs createindex \
  SH_db/SH_db tmp \
  --search-type 3 \
  --create-lookup 1 --threads 8 \
  --strand 1 \
  --spaced-kmer-mode 0

rm -r tmp
```


## Run the pipeline


### Using Minimap2 for finding best matches to the SH database

```bash
nextflow run vmikk/SHMatcher-NG -r main \
  -resume \
  -profile docker \
  --method minimap \
  --input  $(pwd)/test_inp.fasta \
  --ref_db $(pwd)/sanger_refs_sh_full.mmi \
  --shdata $(pwd)/sh_matching/data
```

### Using MMseqs2 for finding best matches to the SH database

```bash
nextflow run vmikk/SHMatcher-NG -r main \
  -resume \
  -profile docker \
  --method     mmseqs \
  --input      $(pwd)/test_inp.fasta \
  --ref_db_dir $(pwd)/SH_db \
  --shdata     $(pwd)/sh_matching/data
```


## Pipeline parameters

| Parameter            | Description                                          | Workflow \* |
| -------------------- | ---------------------------------------------------- | ----------- |
| input                | Path to input FASTA file                             |             |
| ref_db_dir           | Path to MMseqs2-formatted DB directory               | MMseqs2     |
| ref_db               | Path to Minimap2 index file (mmi)                    | Minimap2    |
| shdata               | Path to SH database data directory                   |             |
| method               | Workflow to run: "mmseqs" or "minimap"               |             |
| preclust_id          | Preclustering: minimum sequence identity             |             |
| preclust_cov         | Preclustering: minimum coverage                      |             |
| preclust_kmers       | Preclustering: k-mers per sequence                   |             |
| preclust_kmerscale   | Preclustering: Scale factor for k-mers per sequence  |             |
| search_id            | Minimum sequence identity for search                 | MMseqs2     |
| search_cov           | Minimum alignment coverage                           | MMseqs2     |
| cov_mode             | Coverage calculation mode                            | MMseqs2     |
| max_seqs             | Max sequences kept per query in prefilter            | MMseqs2     |
| max_accept           | Max accepted alignments before stop                  | MMseqs2     |
| exact_kmer           | Use exact k-mer matching                             | MMseqs2     |
| spaced_kmer          | Spaced k-mer mode; must match DB index               | MMseqs2     |
| kmer_length          | k-mer length                                         |             |
| top_hits             | Top N hits per query to retain                       |             |
| minimap_n_candidates | Candidate alignments per query (-N) before filtering | Minimap2    |
| minimap_opts         | Extra Minimap2 options (eg, "-x map-ont")            | Minimap2    |

\*, there are workflow-specific parameters (depending on if MMseqs2 or Minimap2 is used for finding best matches to the SH database)
