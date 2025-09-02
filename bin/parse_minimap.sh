#!/bin/bash


# Parameter parsing
INPUT_FILE=""
MAX_HITS="10"
OUTPUT_FILE=""
USE_STDIN=false

# Parse arguments
while [[ $# -gt 0 ]]; do
    case $1 in
        -o|--output)
            OUTPUT_FILE="$2"
            shift 2
            ;;
        -h|--help)
            echo "Usage: $0 [input_file|-] [max_hits_per_query] [-o|--output output_file]" >&2
            echo "" >&2
            echo "Arguments:" >&2
            echo "  input_file:         input file in minimap2 PAF format (use '-' or omit for stdin)" >&2
            echo "  max_hits_per_query: maximum number of hits to return per query (default: 10)" >&2
            echo "  -o, --output:       output TSV file (default: stdout)" >&2
            echo "" >&2
            exit 0
            ;;
        -)
            INPUT_FILE="-"
            shift
            ;;
        -*)
            echo "Error: Unknown option $1" >&2
            exit 1
            ;;
        *)
            if [[ -z "$INPUT_FILE" ]]; then
                INPUT_FILE="$1"
            elif [[ "$MAX_HITS" == "10" ]]; then
                MAX_HITS="$1"
            else
                echo "Error: Too many positional arguments" >&2
                exit 1
            fi
            shift
            ;;
    esac
done

## Handle input source
if [[ -z "$INPUT_FILE" || "$INPUT_FILE" == "-" ]]; then
    USE_STDIN=true
    INPUT_SOURCE="/dev/stdin"
else
    INPUT_SOURCE="$INPUT_FILE"
    
    ## Validate input file exists
    if [ ! -f "$INPUT_FILE" ]; then
        echo "Error: Input file '$INPUT_FILE' does not exist" >&2
        exit 1
    fi
fi

## Validate max_hits is a positive integer
if ! [[ "$MAX_HITS" =~ ^[1-9][0-9]*$ ]]; then
    echo "Error: max_hits_per_query must be a positive integer, got: '$MAX_HITS'" >&2
    exit 1
fi


## Build DuckDB query
DUCKDB_QUERY="
WITH ranked_data AS (
  SELECT *,
         ROW_NUMBER() OVER (
           PARTITION BY qseqid 
           ORDER BY 
             tpr ASC,                                                     -- Primary alignment before secondary
             CASE 
               WHEN pident_gc = '' OR pident_gc IS NULL THEN pident
               ELSE CAST(pident_gc AS DOUBLE) 
             END DESC,                                                    -- Gap-compressed identity descending, fallback to pident
             qcov DESC,                                                   -- Query coverage descending  
             alnlen DESC,                                                 -- Alignment length descending
             mapq ASC                                                     -- Mapping quality ascending
         ) as rn
  FROM read_csv('/dev/stdin', header=true, sep='\t', 
                columns={
                  'qseqid': 'VARCHAR',
                  'sseqid': 'VARCHAR', 
                  'pident': 'DOUBLE',
                  'pident_gc': 'VARCHAR',                                 -- Can be empty string, handle in CASE
                  'qcov': 'DOUBLE',
                  'alnlen': 'INTEGER', 
                  'mapq': 'INTEGER',
                  'tpr': 'INTEGER'
                })
)
SELECT qseqid, sseqid, pident, pident_gc, qcov, alnlen, mapq, tpr
FROM ranked_data 
WHERE rn <= ${MAX_HITS}
ORDER BY qseqid, rn"

## Build DuckDB command based on output destination
if [[ -n "$OUTPUT_FILE" ]]; then
    ## Output to file
    DUCKDB_CMD="duckdb -c \"COPY (${DUCKDB_QUERY}) TO '${OUTPUT_FILE}' (FORMAT CSV, DELIMITER '\t', HEADER);\""
else
    ## Output to stdout
    DUCKDB_CMD="duckdb -c \"${DUCKDB_QUERY};\""
fi


## Execute parsing and sorting workflow
awk -F "\t" -v OFS='\t' '
      BEGIN{ print "qseqid","sseqid","pident","pident_gc","qcov","alnlen","mapq","tpr" }
      {
      # Reset per-record optional tags
      tp=""; de="";

      pident = ($10>0 && $11>0) ? ($10/$11)*100 : 0;          # percent identity
      qcov   = ($2>0) ? (($4-$3)/$2)*100 : 0;                 # query coverage
      alnlen = ($11>0 ? $11 : 0);                             # alignment block length
      
      for (i=13; i<=NF; i++) {                                # since optional fields are not positionally fixed
        if ($i ~ /^tp:A:/)      tp = substr($i,6,1)           #   from "tp:A:P|S|I"  grab P/S/I
        else if ($i ~ /^de:f:/) de = substr($i,6)             #   from "de:f:0.0123" grab 0.0123
      }
      pident_gc = (de!="") ? (100*(1-de)) : 0;                # gap-compressed identity
      tpr = (tp=="P" ? 0 : 1);                                # type of alignment (0 = primary, 1 = secondary)

      # MAPQ: treat 255 (missing) as 0; secondaries have MAPQ=0 anyway
      mapq_raw = ($12+0);
      mapq     = (mapq_raw==255 ? 0 : mapq_raw);

      print $1,                                               # query name
        $6,                                                   # target name
        sprintf("%.3f", pident),                              # percent identity
        (pident_gc==""?"":sprintf("%.3f",pident_gc)),         # gap-compressed identity
        sprintf("%.3f", qcov),                                # query coverage
        alnlen,                                               # alignment length
        mapq,                                                 # mapping quality
        tpr                                                   # type of alignment (0 = primary, 1 = secondary and other)
    }' "$INPUT_SOURCE" \
    | eval "$DUCKDB_CMD"
