#!/bin/bash

## Script to parse minimap2 output

awk -F "\t" -v OFS='\t' '
      BEGIN{ print "qseqid","sseqid","pident","pident_gc","qcov","mapq","tp" }
      {
      pident = ($10>0 && $11>0) ? ($10/$11)*100 : "";         # percent identity
      qcov   = ($2>0) ? (($4-$3)/$2)*100 : "";                # query coverage
      for (i=13; i<=NF; i++) {                                # since optional fields are not positionally fixed
        if ($i ~ /^tp:A:/)      tp = substr($i,6,1)           #   from "tp:A:P|S|I"  grab P/S/I
        else if ($i ~ /^de:f:/) de = substr($i,6)             #   from "de:f:0.0123" grab 0.0123
      }
      pident_gc = (de!="") ? (100*(1-de)) : "";               # gap-compressed identity
      print $1,                                               # query name
        $6,                                                   # target name
        sprintf("%.3f", pident),                              # percent identity
        (pident_gc==""?"":sprintf("%.2f",pident_gc)),         # gap-compressed identity
        sprintf("%.3f", qcov),                                # query coverage
        $12,                                                  # mapping quality
        tp                                                    # type of alignment (P = primary, S = secondary)
    }' $1
