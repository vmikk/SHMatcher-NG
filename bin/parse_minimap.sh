#!/bin/bash

## Script to parse minimap2 output

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
        (pident_gc==""?"":sprintf("%.2f",pident_gc)),         # gap-compressed identity
        sprintf("%.3f", qcov),                                # query coverage
        alnlen,                                               # alignment length
        mapq,                                                 # mapping quality
        tpr                                                   # type of alignment (0 = primary, 1 = secondary and other)
    }' $1
