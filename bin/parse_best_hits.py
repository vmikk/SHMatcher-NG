#!/usr/bin/env python
"""
Shared module for parsing best hits files from different alignment tools.

Supports:
- MMseqs2 convertalis output (query,target,evalue,bits,alnlen,pident,qcov,tcov)
- minimap2 output (qseqid,sseqid,pident,pident_gc,qcov,alnlen,mapq,tpr)

Usage:
    from parse_best_hits import parse_best_hits
    
    # For MMseqs2 results
    q2best, q2pident = parse_best_hits("best_hits.tsv", "mmseqs")
    
    # For minimap2 results  
    q2best, q2pident = parse_best_hits("best_hits.tsv", "minimap")
"""

import csv


def parse_best_hits_mmseqs(tsv_path: str):
    """Parse MMseqs2 convertalis output to find best hits per query."""
    expected_cols = 8
    with open(tsv_path) as f:
        reader = csv.reader(f, delimiter="\t")
        first = next(reader, None)
        if first is None:
            return {}, {}
        
        # Check for header and column count
        columns = [c.strip().lower() for c in first]
        has_header = "query" in columns and "target" in columns
        
        if len(first) != expected_cols:
            raise ValueError(f"MMseqs2 best hits file has {len(first)} columns, expected {expected_cols} "
                           f"(query,target,evalue,bits,alnlen,pident,qcov,tcov)")
        
        if has_header:
            col_idx = {name: i for i, name in enumerate(columns)}
        else:
            col_idx = {"query": 0, "target": 1, "evalue": 2, "bits": 3, "alnlen": 4, "pident": 5, "qcov": 6, "tcov": 7}
        
        best = {}
        if not has_header and first:
            rows_iter = [first]
        else:
            rows_iter = []
        
        for row in rows_iter + list(reader):
            if len(row) != expected_cols:
                continue
            try:
                qid = row[col_idx["query"]]
                tid = row[col_idx["target"]]
                evalue = float(row[col_idx.get("evalue", 2)])
                bits = float(row[col_idx.get("bits", 3)])
                pident = float(row[col_idx.get("pident", 5)])
            except (ValueError, IndexError):
                continue
            
            prev = best.get(qid)
            key = (evalue, -bits, -pident)  # Lower evalue, higher bits, higher pident
            if prev is None or key < prev[0]:
                best[qid] = (key, tid, pident)
        
        q2best_ref = {q: v[1] for q, v in best.items()}
        q2pident = {q: v[2] for q, v in best.items()}
        return q2best_ref, q2pident


def parse_best_hits_minimap(tsv_path: str):
    """Parse minimap2 output to find best hits per query."""
    expected_cols = 8
    with open(tsv_path) as f:
        reader = csv.reader(f, delimiter="\t")
        first = next(reader, None)
        if first is None:
            return {}, {}
        
        # Check for header and column count
        columns = [c.strip().lower() for c in first]
        has_header = "qseqid" in columns and "sseqid" in columns
        
        if len(first) != expected_cols:
            raise ValueError(f"Minimap2 best hits file has {len(first)} columns, expected {expected_cols} "
                           f"(qseqid,sseqid,pident,pident_gc,qcov,alnlen,mapq,tpr)")
        
        if has_header:
            col_idx = {name: i for i, name in enumerate(columns)}
        else:
            col_idx = {"qseqid": 0, "sseqid": 1, "pident": 2, "pident_gc": 3, "qcov": 4, "alnlen": 5, "mapq": 6, "tpr": 7}
        
        best = {}
        if not has_header and first:
            rows_iter = [first]
        else:
            rows_iter = []
        
        for row in rows_iter + list(reader):
            if len(row) != expected_cols:
                continue
            try:
                qid = row[col_idx["qseqid"]]
                tid = row[col_idx["sseqid"]]
                pident = float(row[col_idx.get("pident", 2)])
                mapq = float(row[col_idx.get("mapq", 6)])
                qcov = float(row[col_idx.get("qcov", 4)])
            except (ValueError, IndexError):
                continue
            
            prev = best.get(qid)
            key = (-pident, -mapq, -qcov)  # Higher pident, higher mapq, higher qcov
            if prev is None or key < prev[0]:
                best[qid] = (key, tid, pident)
        
        q2best_ref = {q: v[1] for q, v in best.items()}
        q2pident = {q: v[2] for q, v in best.items()}
        return q2best_ref, q2pident


def parse_best_hits(tsv_path: str, method: str):
    """Parse best hits file using the specified method."""
    if method == "mmseqs":
        return parse_best_hits_mmseqs(tsv_path)
    elif method == "minimap":
        return parse_best_hits_minimap(tsv_path)
    else:
        raise ValueError(f"Unknown method '{method}'. Supported methods: 'mmseqs', 'minimap'")
