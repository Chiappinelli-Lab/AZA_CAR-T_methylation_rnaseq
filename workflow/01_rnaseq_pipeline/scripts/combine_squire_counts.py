#!/usr/bin/env python3
"""
Merge per-sample SQuIRE count files into DESeq2-ready matrices.

Produces up to three output files in --output_dir:
  squire_locus_counts.tsv     (from _TEcounts.txt files)
  squire_subfamily_counts.tsv (from _subFcounts.txt files)
  squire_gene_counts.tsv      (from _refGenecounts.txt files, if provided)
"""
import argparse
import os
import pandas as pd


def _sample_name(filepath, suffix):
    """Extract sample name by stripping the trailing suffix from the filename."""
    return os.path.basename(filepath)[: -len(suffix)]


# Column layout of each squire Count output type, per squire's own Call.py
# (create_TE_dict / create_subfamily_dict / create_count_dict): (has_header, id_col, count_col)
_LAYOUT = {
    "_TEcounts.txt":      dict(has_header=True,  id_col=3, count_col=15),  # TE_ID, tot_counts
    "_subFcounts.txt":    dict(has_header=True,  id_col=2, count_col=6),   # Subfamily:Family:Class, tot_counts
    "_refGenecounts.txt": dict(has_header=False, id_col=3, count_col=6),   # Gene_ID, count (no header row in this file)
}


def _read_count_file(filepath, suffix):
    """Return a Series (feature -> count) from a squire Count output file."""
    layout = _LAYOUT[suffix]
    df = pd.read_csv(
        filepath, sep="\t", comment="#",
        header=0 if layout["has_header"] else None,
    )
    series = df.iloc[:, layout["count_col"]]
    series.index = df.iloc[:, layout["id_col"]]
    # a feature can appear on multiple lines (e.g. multiple transcripts/strands); sum its counts
    return series.groupby(level=0).sum()


def _combine(file_list, suffix, output_path):
    frames = {}
    for fp in file_list:
        sample = _sample_name(fp, suffix)
        frames[sample] = _read_count_file(fp, suffix)
    matrix = pd.DataFrame(frames).fillna(0).astype(int)
    matrix.index.name = "feature"
    matrix.to_csv(output_path, sep="\t")


def main():
    parser = argparse.ArgumentParser()
    parser.add_argument("--output_dir", required=True)
    parser.add_argument("--locus_files", nargs="+", required=True,
                        help="Per-sample *_TEcounts.txt files from squire Count")
    parser.add_argument("--subfamily_files", nargs="+", required=True,
                        help="Per-sample *_subFcounts.txt files from squire Count")
    parser.add_argument("--gene_files", nargs="+", required=False,
                        help="Per-sample *_refGenecounts.txt files from squire Count")
    args = parser.parse_args()

    os.makedirs(args.output_dir, exist_ok=True)

    _combine(args.locus_files,    "_TEcounts.txt",     os.path.join(args.output_dir, "squire_locus_counts.tsv"))
    _combine(args.subfamily_files, "_subFcounts.txt",  os.path.join(args.output_dir, "squire_subfamily_counts.tsv"))
    if args.gene_files:
        _combine(args.gene_files, "_refGenecounts.txt", os.path.join(args.output_dir, "squire_gene_counts.tsv"))


if __name__ == "__main__":
    main()
