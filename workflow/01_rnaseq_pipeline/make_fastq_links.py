#!/usr/bin/env python3
"""Create the <sample>_R{1,2}.fastq.gz links the Snakefile expects, from the FASTQs as delivered.

Delivered names are <N>:<external_event_id>.<LIB>.<LIB>_<flowcell>_S<n>_L006_R{1,2}_001.fastq.gz; the
sample manifest (metadata/haydar_rnaseq_sample_manifest.tsv) lists each file without the "<N>:" prefix,
with its sample (external_event_id) and read (read_pair_number). Sample names use "_" for "-"
(NS-07-0 -> NS_07_0), as in the count matrix and the DESeq2 scripts.

usage: make_fastq_links.py <fastq_dir> <manifest.tsv> <out_dir>
"""
import csv
import os
import sys
from pathlib import Path

fastq_dir, manifest, out_dir = Path(sys.argv[1]), Path(sys.argv[2]), Path(sys.argv[3])
out_dir.mkdir(parents=True, exist_ok=True)
by_name = {}
for p in fastq_dir.glob("*.fastq.gz"):
    by_name[p.name.split(":", 1)[-1]] = p

n = 0
with open(manifest) as fh:
    for row in csv.DictReader(fh, delimiter="\t"):
        src = by_name.get(row["file_name"])
        if src is None:
            sys.exit(f"FASTQ not found in {fastq_dir}: {row['file_name']}")
        link = out_dir / f"{row['external_event_id'].replace('-', '_')}_{row['read_pair_number']}.fastq.gz"
        if link.is_symlink() or link.exists():
            if os.path.realpath(link) == str(src.resolve()):
                n += 1
                continue
            sys.exit(f"{link} exists and points elsewhere")
        link.symlink_to(src.resolve())
        n += 1
if n != 2 * len({p.name.rsplit("_R", 1)[0] for p in out_dir.glob("*_R1.fastq.gz")}):
    sys.exit("R1/R2 links don't pair up")
print(f"{n} FASTQ links in {out_dir}")
