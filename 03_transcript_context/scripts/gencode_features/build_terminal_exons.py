#!/usr/bin/env python3
"""
03_transcript_context step 4 (effective order 11): transcript-aware 5'/3' UTR
split and terminal (last) exon per transcript.

Reimplements, rule for rule, the logic embedded in
scripts/context/07_overlap_ltrs_vs_features.sh (the UTR split is also verbatim
in 06_overlap_domains_vs_features.sh), which produced the historical
context/assets/gencode.v48.{utr_5p3p,utr3,last_exon_per_tx}.sorted.bed.gz.
Everything is transcript-resolved; nothing is collapsed to gene level.

Inputs (accepted step-3 sorted step-1 assets, BED 0-based half-open):
  exons  transcript_id|gene_id|gene_name|exon|gene_type|transcript_type
  CDS    transcript_id|gene_id|gene_name|CDS|gene_type|transcript_type
  UTR    transcript_id|gene_id|gene_name|UTR|gene_type|transcript_type
         (GENCODE v48 emits only the generic "UTR" feature)

CDS span, per transcript: cds_min = min(CDS start), cds_max = max(CDS end).

UTR label, for each generic UTR segment of transcript t (interval [s, e)):
  t has no CDS                      -> UTR_nocds
  strand +:  e <= cds_min           -> five_prime_UTR
             s >= cds_max           -> three_prime_UTR
             otherwise              -> UTR_internal
  strand -:  e <= cds_min           -> three_prime_UTR
             s >= cds_max           -> five_prime_UTR
             otherwise              -> UTR_internal
The label replaces the first "|UTR|" token of the name. Every UTR segment is
emitted once, coordinates unchanged; nothing is split or dropped.

Terminal exon, per transcript, scanning the genome-sorted exon file in order:
  strand +: the exon with the greatest end   (replace only if end > current)
  strand -: the exon with the smallest start (replace only if start < current)
A single-exon transcript's only exon is its terminal exon. The first "|exon|"
token of the name becomes "|last_exon|"; coordinates are the source exon's.
No first-exon output exists historically, so none is produced.

Outputs (features/sorted/, genome order via bedtools sort -g, as historically):
  gencode.v48.utr_5p3p.sorted.bed.gz        all UTR segments, labelled
  gencode.v48.utr3.sorted.bed.gz            the three_prime_UTR subset
  gencode.v48.last_exon_per_tx.sorted.bed.gz one row per transcript
Ordering reproduces history: rows are first sorted by (chrom, start, end,
strand, name) -- lexical chrom, and for last exons the name BEFORE relabelling
-- then passed through bedtools sort -g, whose tie order depends on that input
order. gzip -n makes the files byte-reproducible.

Neither the gene-level step-2 assets (D022) nor any TSS file (D023) is read.
"""
from __future__ import annotations

import argparse
import gzip
import os
import re
import subprocess
import sys
import tempfile

UTR3_RE = re.compile(r"\|three_prime_UTR\|")


def die(msg):
    sys.exit(f"[error] {msg}")


def read_rows(path):
    with gzip.open(path, "rt") as fh:
        for line in fh:
            if not line.strip():
                continue
            f = line.rstrip("\n").split("\t")[:6]
            if len(f) != 6:
                die(f"{path}: fewer than 6 columns")
            yield f[0], int(f[1]), int(f[2]), f[3], f[4], f[5]


def utr_split(cds_path, utr_path):
    bounds = {}
    for chrom, s, e, name, score, strand in read_rows(cds_path):
        tx = name.split("|")[0]
        if tx not in bounds:
            bounds[tx] = [chrom, strand, s, e]
        else:
            bounds[tx][2] = min(bounds[tx][2], s)
            bounds[tx][3] = max(bounds[tx][3], e)
    out = []
    for chrom, s, e, name, score, strand in read_rows(utr_path):
        tx = name.split("|")[0]
        if tx not in bounds:
            label = "UTR_nocds"
        else:
            cmin, cmax = bounds[tx][2], bounds[tx][3]
            if strand == "+":
                label = "five_prime_UTR" if e <= cmin else "three_prime_UTR" if s >= cmax else "UTR_internal"
            else:
                label = "three_prime_UTR" if e <= cmin else "five_prime_UTR" if s >= cmax else "UTR_internal"
        new = name.replace("|UTR|", f"|{label}|", 1) if "|UTR|" in name else name + f"|{label}"
        out.append((chrom, s, e, new, score, strand))
    out.sort(key=lambda x: (x[0], x[1], x[2], x[5], x[3]))
    return out


def last_exons(exons_path):
    best = {}
    for row in read_rows(exons_path):
        chrom, s, e, name, score, strand = row
        tx = name.split("|")[0]
        if tx not in best:
            best[tx] = row
            continue
        _, s0, e0, _, _, _ = best[tx]
        if strand == "+":
            if e > e0:
                best[tx] = row
        elif s < s0:
            best[tx] = row
    rows = sorted(best.values(), key=lambda x: (x[0], x[1], x[2], x[5], x[3]))
    return [(c, s, e, n.replace("|exon|", "|last_exon|", 1) if "|exon|" in n else n + "|last_exon", sc, st)
            for c, s, e, n, sc, st in rows]


def bed_text(rows):
    return "".join(f"{c}\t{s}\t{e}\t{n}\t{sc}\t{st}\n" for c, s, e, n, sc, st in rows)


def genome_sort_gz(bedtools, genome, text, out_path):
    """bedtools sort -g on `text`, gzip -n to out_path; returns sorted text."""
    srt = subprocess.run([bedtools, "sort", "-g", genome, "-i", "stdin"], input=text.encode(),
                         stdout=subprocess.PIPE, stderr=subprocess.PIPE, check=True).stdout
    with open(out_path, "wb") as fh:
        subprocess.run(["gzip", "-n", "-6"], input=srt, stdout=fh, check=True)
    os.chmod(out_path, 0o644)
    return srt.decode()


def main():
    ap = argparse.ArgumentParser(description="transcript-aware UTR split and terminal exons")
    ap.add_argument("--exons", required=True)
    ap.add_argument("--cds", required=True)
    ap.add_argument("--utr", required=True)
    ap.add_argument("--genome", required=True)
    ap.add_argument("--outdir", required=True)
    ap.add_argument("--bedtools", default=os.environ.get("NCRNA_BEDTOOLS", "bedtools"))
    a = ap.parse_args()

    for p in (a.exons, a.cds, a.utr):
        base = os.path.basename(p)
        if "gene_exons_merged" in base or "gene_introns" in base or "tss" in base.lower():
            die(f"refusing {base}: step 4 is transcript-resolved and must not read D022/D023 assets")
        if "historical_reproduction" in p:
            die("refusing a historical-reproduction input")
    ver = subprocess.run([a.bedtools, "--version"], stdout=subprocess.PIPE, text=True).stdout.strip()
    if ver != "bedtools v2.31.1":
        die(f"bedtools version {ver!r}, release pins v2.31.1")

    split = utr_split(a.cds, a.utr)
    last = last_exons(a.exons)

    os.makedirs(a.outdir, exist_ok=True)
    stage = tempfile.mkdtemp(prefix=".stage_terminal.", dir=a.outdir)
    names = {"split": "gencode.v48.utr_5p3p.sorted.bed.gz",
             "utr3": "gencode.v48.utr3.sorted.bed.gz",
             "last": "gencode.v48.last_exon_per_tx.sorted.bed.gz"}
    try:
        split_sorted = genome_sort_gz(a.bedtools, a.genome, bed_text(split), os.path.join(stage, names["split"]))
        utr3_text = "".join(l + "\n" for l in split_sorted.split("\n") if l and UTR3_RE.search(l.split("\t")[3]))
        genome_sort_gz(a.bedtools, a.genome, utr3_text, os.path.join(stage, names["utr3"]))
        genome_sort_gz(a.bedtools, a.genome, bed_text(last), os.path.join(stage, names["last"]))
        for n in names.values():
            if os.path.getsize(os.path.join(stage, n)) == 0:
                die(f"empty output {n}")
        for n in names.values():
            os.replace(os.path.join(stage, n), os.path.join(a.outdir, n))
    finally:
        for f in os.listdir(stage):
            os.remove(os.path.join(stage, f))
        os.rmdir(stage)
    print(f"[ok] {ver}")
    print(f"[ok] utr_5p3p rows {len(split)}; utr3 rows {utr3_text.count(chr(10))}; last_exon rows {len(last)}")


if __name__ == "__main__":
    main()
