#!/usr/bin/env python3

import re
import csv
import sys
from pathlib import Path
from collections import defaultdict

# ============================================================
# CONFIG
# ============================================================

LNC_TSV = "../../../results/lncRNA/lncRNA_overlapping_HERVs.gene_level.tsv"
GTF_FILE = "../../../data/gencode.v48.primary_assembly.annotation.gtf"
OUT_TSV = "../../../results/lncRNA/lncRNA_herv_coverage.tsv"

# ============================================================
# HELPERS
# ============================================================

def parse_attrs(attr_string):
    """
    Parse GTF attribute field into a dict.
    """
    attrs = {}
    for m in re.finditer(r'(\S+)\s+"([^"]+)"', attr_string):
        attrs[m.group(1)] = m.group(2)
    return attrs


def parse_coord(coord_str):
    """
    Parse coordinates like:
      chr15:24163598-24164700(-)
      chr19:42574126-42575513(+)
      KI270743.1:122714-122948(+)

    Returns: (chrom, start, end, strand)
    Using 1-based inclusive coordinates, matching GTF.
    """
    coord_str = coord_str.strip()
    m = re.match(r'^([^:]+):(\d+)-(\d+)\(([+-])\)$', coord_str)
    if not m:
        raise ValueError(f"Could not parse coordinate: {coord_str}")
    chrom = m.group(1)
    start = int(m.group(2))
    end = int(m.group(3))
    strand = m.group(4)
    if start > end:
        start, end = end, start
    return chrom, start, end, strand


def merge_intervals(intervals):
    """
    Merge a list of 1-based inclusive intervals [(start, end), ...].
    Returns merged list.
    """
    if not intervals:
        return []

    intervals = sorted(intervals, key=lambda x: (x[0], x[1]))
    merged = [list(intervals[0])]

    for s, e in intervals[1:]:
        last_s, last_e = merged[-1]
        if s <= last_e + 1:
            merged[-1][1] = max(last_e, e)
        else:
            merged.append([s, e])

    return [(s, e) for s, e in merged]


def interval_length(intervals):
    """
    Sum lengths of 1-based inclusive intervals.
    """
    return sum(e - s + 1 for s, e in intervals)


def intersect_two(a, b):
    """
    Intersect two 1-based inclusive intervals.
    a=(s1,e1), b=(s2,e2)
    Returns None or (s,e)
    """
    s = max(a[0], b[0])
    e = min(a[1], b[1])
    if s <= e:
        return (s, e)
    return None


def intersect_interval_lists(intervals_a, intervals_b):
    """
    Intersect two merged interval lists.
    Returns merged intersections.
    """
    out = []
    i = 0
    j = 0

    a = merge_intervals(intervals_a)
    b = merge_intervals(intervals_b)

    while i < len(a) and j < len(b):
        ia = a[i]
        ib = b[j]

        inter = intersect_two(ia, ib)
        if inter is not None:
            out.append(inter)

        if ia[1] < ib[1]:
            i += 1
        else:
            j += 1

    return merge_intervals(out)


# ============================================================
# LOAD GTF: gene spans + exon spans per gene
# ============================================================

gene_info = {}
gene_exons = defaultdict(list)

with open(GTF_FILE, "r") as fh:
    for line in fh:
        if line.startswith("#"):
            continue

        fields = line.rstrip("\n").split("\t")
        if len(fields) != 9:
            continue

        chrom, source, feature, start, end, score, strand, frame, attrs_str = fields
        start = int(start)
        end = int(end)
        attrs = parse_attrs(attrs_str)

        gene_id = attrs.get("gene_id")
        gene_name = attrs.get("gene_name", gene_id)
        gene_type = attrs.get("gene_type")

        if gene_id is None:
            continue

        if feature == "gene":
            gene_info[gene_id] = {
                "chrom": chrom,
                "start": start,
                "end": end,
                "strand": strand,
                "gene_name": gene_name,
                "gene_type": gene_type
            }

        elif feature == "exon":
            gene_exons[gene_id].append((start, end))

# merge exon intervals per gene
gene_exons_merged = {
    gene_id: merge_intervals(exons)
    for gene_id, exons in gene_exons.items()
}

# ============================================================
# PROCESS lncRNA/HERV overlap table
# ============================================================

rows_out = []

with open(LNC_TSV, "r") as fh:
    reader = csv.DictReader(fh, delimiter="\t")

    for row in reader:
        gene_id = row["gene_stable_id"]

        if gene_id not in gene_info:
            # skip genes absent from GTF
            continue

        g = gene_info[gene_id]
        chrom = g["chrom"]
        strand = g["strand"]
        gene_interval = [(g["start"], g["end"])]
        gene_len = interval_length(gene_interval)

        exon_intervals = gene_exons_merged.get(gene_id, [])
        exon_len = interval_length(exon_intervals)

        # Parse HERV coords from semicolon-separated field
        herv_coords_raw = row["overlapping_coords"].strip()
        herv_intervals_same_gene = []

        if herv_coords_raw:
            for token in herv_coords_raw.split(";"):
                token = token.strip()
                if not token:
                    continue
                h_chrom, h_start, h_end, h_strand = parse_coord(token)

                # keep only same chromosome and strand as gene
                if h_chrom != chrom:
                    continue
                if h_strand != strand:
                    continue

                herv_intervals_same_gene.append((h_start, h_end))

        herv_intervals_merged = merge_intervals(herv_intervals_same_gene)

        # overlap with gene body
        herv_in_gene = intersect_interval_lists(herv_intervals_merged, gene_interval)
        herv_in_gene_bp = interval_length(herv_in_gene)

        # overlap with exons
        herv_in_exons = intersect_interval_lists(herv_intervals_merged, exon_intervals)
        herv_in_exons_bp = interval_length(herv_in_exons)

        gene_body_fraction = herv_in_gene_bp / gene_len if gene_len > 0 else 0.0
        exonic_fraction = herv_in_exons_bp / exon_len if exon_len > 0 else 0.0

        rows_out.append({
            "gene_stable_id": gene_id,
            "gene_name": row["gene_name"],
            "gene_type": row["gene_type"],
            "chrom": chrom,
            "strand": strand,
            "gene_start": g["start"],
            "gene_end": g["end"],
            "gene_body_bp": gene_len,
            "exonic_bp": exon_len,
            "n_overlapping_hervs": row["n_overlapping_hervs"],
            "herv_bp_in_gene_body": herv_in_gene_bp,
            "herv_bp_in_exons": herv_in_exons_bp,
            "gene_body_herv_fraction": round(gene_body_fraction, 4),
            "gene_body_herv_percent": round(gene_body_fraction * 100, 2),
            "exonic_herv_fraction": round(exonic_fraction, 4),
            "exonic_herv_percent": round(exonic_fraction * 100, 2),
            "overlapping_subfamilies": row["overlapping_subfamilies"]
        })

# Sort descending by exonic percent, then gene body percent
rows_out = sorted(
    rows_out,
    key=lambda x: (x["exonic_herv_percent"], x["gene_body_herv_percent"]),
    reverse=True
)

# ============================================================
# WRITE OUTPUT
# ============================================================

fieldnames = [
    "gene_stable_id", "gene_name", "gene_type",
    "chrom", "strand", "gene_start", "gene_end",
    "gene_body_bp", "exonic_bp", "n_overlapping_hervs",
    "herv_bp_in_gene_body", "herv_bp_in_exons",
    "gene_body_herv_fraction", "gene_body_herv_percent",
    "exonic_herv_fraction", "exonic_herv_percent",
    "overlapping_subfamilies"
]

with open(OUT_TSV, "w", newline="") as outfh:
    writer = csv.DictWriter(outfh, fieldnames=fieldnames, delimiter="\t")
    writer.writeheader()
    writer.writerows(rows_out)

print(f"Wrote: {OUT_TSV}")
print(f"Genes processed: {len(rows_out)}")
