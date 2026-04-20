#!/usr/bin/env python3
"""
Parse FIMO TSV output and generate:
1) A TSV sorted by genomic coordinates
2) A BED with genomic coordinates of each motif hit

IMPORTANT:
- sequence_name encodes genomic interval and genomic strand:
    <subfamily>-int_pos_<chrom>_<start>_<end>_strand_<+|->

- FIMO start/stop are 1-based inclusive relative to the *scanned FASTA sequence*.
- For antisense runs, the FASTA sequence is reverse-complemented compared to the sense FASTA,
  but the header strand remains the genomic strand. Therefore you MUST specify --run-orientation.

Usage:
  python 05_parse_fimo_internal_rbp_to_bed.py \
    --input fimo.tsv \
    --output_tsv fimo_sorted.tsv \
    --output_bed fimo_hits.bed \
    --run-orientation sense
"""

import pandas as pd
import argparse
import re
from tqdm import tqdm

tqdm.pandas()


def flip_strand(s: str) -> str:
    if s == '+':
        return '-'
    if s == '-':
        return '+'
    return s


def parse_coordinates_from_name(name: str):
    """
    Extract chrom, start, end, strand from names like:
      HERVK14-int_pos_GL000008.2_57049_58824_strand_+
      HERVK13-int_pos_chr1_3801949_3803126_strand_-_divergence_6.3_1
    """
    pattern = r'_(?P<chrom>(chr)?[^_]+)_(?P<start>\d+)_(?P<end>\d+)_strand_(?P<strand>[+-])'
    m = re.search(pattern, name)
    if not m:
        raise ValueError(f"Invalid sequence name format: {name}")
    chrom = m.group('chrom')
    start = int(m.group('start'))
    end = int(m.group('end'))
    strand = m.group('strand')
    return chrom, start, end, strand


def main(input_file, output_tsv, output_bed, run_orientation: str):
    print(f"📥 Reading: {input_file}")
    df = pd.read_csv(input_file, sep='\t', comment='#')

    if df.empty:
        raise SystemExit("No rows found in input TSV (is it empty or only comments?).")

    print(f"🔍 Parsing coordinates from {len(df)} motif hits...")
    coord_data = [parse_coordinates_from_name(name) for name in tqdm(df['sequence_name'], desc="Parsing coords")]
    df[['chrom', 'internal_start', 'internal_end', 'internal_strand']] = pd.DataFrame(coord_data, index=df.index)

    # Decide if the scanned sequence is reverse-complemented relative to genomic coordinates.
    # If your sense FASTA is transcript-oriented (5'->3' of genomic strand),
    # then antisense FASTA flips that orientation for all loci.
    seq_is_rc = (run_orientation == "antisense")
    gen_is_minus = (df['internal_strand'] == '-')

    # XOR: reverse mapping needed when genomic strand is '-' (transcript-oriented FASTA),
    # but flips again if we're scanning the antisense FASTA
    needs_reverse = gen_is_minus ^ seq_is_rc

    print("🧬 Generating genomic coordinates for BED...")
    # FIMO start/stop are 1-based inclusive.
    # We output BED: 0-based start, end-exclusive.

    df['genomic_start'] = pd.NA
    df['genomic_end'] = pd.NA

    # Forward mapping: genomic_start + (start-1), genomic_start + stop
    fwd = ~needs_reverse
    df.loc[fwd, 'genomic_start'] = df.loc[fwd, 'internal_start'] + (df.loc[fwd, 'start'] - 1)
    df.loc[fwd, 'genomic_end']   = df.loc[fwd, 'internal_start'] + df.loc[fwd, 'stop']

    # Reverse mapping: use internal_end
    rev = needs_reverse
    df.loc[rev, 'genomic_start'] = df.loc[rev, 'internal_end'] - df.loc[rev, 'stop']
    df.loc[rev, 'genomic_end']   = df.loc[rev, 'internal_end'] - (df.loc[rev, 'start'] - 1)

    print("🧪 Sanity-checking coordinates...")
    bad = df.progress_apply(lambda r: r["genomic_end"] <= r["genomic_start"], axis=1)
    n_bad = int(bad.sum())
    if n_bad > 0:
        raise ValueError(f"Found {n_bad} hits with genomic_end <= genomic_start. Check coordinate conventions.")

    # Motif hit strand: FIMO reports strand relative to scanned sequence.
    # If scanned sequence was reverse relative to genome, flip it.
    df['match_strand'] = df['strand']
    df.loc[rev, 'match_strand'] = df.loc[rev, 'strand'].map(flip_strand)

    print("📊 Sorting TSV by chromosome, genomic_start...")
    df_sorted = df.sort_values(by=['chrom', 'genomic_start', 'genomic_end'])

    print(f"💾 Saving sorted TSV to {output_tsv}")

    tsv_cols = [
        "motif_id", "motif_alt_id", "sequence_name",
        "start", "stop", "strand",
        "score", "p-value", "q-value", "matched_sequence",
        "chrom", "internal_start", "internal_end", "internal_strand", "match_strand",
    ]

    # keep only columns that exist (sometimes q-value is missing depending on FIMO options)
    tsv_cols = [c for c in tsv_cols if c in df_sorted.columns]

    df_sorted.loc[:, tsv_cols].to_csv(output_tsv, sep="\t", index=False)


    print("🏷️ Building BED names...")
    df_sorted["bed_name"] = df_sorted.progress_apply(
        lambda r: f"{r['sequence_name']}_{r['motif_id']}_{r['motif_alt_id']}",
        axis=1
    )

    print(f"💾 Saving BED file to {output_bed}")
    bed_df = df_sorted[["chrom", "genomic_start", "genomic_end", "bed_name", "score", "match_strand"]]
    bed_df.sort_values(by=["chrom", "genomic_start"]).to_csv(output_bed, sep="\t", header=False, index=False)

    print("✅ Done!")


if __name__ == '__main__':
    parser = argparse.ArgumentParser(description='Parse FIMO TSV and produce sorted TSV and BED.')
    parser.add_argument('--input', required=True, help='Input FIMO TSV file')
    parser.add_argument('--output_tsv', required=True, help='Output sorted TSV file')
    parser.add_argument('--output_bed', required=True, help='Output BED file')
    parser.add_argument('--run-orientation', choices=['sense', 'antisense'], required=True,
                        help='Which FASTA orientation was scanned to produce this FIMO TSV.')
    args = parser.parse_args()
    main(args.input, args.output_tsv, args.output_bed, args.run_orientation)
