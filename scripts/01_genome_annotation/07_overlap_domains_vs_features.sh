#!/usr/bin/env bash
set -euo pipefail

# ----------------
# INPUTS
# ----------------
GENOME="input/genome/grch38.primary.genome"

DOMAINS_BED="input/annotations/ERV_GyDB_v6_domains.bed"

# Pre-built, transcript-aware, sorted features (your existing files)
BED_DIR="bed"
TX_EXONS_BED="${BED_DIR}/gencode.v48.exons.sorted.bed.gz"
CDS_BED="${BED_DIR}/gencode.v48.cds.sorted.bed.gz"
UTR_BED="${BED_DIR}/gencode.v48.utr.sorted.bed.gz"

# Optional gene-merged exons (already exists in your bed dir)
GENE_EXONS_MERGED_BED="${BED_DIR}/gencode.v48.gene_exons_merged.bed.gz"

# ----------------
# OUTPUTS
# ----------------
OUTDIR="context/results/overlap"
LOGDIR="context/logs"
mkdir -p "$OUTDIR" "$LOGDIR" "context/assets"

LOG="${LOGDIR}/05_overlap_domains_vs_gencode_transcript_features.log"
: > "$LOG"

# ----------------
# Checks
# ----------------
req() { [[ -s "$1" ]] || { echo "[ERROR] Missing/empty: $1" | tee -a "$LOG"; exit 1; }; }

req "$GENOME"
req "$DOMAINS_BED"
req "$TX_EXONS_BED"
req "$CDS_BED"
req "$UTR_BED"

# ----------------
# Helpers
# ----------------
stream_bed() {
  local f="$1"
  if [[ "$f" == *.gz ]]; then
    gzip -dc "$f"
  else
    cat "$f"
  fi
}

# Sort domains once (genome order) for chromsweep speed
DOMAINS_SORTED="context/assets/domains.sorted.bed.gz"
echo "[INFO] Sorting domains with genome order -> $DOMAINS_SORTED" | tee -a "$LOG"
bedtools sort -g "$GENOME" -i <(stream_bed "$DOMAINS_BED") | gzip -c > "$DOMAINS_SORTED"
req "$DOMAINS_SORTED"

# Generic intersect writer (no -wa/-wb; -wo already outputs both)
run_intersect_wo() {
  local b_bed_gz="$1"
  local out="$2"
  local label="$3"
  echo "[INFO] Intersect: ${label}" | tee -a "$LOG"
  echo "       A=$DOMAINS_SORTED" | tee -a "$LOG"
  echo "       B=$b_bed_gz" | tee -a "$LOG"
  echo "       OUT=$out" | tee -a "$LOG"

  {
    echo -e "#a_chrom\ta_start\ta_end\ta_name\ta_score\ta_strand\tb_chrom\tb_start\tb_end\tb_name\tb_score\tb_strand\toverlap_bp"
    bedtools intersect -sorted -g "$GENOME" -s -wo \
      -a <(gzip -dc "$DOMAINS_SORTED") \
      -b <(gzip -dc "$b_bed_gz")
  } > "$out"
}

# ----------------
# 1) Domains vs transcript exons (transcript-aware)
# ----------------
OUT_TX_EXONS="${OUTDIR}/domains_vs_gencode_tx_exons.same_strand.wo.tsv"
run_intersect_wo "$TX_EXONS_BED" "$OUT_TX_EXONS" "domains vs transcript exons"

# ----------------
# 2) Domains vs CDS (transcript-aware)
# ----------------
OUT_CDS="${OUTDIR}/domains_vs_gencode_cds.same_strand.wo.tsv"
run_intersect_wo "$CDS_BED" "$OUT_CDS" "domains vs CDS"

# ----------------
# 3) Domains vs UTR (transcript-aware, but UTR not split 5'/3' yet)
# ----------------
OUT_UTR="${OUTDIR}/domains_vs_gencode_utr.same_strand.wo.tsv"
run_intersect_wo "$UTR_BED" "$OUT_UTR" "domains vs UTR (unsplit)"

# ----------------
# 4) OPTIONAL: Domains vs merged gene exons (gene-level sanity layer)
# ----------------
if [[ -s "$GENE_EXONS_MERGED_BED" ]]; then
  OUT_GENE_EXONS="${OUTDIR}/domains_vs_gencode_gene_exons_merged.same_strand.wo.tsv"
  run_intersect_wo "$GENE_EXONS_MERGED_BED" "$OUT_GENE_EXONS" "domains vs merged gene exons"
else
  echo "[WARN] Skipping merged gene exons; file missing: $GENE_EXONS_MERGED_BED" | tee -a "$LOG"
fi

# ----------------
# 5) OPTIONAL (recommended): split UTR into 5' and 3' per transcript
#     using transcript strand + CDS bounds
#
#     Logic:
#       For each transcript:
#         cds_min_start = min(CDS start), cds_max_end = max(CDS end)
#         On '+' strand:
#            UTR segments with end <= cds_min_start => 5'UTR
#            UTR segments with start >= cds_max_end => 3'UTR
#         On '-' strand:
#            UTR segments with end <= cds_min_start => 3'UTR   (genomic left is transcript 3')
#            UTR segments with start >= cds_max_end => 5'UTR
#
#     Anything overlapping the CDS span is rare/odd; we label it "UTR_internal"
# ----------------
UTR_SPLIT_BED="context/assets/gencode.v48.utr_5p3p.sorted.bed.gz"
echo "[INFO] Building transcript-aware 5'UTR/3'UTR split -> $UTR_SPLIT_BED" | tee -a "$LOG"

python3 - <<'PY'
import gzip
import sys
from collections import defaultdict

cds_path = "bed/gencode.v48.cds.sorted.bed.gz"
utr_path = "bed/gencode.v48.utr.sorted.bed.gz"
out_path = "context/assets/gencode.v48.utr_5p3p.sorted.bed.gz"

# Parse CDS bounds per transcript
# name looks like: ENST|ENSG|GENE|CDS|...
tx_bounds = {}  # tx -> (chrom, strand, min_start, max_end)
with gzip.open(cds_path, "rt") as f:
    for line in f:
        if not line.strip():
            continue
        chrom, start, end, name, score, strand = line.rstrip("\n").split("\t")[:6]
        start = int(start); end = int(end)
        tx = name.split("|")[0]
        if tx not in tx_bounds:
            tx_bounds[tx] = [chrom, strand, start, end]
        else:
            tx_bounds[tx][2] = min(tx_bounds[tx][2], start)
            tx_bounds[tx][3] = max(tx_bounds[tx][3], end)

# Split UTR segments using CDS span
out_lines = []
with gzip.open(utr_path, "rt") as f:
    for line in f:
        if not line.strip():
            continue
        chrom, start, end, name, score, strand = line.rstrip("\n").split("\t")[:6]
        start = int(start); end = int(end)
        parts = name.split("|")
        tx = parts[0]

        # If transcript has no CDS (non-coding), we can’t define 5' vs 3' from CDS bounds.
        # Label as UTR_nocds (still transcript-aware).
        if tx not in tx_bounds:
            new_name = name.replace("|UTR|", "|UTR_nocds|", 1)
            out_lines.append((chrom, start, end, new_name, score, strand))
            continue

        cds_min = tx_bounds[tx][2]
        cds_max = tx_bounds[tx][3]

        # classify
        if strand == "+":
            if end <= cds_min:
                utr_type = "five_prime_UTR"
            elif start >= cds_max:
                utr_type = "three_prime_UTR"
            else:
                utr_type = "UTR_internal"
        else:  # strand == "-"
            if end <= cds_min:
                utr_type = "three_prime_UTR"
            elif start >= cds_max:
                utr_type = "five_prime_UTR"
            else:
                utr_type = "UTR_internal"

        # replace the 4th field token "UTR" with utr_type if present, else append
        if "|UTR|" in name:
            new_name = name.replace("|UTR|", f"|{utr_type}|", 1)
        else:
            new_name = name + f"|{utr_type}"
        out_lines.append((chrom, start, end, new_name, score, strand))

# Write sorted (they should already be sorted, but we keep stable output)
# We'll sort lexicographically here; bedtools sort will enforce genome order later if needed.
out_lines.sort(key=lambda x: (x[0], x[1], x[2], x[5], x[3]))

with gzip.open(out_path, "wt") as out:
    for chrom, start, end, name, score, strand in out_lines:
        out.write(f"{chrom}\t{start}\t{end}\t{name}\t{score}\t{strand}\n")
PY

# Enforce genome sort order for chromsweep consistency
bedtools sort -g "$GENOME" -i <(gzip -dc "$UTR_SPLIT_BED") | gzip -c > "${UTR_SPLIT_BED}.tmp"
mv "${UTR_SPLIT_BED}.tmp" "$UTR_SPLIT_BED"
req "$UTR_SPLIT_BED"

OUT_UTR_SPLIT="${OUTDIR}/domains_vs_gencode_utr_5p3p.same_strand.wo.tsv"
run_intersect_wo "$UTR_SPLIT_BED" "$OUT_UTR_SPLIT" "domains vs UTR (5'/3' split)"

echo "[OK] Done. Outputs:" | tee -a "$LOG"
ls -lh "$OUTDIR"/domains_vs_gencode_* | tee -a "$LOG"
echo "[OK] Log: $LOG" | tee -a "$LOG"