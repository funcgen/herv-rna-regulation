#!/usr/bin/env bash
set -euo pipefail

# ----------------
# INPUTS
# ----------------
GENOME="input/genome/grch38.primary.genome"

# LTRs (plain bed)
LTRS_BED="bed/ltrs_merged.bed"

# Transcript-aware, sorted features (same as you used)
BED_DIR="bed"
TX_EXONS_BED="${BED_DIR}/gencode.v48.exons.sorted.bed.gz"
CDS_BED="${BED_DIR}/gencode.v48.cds.sorted.bed.gz"
UTR_BED="${BED_DIR}/gencode.v48.utr.sorted.bed.gz"

# ----------------
# OUTPUTS
# ----------------
OUTDIR="context/results/overlap_ltrs"
ASSETDIR="context/assets"
LOGDIR="context/logs"
mkdir -p "$OUTDIR" "$ASSETDIR" "$LOGDIR"

LOG="${LOGDIR}/05_overlap_ltrs_vs_gencode_transcript_features.log"
: > "$LOG"

# ----------------
# Checks
# ----------------
req() { [[ -s "$1" ]] || { echo "[ERROR] Missing/empty: $1" | tee -a "$LOG"; exit 1; }; }

req "$GENOME"
req "$LTRS_BED"
req "$TX_EXONS_BED"
req "$CDS_BED"
req "$UTR_BED"

# ----------------
# Sort LTRs once
# ----------------
LTRS_SORTED="${ASSETDIR}/ltrs.sorted.bed.gz"
echo "[INFO] Sorting LTRs with genome order -> $LTRS_SORTED" | tee -a "$LOG"
bedtools sort -g "$GENOME" -i "$LTRS_BED" | gzip -c > "$LTRS_SORTED"
req "$LTRS_SORTED"

# ----------------
# Split UTR into 5'/3' per transcript (same logic you used)
# ----------------
UTR_SPLIT_BED="${ASSETDIR}/gencode.v48.utr_5p3p.sorted.bed.gz"
echo "[INFO] Building transcript-aware 5'UTR/3'UTR split -> $UTR_SPLIT_BED" | tee -a "$LOG"

python3 - <<'PY'
import gzip

cds_path = "bed/gencode.v48.cds.sorted.bed.gz"
utr_path = "bed/gencode.v48.utr.sorted.bed.gz"
out_path = "context/assets/gencode.v48.utr_5p3p.sorted.bed.gz"

# tx -> (chrom, strand, cds_min_start, cds_max_end)
tx_bounds = {}
with gzip.open(cds_path, "rt") as f:
    for line in f:
        if not line.strip(): continue
        chrom, start, end, name, score, strand = line.rstrip("\n").split("\t")[:6]
        start = int(start); end = int(end)
        tx = name.split("|")[0]
        if tx not in tx_bounds:
            tx_bounds[tx] = [chrom, strand, start, end]
        else:
            tx_bounds[tx][2] = min(tx_bounds[tx][2], start)
            tx_bounds[tx][3] = max(tx_bounds[tx][3], end)

out_lines = []
with gzip.open(utr_path, "rt") as f:
    for line in f:
        if not line.strip(): continue
        chrom, start, end, name, score, strand = line.rstrip("\n").split("\t")[:6]
        start = int(start); end = int(end)
        tx = name.split("|")[0]

        if tx not in tx_bounds:
            # non-coding transcript
            if "|UTR|" in name:
                new_name = name.replace("|UTR|", "|UTR_nocds|", 1)
            else:
                new_name = name + "|UTR_nocds"
            out_lines.append((chrom, start, end, new_name, score, strand))
            continue

        cds_min = tx_bounds[tx][2]
        cds_max = tx_bounds[tx][3]

        if strand == "+":
            if end <= cds_min:
                utr_type = "five_prime_UTR"
            elif start >= cds_max:
                utr_type = "three_prime_UTR"
            else:
                utr_type = "UTR_internal"
        else:
            if end <= cds_min:
                utr_type = "three_prime_UTR"
            elif start >= cds_max:
                utr_type = "five_prime_UTR"
            else:
                utr_type = "UTR_internal"

        if "|UTR|" in name:
            new_name = name.replace("|UTR|", f"|{utr_type}|", 1)
        else:
            new_name = name + f"|{utr_type}"
        out_lines.append((chrom, start, end, new_name, score, strand))

out_lines.sort(key=lambda x: (x[0], x[1], x[2], x[5], x[3]))
with gzip.open(out_path, "wt") as out:
    for chrom, start, end, name, score, strand in out_lines:
        out.write(f"{chrom}\t{start}\t{end}\t{name}\t{score}\t{strand}\n")
PY

bedtools sort -g "$GENOME" -i <(gzip -dc "$UTR_SPLIT_BED") | gzip -c > "${UTR_SPLIT_BED}.tmp"
mv "${UTR_SPLIT_BED}.tmp" "$UTR_SPLIT_BED"
req "$UTR_SPLIT_BED"

UTR3_BED="${ASSETDIR}/gencode.v48.utr3.sorted.bed.gz"
echo "[INFO] Extracting 3'UTR only -> $UTR3_BED" | tee -a "$LOG"
gzip -dc "$UTR_SPLIT_BED" | awk -F'\t' 'BEGIN{OFS="\t"} $4 ~ /\|three_prime_UTR\|/ {print $0}' \
  | bedtools sort -g "$GENOME" -i - | gzip -c > "$UTR3_BED"
req "$UTR3_BED"

# ----------------
# Build last exon per transcript (from exon bed)
# ----------------
LAST_EXON_BED="${ASSETDIR}/gencode.v48.last_exon_per_tx.sorted.bed.gz"
echo "[INFO] Building last exon per transcript -> $LAST_EXON_BED" | tee -a "$LOG"

python3 - <<'PY'
import gzip

exons_path = "bed/gencode.v48.exons.sorted.bed.gz"
out_path = "context/assets/gencode.v48.last_exon_per_tx.sorted.bed.gz"

best = {}  # tx -> (chrom, start, end, name, score, strand)

with gzip.open(exons_path, "rt") as f:
    for line in f:
        if not line.strip(): continue
        chrom, start, end, name, score, strand = line.rstrip("\n").split("\t")[:6]
        start = int(start); end = int(end)
        tx = name.split("|")[0]

        if tx not in best:
            best[tx] = (chrom, start, end, name, score, strand)
            continue

        chrom0, s0, e0, n0, sc0, st0 = best[tx]
        # '+' => last exon has max end; '-' => last exon has min start
        if strand == "+":
            if end > e0:
                best[tx] = (chrom, start, end, name, score, strand)
        else:
            if start < s0:
                best[tx] = (chrom, start, end, name, score, strand)

rows = list(best.values())
rows.sort(key=lambda x: (x[0], x[1], x[2], x[5], x[3]))

with gzip.open(out_path, "wt") as w:
    for chrom, start, end, name, score, strand in rows:
        if "|exon|" in name:
            new_name = name.replace("|exon|", "|last_exon|", 1)
        else:
            new_name = name + "|last_exon"
        w.write(f"{chrom}\t{start}\t{end}\t{new_name}\t{score}\t{strand}\n")
PY

bedtools sort -g "$GENOME" -i <(gzip -dc "$LAST_EXON_BED") | gzip -c > "${LAST_EXON_BED}.tmp"
mv "${LAST_EXON_BED}.tmp" "$LAST_EXON_BED"
req "$LAST_EXON_BED"

# ----------------
# Generic intersect writer: ANTISENSE (-S)
# ----------------
run_intersect_wo_antisense() {
  local b_bed_gz="$1"
  local out="$2"
  local label="$3"
  echo "[INFO] Intersect antisense: ${label}" | tee -a "$LOG"
  echo "       A=$LTRS_SORTED" | tee -a "$LOG"
  echo "       B=$b_bed_gz" | tee -a "$LOG"
  echo "       OUT=$out" | tee -a "$LOG"
  {
    echo -e "#a_chrom\ta_start\ta_end\ta_name\ta_score\ta_strand\tb_chrom\tb_start\tb_end\tb_name\tb_score\tb_strand\toverlap_bp"
    bedtools intersect -sorted -g "$GENOME" -S -wo \
      -a <(gzip -dc "$LTRS_SORTED") \
      -b <(gzip -dc "$b_bed_gz")
  } > "$out"
}

# ----------------
# Write overlap TSVs
# ----------------
run_intersect_wo_antisense "$TX_EXONS_BED"  "${OUTDIR}/ltrs_vs_gencode_tx_exons.antisense.wo.tsv" "LTRs vs transcript exons"
run_intersect_wo_antisense "$CDS_BED"      "${OUTDIR}/ltrs_vs_gencode_cds.antisense.wo.tsv"      "LTRs vs CDS"
run_intersect_wo_antisense "$UTR_SPLIT_BED" "${OUTDIR}/ltrs_vs_gencode_utr_5p3p.antisense.wo.tsv" "LTRs vs UTR (5'/3' split)"
run_intersect_wo_antisense "$UTR3_BED"     "${OUTDIR}/ltrs_vs_gencode_utr3.antisense.wo.tsv"     "LTRs vs 3'UTR only"
run_intersect_wo_antisense "$LAST_EXON_BED" "${OUTDIR}/ltrs_vs_gencode_last_exon.antisense.wo.tsv" "LTRs vs last exon (per tx)"

echo "[OK] Done. Outputs:" | tee -a "$LOG"
ls -lh "$OUTDIR"/ltrs_vs_gencode_*.tsv | tee -a "$LOG"
echo "[OK] Log: $LOG" | tee -a "$LOG"