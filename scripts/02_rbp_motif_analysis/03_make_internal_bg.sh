#!/usr/bin/env bash
# 03_make_internal_bg.sh — build Markov background for FIMO on ERV internal FASTA (DNA letters, RNA-like scanning)
#
# Usage:
#   ./03_make_internal_bg.sh <input_fasta> <output_bg.txt> <markov_order>
#
# Example:
#   ./03_make_internal_bg.sh ../../fasta/ERV_internal.fasta ../../rbp/background/bg_internal_m1.sense.txt 1
#
# Notes:
# - Requires fasta-get-markov (MEME suite) in PATH.
# - Although we scan DNA sequences (A/C/G/T), we treat them RNA-like:
#     * scan one orientation at a time (FIMO --norc)
#     * therefore background should also be strand-specific: fasta-get-markov -norc
# - Recommended Markov order for RBP scanning: 1

set -euo pipefail

usage() {
  cat <<EOF
Usage:
  $(basename "$0") <input_fasta> <output_bg.txt> <markov_order>

Arguments:
  input_fasta     FASTA file used to estimate background (e.g. ERV_internal.fasta)
  output_bg.txt   Output background file (e.g. bg_internal_m1.sense.txt)
  markov_order    Markov model order (0 = mono, 1 = di, 2 = tri, ...)

Example:
  $(basename "$0") ../../fasta/ERV_internal.fasta ../../rbp/background/bg_internal_m1.sense.txt 1

Recommended:
  Use markov_order = 1 for ERV internal regions.

This script uses DNA alphabet but strand-specific estimation:
  -dna    force DNA alphabet (A,C,G,T)
  -norc   do NOT combine forward and reverse-complement frequencies
EOF
}

# --- check args ---
if [[ $# -ne 3 ]]; then
  echo "❌ ERROR: wrong number of arguments" >&2
  usage
  exit 1
fi

IN="$1"
OUT="$2"
ORDER="$3"

# --- validate ---
command -v fasta-get-markov >/dev/null 2>&1 || {
  echo "❌ ERROR: fasta-get-markov not found in PATH" >&2
  exit 1
}

[[ -f "$IN" ]] || {
  echo "❌ ERROR: input FASTA not found: $IN" >&2
  exit 1
}

[[ "$ORDER" =~ ^[0-9]+$ ]] || {
  echo "❌ ERROR: markov_order must be a non-negative integer" >&2
  exit 1
}

# --- run ---
echo "[INFO] Generating Markov background model (DNA alphabet, strand-specific)"
echo "       Input FASTA : $IN"
echo "       Output file : $OUT"
echo "       Markov order: $ORDER"
echo "       Flags       : -dna -norc"

# Key fix: -norc prevents strand-combining; -dna keeps A/C/G/T alphabet
fasta-get-markov -dna -norc -m "$ORDER" "$IN" "$OUT"

# --- report ---
echo "[INFO] Done."
echo "[INFO] Background file preview:"
head -n 16 "$OUT" || true

echo "[INFO] Mono-nucleotide frequencies (if present):"
grep -E '^(A|C|G|T)\s' "$OUT" || true

# --- optional sanity note (non-fatal) ---
# If A==T and C==G exactly, you might still be seeing strand-symmetry (can happen naturally, but uncommon).
A=$(awk '$1=="A"{print $2; exit}' "$OUT" 2>/dev/null || echo "")
C=$(awk '$1=="C"{print $2; exit}' "$OUT" 2>/dev/null || echo "")
G=$(awk '$1=="G"{print $2; exit}' "$OUT" 2>/dev/null || echo "")
T=$(awk '$1=="T"{print $2; exit}' "$OUT" 2>/dev/null || echo "")

if [[ -n "$A" && -n "$C" && -n "$G" && -n "$T" ]]; then
  WARN=$(awk -v A="$A" -v T="$T" -v C="$C" -v G="$G" 'BEGIN{
      d1 = (A>T)? A-T : T-A;
      d2 = (C>G)? C-G : G-C;
      if (d1 < 1e-12 && d2 < 1e-12) print "1"; else print "0";
  }')
  if [[ "$WARN" == "1" ]]; then
    echo "[WARN] A==T and C==G in order-0 frequencies."
    echo "       This can happen naturally for large sequence sets, but it also resembles strand-combined models."
    echo "       Confirm your fasta-get-markov version respects -norc (it should) and that the FASTA is as expected."
  fi
fi
