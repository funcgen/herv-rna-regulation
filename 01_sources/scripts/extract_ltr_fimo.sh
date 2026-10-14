#!/usr/bin/env bash
# =============================================================================
# 01_sources : verify and extract the published LTR FIMO archive (SRC09)
# =============================================================================
# The published archive HERV_LTR_fimo_results_tsv.tar.xz contains exactly one
# member, merged_fimo_parsed_v4.tsv, proven byte-identical to the table the
# historical analyses consumed (decision D013).
#
#   published archive -> verify archive md5 -> extract -> verify member md5
#
# The 6.3 GB member is extracted on demand rather than stored twice.
#
# Concurrency and atomicity:
#   * an flock guard makes concurrent runs impossible (a second run exits
#     immediately rather than racing the first into the same output file);
#   * extraction goes to a private staging directory and is verified there;
#   * only a fully verified member is moved into its final location, so a
#     partial or interrupted extraction can never be mistaken for a good one.
#
# Usage:  ./extract_ltr_fimo.sh
# =============================================================================
set -euo pipefail
HERE="$(cd "$(dirname "${BASH_SOURCE[0]}")" && pwd)"
source "${HERE}/../../config/lib/common.sh"

MANIFEST="${SOURCES_DIR}/manifests/archive_members.tsv"
mkdir -p "${LOG_DIR}"
export NCRNA_LOG_FILE="${LOG_DIR}/extract_ltr_fimo.$(date +%Y%m%dT%H%M%S).log"

read -r _ archive archive_md5 member member_md5 member_size _ < <(sed -n '2p' "$MANIFEST")

ARCHIVE_PATH="${SOURCES_DIR}/frozen/ltr/${archive}"
EXTRACT_DIR="${SOURCES_DIR}/frozen/ltr/extracted"
MEMBER_PATH="${EXTRACT_DIR}/${member}"
LOCK_FILE="${SOURCES_DIR}/frozen/ltr/.extract.lock"

mkdir -p "$EXTRACT_DIR"
exec 9>"$LOCK_FILE"
if ! flock -n 9; then
  die "another extraction is already running (lock: ${LOCK_FILE}). Refusing to run concurrently."
fi

[[ -f "$ARCHIVE_PATH" ]] || die "archive not present: ${ARCHIVE_PATH} (run fetch_upstream_releases.sh first)"

# Already done and verified?
if [[ -f "$MEMBER_PATH" ]] && [[ "$(md5sum "$MEMBER_PATH" | cut -d' ' -f1)" == "$member_md5" ]]; then
  log "member already extracted and verified: ${member}"
  exit 0
fi
[[ -f "$MEMBER_PATH" ]] && log "WARN existing member failed verification (partial or stale); it will be replaced"

log "verifying archive ${archive}"
verify_md5 "$ARCHIVE_PATH" "$archive_md5"

# One decompression pass only. The "exactly one member" assertion is made on the
# staging directory after extraction, which is equivalent to listing the archive
# first but avoids decompressing 596 MB twice.
STAGE="$(mktemp -d "${EXTRACT_DIR}/.stage.XXXXXX")"
cleanup() { rm -rf "$STAGE"; }
trap cleanup EXIT

log "extracting ${member} (${member_size} bytes) to staging"
tar xJf "$ARCHIVE_PATH" -C "$STAGE"

n_members="$(find "$STAGE" -type f | wc -l)"
[[ "$n_members" == "1" ]] || die "expected exactly 1 archive member, found ${n_members}"

STAGED="${STAGE}/${member}"
[[ -f "$STAGED" ]] || die "expected member not found after extraction: ${member}"
actual_size="$(stat -c%s "$STAGED")"
[[ "$actual_size" == "$member_size" ]] || die "member size mismatch: expected ${member_size}, got ${actual_size}"
log "verifying extracted member"
verify_md5 "$STAGED" "$member_md5"

log "publishing verified member to ${MEMBER_PATH}"
mv -f "$STAGED" "$MEMBER_PATH"
log "done"
