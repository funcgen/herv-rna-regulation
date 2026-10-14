#!/usr/bin/env bash
# =============================================================================
# 01_sources : fetch and verify published upstream release files
# =============================================================================
# Downloads every published source the canonical workflow depends on, directly
# from its Zenodo record, and verifies each against the md5 published by Zenodo.
#
# Files are stored under their PUBLICATION-FACING names. Historical working
# aliases (ERV_internal.fasta, ERV_GyDB_v6_domains.bed, ...) are recorded in
# manifests/upstream_releases.tsv for provenance only and are never exposed as
# the release interface.
#
# Fails loudly and stops on any checksum mismatch. Never overwrites a file that
# already verifies.
#
# Usage:  ./fetch_upstream_releases.sh [--dry-run]
# Env:    NCRNA_SYSTEM, NCRNA_RELEASE_ROOT (see config/paths.yaml)
# =============================================================================
set -euo pipefail
HERE="$(cd "$(dirname "${BASH_SOURCE[0]}")" && pwd)"
source "${HERE}/../../config/lib/common.sh"

DRY_RUN=0
[[ "${1:-}" == "--dry-run" ]] && DRY_RUN=1

MANIFEST="${SOURCES_DIR}/manifests/upstream_releases.tsv"
mkdir -p "${LOG_DIR}"
export NCRNA_LOG_FILE="${LOG_DIR}/fetch_upstream_releases.$(date +%Y%m%dT%H%M%S).log"

log "system=${NCRNA_SYSTEM} release_root=${RELEASE_ROOT}"
log "manifest=${MANIFEST}"

dest_dir_for() {
  case "$1" in
    SRC03|SRC04|SRC05|SRC06) echo "${SOURCES_DIR}/frozen/herv" ;;
    SRC07|SRC08|SRC09|SRC10) echo "${SOURCES_DIR}/frozen/ltr" ;;
    SRC12a|SRC12b)           echo "${SOURCES_DIR}/frozen/attract" ;;
    SRC13)                   echo "${SOURCES_DIR}/derived" ;;
    *) die "no destination configured for source id: $1" ;;
  esac
}

n_ok=0; n_fetched=0
while IFS=$'\t' read -r sid rec doi fname md5 size alias target; do
  [[ "$sid" == "source_id" ]] && continue
  [[ -z "$sid" ]] && continue
  if [[ "$target" != "$NCRNA_SYSTEM" ]]; then
    log "SKIP ${sid} ${fname} (targets '${target}', this system is '${NCRNA_SYSTEM}')"
    continue
  fi
  dest_dir="$(dest_dir_for "$sid")"; mkdir -p "$dest_dir"
  dest="${dest_dir}/${fname}"
  url="https://zenodo.org/api/records/${rec}/files/${fname}/content"
  if [[ "$DRY_RUN" == "1" ]]; then
    log "DRY-RUN would fetch ${sid} ${fname} (${size} bytes) -> ${dest}"
    continue
  fi
  if [[ -f "$dest" ]] && [[ "$(md5sum "$dest" | cut -d' ' -f1)" == "$md5" ]]; then
    log "OK   ${sid} ${fname} already present and verified"
    n_ok=$((n_ok+1)); continue
  fi
  safe_fetch "$url" "$dest" "$md5"
  actual_size="$(stat -c%s "$dest")"
  [[ "$actual_size" == "$size" ]] || die "size mismatch for ${fname}: expected ${size}, got ${actual_size}"
  n_fetched=$((n_fetched+1))
done < "$MANIFEST"

log "done: ${n_fetched} fetched, ${n_ok} already verified"
