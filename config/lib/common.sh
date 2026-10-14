#!/usr/bin/env bash
# Shared helpers for release scripts. Source this; do not execute it.
# Provides: release_root, log, die, verify_md5, verify_sha256, record_provenance

set -euo pipefail

_detect_system() {
  if [[ -n "${NCRNA_SYSTEM:-}" ]]; then echo "$NCRNA_SYSTEM"; return; fi
  case "$(hostname -f 2>/dev/null || hostname)" in
    *cnag*|*hpc*|login*|cn*) echo cnag ;;
    *) echo local ;;
  esac
}

NCRNA_SYSTEM="$(_detect_system)"
if [[ -n "${NCRNA_RELEASE_ROOT:-}" ]]; then
  RELEASE_ROOT="$NCRNA_RELEASE_ROOT"
elif [[ "$NCRNA_SYSTEM" == "cnag" ]]; then
  RELEASE_ROOT="/path/to/herv-rna-regulation"
else
  RELEASE_ROOT="$HOME/ncRNA"
fi
export NCRNA_SYSTEM RELEASE_ROOT

SOURCES_DIR="${RELEASE_ROOT}/01_sources"
LOG_DIR="${SOURCES_DIR}/logs"

log()  { printf '[%s] %s\n' "$(date -Is)" "$*" | tee -a "${NCRNA_LOG_FILE:-/dev/null}"; }
die()  { printf '[%s] ERROR: %s\n' "$(date -Is)" "$*" >&2 | tee -a "${NCRNA_LOG_FILE:-/dev/null}" >&2; exit 1; }

# verify_md5 <file> <expected_md5>  -- fails loudly on mismatch
verify_md5() {
  local f="$1" want="$2" got
  [[ -f "$f" ]] || die "missing file: $f"
  got="$(md5sum "$f" | cut -d' ' -f1)"
  if [[ "$got" != "$want" ]]; then
    die "CHECKSUM MISMATCH for $f
       expected md5 $want
       observed md5 $got
     The file does NOT match the published release. Refusing to continue."
  fi
  log "OK md5 $want  $(basename "$f")"
}

# verify_sha256 <file> <expected_sha256>
verify_sha256() {
  local f="$1" want="$2" got
  [[ -f "$f" ]] || die "missing file: $f"
  got="$(sha256sum "$f" | cut -d' ' -f1)"
  [[ "$got" == "$want" ]] || die "SHA256 MISMATCH for $f (expected $want, got $got)"
  log "OK sha256 ${want:0:16}...  $(basename "$f")"
}

# safe_fetch <url> <dest> <expected_md5> -- never overwrites a verified file
safe_fetch() {
  local url="$1" dest="$2" want="$3"
  if [[ -f "$dest" ]] && [[ "$(md5sum "$dest" | cut -d' ' -f1)" == "$want" ]]; then
    log "SKIP already present and verified: $(basename "$dest")"
    return 0
  fi
  [[ -f "$dest" ]] && log "WARN existing file failed verification, re-fetching: $dest"
  log "FETCH $url"
  curl -sSL --fail --retry 3 -o "${dest}.part" "$url" || die "download failed: $url"
  mv "${dest}.part" "$dest"
  verify_md5 "$dest" "$want"
}
