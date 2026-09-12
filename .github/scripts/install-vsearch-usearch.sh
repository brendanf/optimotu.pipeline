#!/usr/bin/env bash
# Download platform binaries for vsearch and USEARCH 11 into a directory on
# PATH for GitHub Actions (and local smoke-testing).
#
# USEARCH 11 (not 12) is required: our helpers and optimotu still call
# calc_distmx, allpairs_global, makeudb_sintax, and fastx_getseqs, which were
# removed in usearch12. Binaries are CC0 from:
# https://github.com/rcedgar/usearch_old_binaries
set -euo pipefail

DEST="${OPTIMOTU_TOOLS_DIR:-${HOME}/.local/optimotu-tools}"
mkdir -p "${DEST}"

VSEARCH_VERSION="${VSEARCH_VERSION:-2.31.0}"
USEARCH_VERSION="${USEARCH_VERSION:-11.0.667}"
USEARCH_BASE_URL="${USEARCH_BASE_URL:-https://raw.githubusercontent.com/rcedgar/usearch_old_binaries/main/bin}"

OS="$(uname -s)"
ARCH="$(uname -m)"

is_windows() {
  case "${OS}" in
    MINGW*|MSYS*|CYGWIN*) return 0 ;;
    *) return 1 ;;
  esac
}

download() {
  local url="$1"
  local out="$2"
  curl -fsSL --retry 3 --retry-delay 2 -o "${out}" "${url}"
}

# --- vsearch -----------------------------------------------------------------
case "${OS}" in
  Linux)
    case "${ARCH}" in
      x86_64|amd64)
        VS_ASSET="vsearch-${VSEARCH_VERSION}-linux-x86_64.tar.gz"
        ;;
      aarch64|arm64)
        VS_ASSET="vsearch-${VSEARCH_VERSION}-linux-aarch64.tar.gz"
        ;;
      *)
        echo "Unsupported Linux architecture for vsearch: ${ARCH}" >&2
        exit 1
        ;;
    esac
    VS_KIND="tar"
    ;;
  Darwin)
    VS_ASSET="vsearch-${VSEARCH_VERSION}-macos-universal.tar.gz"
    VS_KIND="tar"
    ;;
  MINGW*|MSYS*|CYGWIN*)
    VS_ASSET="vsearch-${VSEARCH_VERSION}-win-x86_64.zip"
    VS_KIND="zip"
    ;;
  *)
    echo "Unsupported OS for vsearch: ${OS}" >&2
    exit 1
    ;;
esac

VS_URL="https://github.com/torognes/vsearch/releases/download/v${VSEARCH_VERSION}/${VS_ASSET}"
VS_TMP="$(mktemp -d)"
trap 'rm -rf "${VS_TMP}"' EXIT

echo "Installing vsearch ${VSEARCH_VERSION} from ${VS_URL}"
download "${VS_URL}" "${VS_TMP}/${VS_ASSET}"
if [[ "${VS_KIND}" == "tar" ]]; then
  tar -xzf "${VS_TMP}/${VS_ASSET}" -C "${VS_TMP}"
  VS_BIN="$(find "${VS_TMP}" -type f -path '*/bin/vsearch' | head -n 1)"
  if [[ -z "${VS_BIN}" ]]; then
    echo "Could not locate vsearch binary in ${VS_ASSET}" >&2
    exit 1
  fi
  cp -f "${VS_BIN}" "${DEST}/vsearch"
  chmod +x "${DEST}/vsearch"
else
  unzip -q "${VS_TMP}/${VS_ASSET}" -d "${VS_TMP}"
  # Windows needs zlib/libbz2 DLLs beside the executable.
  VS_BIN_DIR="$(find "${VS_TMP}" -type d -path '*/bin' | head -n 1)"
  if [[ -z "${VS_BIN_DIR}" ]]; then
    echo "Could not locate vsearch bin directory in ${VS_ASSET}" >&2
    exit 1
  fi
  cp -f "${VS_BIN_DIR}"/* "${DEST}/"
  chmod +x "${DEST}/vsearch.exe" 2>/dev/null || true
fi

# --- usearch 11 --------------------------------------------------------------
# Only x86 binaries were published for v11. On Apple Silicon, the osx64 binary
# runs under Rosetta 2 (available on GitHub Actions macos runners).
US_NAME="usearch"
US_ASSET=""
case "${OS}" in
  Linux)
    case "${ARCH}" in
      x86_64|amd64)
        US_ASSET="usearch${USEARCH_VERSION}_i86linux64"
        ;;
      *)
        echo "WARNING: no USEARCH ${USEARCH_VERSION} binary for Linux ${ARCH}; skipping" >&2
        ;;
    esac
    ;;
  Darwin)
    # Prefer 64-bit OSX binary (works natively on Intel; via Rosetta on ARM).
    US_ASSET="usearch${USEARCH_VERSION}_i86osx64"
    ;;
  MINGW*|MSYS*|CYGWIN*)
    US_ASSET="usearch${USEARCH_VERSION}_win64.exe"
    US_NAME="usearch.exe"
    ;;
  *)
    echo "WARNING: no USEARCH ${USEARCH_VERSION} binary for OS ${OS}; skipping" >&2
    ;;
esac

if [[ -n "${US_ASSET}" ]]; then
  US_URL="${USEARCH_BASE_URL}/${US_ASSET}"
  echo "Installing usearch ${USEARCH_VERSION} from ${US_URL}"
  download "${US_URL}" "${DEST}/${US_NAME}"
  chmod +x "${DEST}/${US_NAME}"
fi

# --- PATH --------------------------------------------------------------------
# On Windows, GITHUB_PATH must be a native path (C:\...) so later R/cmd steps
# see the tools. Bash MSYS paths like /c/Users/... are ignored there.
path_for_github() {
  local dir="$1"
  if is_windows && command -v cygpath >/dev/null 2>&1; then
    cygpath -w "${dir}"
  else
    printf '%s\n' "${dir}"
  fi
}

if [[ -n "${GITHUB_PATH:-}" ]]; then
  path_for_github "${DEST}" >> "${GITHUB_PATH}"
else
  export PATH="${DEST}:${PATH}"
fi

echo "Installed tools in ${DEST}:"
ls -l "${DEST}"
if is_windows; then
  "${DEST}/vsearch.exe" --version
  if [[ -x "${DEST}/usearch.exe" ]]; then
    "${DEST}/usearch.exe" -version
  fi
else
  "${DEST}/vsearch" --version
  if [[ -x "${DEST}/usearch" ]]; then
    "${DEST}/usearch" -version
  fi
fi
