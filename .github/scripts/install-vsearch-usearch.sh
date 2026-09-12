#!/usr/bin/env bash
# Download platform binaries for vsearch and open-source usearch12 into a
# directory on PATH for GitHub Actions (and local smoke-testing).
set -euo pipefail

DEST="${OPTIMOTU_TOOLS_DIR:-${HOME}/.local/optimotu-tools}"
mkdir -p "${DEST}"

VSEARCH_VERSION="${VSEARCH_VERSION:-2.31.0}"
# Open-source USEARCH 12; currently only a beta release with binaries.
# https://github.com/rcedgar/usearch12
USEARCH_TAG="${USEARCH_TAG:-v12.0-beta1}"
USEARCH_VERSION="${USEARCH_VERSION:-12.0-beta}"

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

# --- usearch12 ---------------------------------------------------------------
case "${OS}" in
  Linux)
    case "${ARCH}" in
      x86_64|amd64)
        US_ASSET="usearch_linux_x86_${USEARCH_VERSION}"
        US_NAME="usearch"
        ;;
      aarch64|arm64)
        US_ASSET="usearch_linux_arch64_${USEARCH_VERSION}"
        US_NAME="usearch"
        ;;
      *)
        echo "Unsupported Linux architecture for usearch: ${ARCH}" >&2
        exit 1
        ;;
    esac
    ;;
  Darwin)
    case "${ARCH}" in
      x86_64|amd64)
        US_ASSET="usearch_osx_x86_${USEARCH_VERSION}"
        US_NAME="usearch"
        ;;
      arm64|aarch64)
        US_ASSET="usearch_osx_m_${USEARCH_VERSION}"
        US_NAME="usearch"
        ;;
      *)
        echo "Unsupported macOS architecture for usearch: ${ARCH}" >&2
        exit 1
        ;;
    esac
    ;;
  MINGW*|MSYS*|CYGWIN*)
    US_ASSET="usearch_win_${USEARCH_VERSION}.exe"
    US_NAME="usearch.exe"
    ;;
  *)
    echo "Unsupported OS for usearch: ${OS}" >&2
    exit 1
    ;;
esac

US_URL="https://github.com/rcedgar/usearch12/releases/download/${USEARCH_TAG}/${US_ASSET}"
echo "Installing usearch ${USEARCH_TAG} from ${US_URL}"
download "${US_URL}" "${DEST}/${US_NAME}"
chmod +x "${DEST}/${US_NAME}"

# --- PATH --------------------------------------------------------------------
if [[ -n "${GITHUB_PATH:-}" ]]; then
  echo "${DEST}" >> "${GITHUB_PATH}"
else
  export PATH="${DEST}:${PATH}"
fi

echo "Installed tools in ${DEST}:"
ls -l "${DEST}"
if is_windows; then
  "${DEST}/vsearch.exe" --version
  test -x "${DEST}/usearch.exe"
else
  "${DEST}/vsearch" --version
  test -x "${DEST}/usearch"
fi
# usearch12 does not accept a bare -version/--version flag the way older
# binaries did; presence + execute bit is enough for CI PATH setup.
