#!/usr/bin/env bash
# Run the TCG executable over all images in a folder.
#
# Expects each image stem to have matching <stem>.edg and <stem>.cem in the
# same folder (or in --edg-dir / --cem-dir if provided).
#
# Usage:
#   ./run_tcg_folder.sh <image_dir> [options]
#
# Options:
#   -o, --output-dir DIR   Output directory (default: <image_dir>/tcg_cpp_out)
#   -f, --format FMT       cem | cemv | both (default: both)
#   -b, --binary PATH      Path to TCG executable
#                          (default: <repo>/CPP/build/TCG)
#       --edg-dir DIR      Directory for .edg files (default: image_dir)
#       --cem-dir DIR      Directory for .cem files (default: image_dir)
#   -n, --dry-run          Print commands without running them
#   -h, --help             Show this help

set -euo pipefail

SCRIPT_DIR="$(cd "$(dirname "${BASH_SOURCE[0]}")" && pwd)"
DEFAULT_BINARY="${SCRIPT_DIR}/CPP/build/TCG"

IMAGE_EXTS=("jpg" "jpeg" "png" "bmp" "tif" "tiff" "JPEG" "JPG" "PNG" "BMP" "TIF" "TIFF")

verbose_usage() {
  cat <<'EOF'
Run the TCG executable over all .edg, .cem, and image files in a folder.

Expects each image stem to have matching .edg and .cem files in the same folder
(or in --edg-dir / --cem-dir if provided).

Usage:
  ./run_tcg_batch.sh <image_dir> [options]

Options:
  -o, --output-dir DIR   Output directory         (default: <image_dir>/tcg_cpp_out)
  -f, --format FMT       cem | cemv | both        (default: both)
  -b, --binary PATH      Path to TCG executable   (default: <repo>/CPP/build/TCG)
      --edg-dir DIR      Directory for .edg files (default: image_dir)
      --cem-dir DIR      Directory for .cem files (default: image_dir)
  -n, --dry-run          Print commands without running them
  -h, --help             Show this help
EOF
}

die() {
  echo "error: $*" >&2
  exit 1
}

IMAGE_DIR=""
OUTPUT_DIR=""
FORMAT="both"
BINARY="${DEFAULT_BINARY}"
EDG_DIR=""
CEM_DIR=""
DRY_RUN=0

while [[ $# -gt 0 ]]; do
  case "$1" in
    -o|--output-dir)
      [[ $# -ge 2 ]] || die "$1 requires an argument"
      OUTPUT_DIR="$2"
      shift 2
      ;;
    -f|--format)
      [[ $# -ge 2 ]] || die "$1 requires an argument"
      FORMAT="$2"
      shift 2
      ;;
    -b|--binary)
      [[ $# -ge 2 ]] || die "$1 requires an argument"
      BINARY="$2"
      shift 2
      ;;
    --edg-dir)
      [[ $# -ge 2 ]] || die "$1 requires an argument"
      EDG_DIR="$2"
      shift 2
      ;;
    --cem-dir)
      [[ $# -ge 2 ]] || die "$1 requires an argument"
      CEM_DIR="$2"
      shift 2
      ;;
    -n|--dry-run)
      DRY_RUN=1
      shift
      ;;
    -h|--help)
      verbose_usage
      exit 0
      ;;
    -*)
      die "unknown option: $1"
      ;;
    *)
      if [[ -z "${IMAGE_DIR}" ]]; then
        IMAGE_DIR="$1"
      else
        die "unexpected argument: $1"
      fi
      shift
      ;;
  esac
done

[[ -n "${IMAGE_DIR}" ]] || { verbose_usage; die "image directory is required"; }
[[ -d "${IMAGE_DIR}" ]] || die "image directory not found: ${IMAGE_DIR}"

IMAGE_DIR="$(cd "${IMAGE_DIR}" && pwd)"
EDG_DIR="$(cd "${EDG_DIR:-${IMAGE_DIR}}" && pwd)"
CEM_DIR="$(cd "${CEM_DIR:-${IMAGE_DIR}}" && pwd)"
OUTPUT_DIR="${OUTPUT_DIR:-${IMAGE_DIR}/tcg_cpp_out}"

case "${FORMAT}" in
  cem|cemv|both) ;;
  *) die "format must be cem, cemv, or both (got: ${FORMAT})" ;;
esac

if [[ "${DRY_RUN}" -eq 0 ]]; then
  [[ -x "${BINARY}" ]] || die "TCG binary not found or not executable: ${BINARY}
Build it with:
  cd ${SCRIPT_DIR}/CPP && mkdir -p build && cd build && cmake .. -DCMAKE_BUILD_TYPE=Release && cmake --build . -j"
fi

mkdir -p "${OUTPUT_DIR}"
OUTPUT_DIR="$(cd "${OUTPUT_DIR}" && pwd)"

shopt -s nullglob
images=()
for ext in "${IMAGE_EXTS[@]}"; do
  for f in "${IMAGE_DIR}"/*."${ext}"; do
    images+=("$f")
  done
done

# De-duplicate (e.g. if filesystem is case-insensitive)
declare -A seen=()
unique_images=()
for img in "${images[@]}"; do
  key="$(realpath -m "$img")"
  if [[ -z "${seen[$key]+x}" ]]; then
    seen[$key]=1
    unique_images+=("$img")
  fi
done
images=("${unique_images[@]}")

if [[ ${#images[@]} -eq 0 ]]; then
  die "no images found in ${IMAGE_DIR} (extensions: ${IMAGE_EXTS[*]})"
fi

echo "TCG binary : ${BINARY}"
echo "Image dir  : ${IMAGE_DIR}"
echo "Edg dir    : ${EDG_DIR}"
echo "Cem dir    : ${CEM_DIR}"
echo "Output dir : ${OUTPUT_DIR}"
echo "Format     : ${FORMAT}"
echo "Images     : ${#images[@]}"
echo

ok=0
skip=0
fail=0

for img in "${images[@]}"; do
  base="$(basename "${img}")"
  stem="${base%.*}"
  edg="${EDG_DIR}/${stem}.edg"
  cem="${CEM_DIR}/${stem}.cem"
  out_stem="${OUTPUT_DIR}/${stem}_tcg_cpp"

  if [[ ! -f "${edg}" ]]; then
    echo "[SKIP] ${base}: missing ${edg}"
    skip=$((skip + 1))
    continue
  fi
  if [[ ! -f "${cem}" ]]; then
    echo "[SKIP] ${base}: missing ${cem}"
    skip=$((skip + 1))
    continue
  fi

  case "${FORMAT}" in
    cem)  out_arg="${out_stem}.cem" ;;
    cemv) out_arg="${out_stem}.cemv" ;;
    both) out_arg="${out_stem}.cem" ;;
  esac

  cmd=("${BINARY}" "${edg}" "${cem}" "${img}" "${FORMAT}" "${out_arg}")
  echo "[RUN] ${stem}"
  if [[ "${DRY_RUN}" -eq 1 ]]; then
    printf '  %q' "${cmd[@]}"
    echo
    ok=$((ok + 1))
    continue
  fi

  if "${cmd[@]}"; then
    ok=$((ok + 1))
  else
    echo "[FAIL] ${stem}" >&2
    fail=$((fail + 1))
  fi
  echo
done

echo "Done. ok=${ok} skip=${skip} fail=${fail}"
[[ "${fail}" -eq 0 ]]
