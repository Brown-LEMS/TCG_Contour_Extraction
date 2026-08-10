#!/usr/bin/env bash
# Run the TCG executable over all .JPEG images under a folder (recursively).
#
# Expects each image stem to have matching <stem>.edg and <stem>.cem (or
# <stem>_se_tcg.cem) in the same folder as the image (or in --edg-dir /
# --cem-dir if provided).
#
# Usage:
#   ./run_tcg_batch.sh <image_dir> [options]
#
# Example (imagenette class folders):
#   ./run_tcg_batch.sh /gpfs/data/bkimia/imagenette/images/
#
# Options:
#   -o, --output-dir DIR   Output directory (default: <image_dir>/tcg_cpp_out)
#                          Mirrors the image subdirectory layout.
#   -f, --format FMT       cem | cemv | both (default: both)
#   -b, --binary PATH      Path to TCG executable
#                          (default: <repo>/CPP/build/TCG)
#       --edg-dir DIR      Directory for .edg files (default: each image's dir)
#       --cem-dir DIR      Directory for .cem files (default: each image's dir)
#   -n, --dry-run          Print commands without running them
#   -h, --help             Show this help

set -euo pipefail

SCRIPT_DIR="$(cd "$(dirname "${BASH_SOURCE[0]}")" && pwd)"
DEFAULT_BINARY="${SCRIPT_DIR}/CPP/build/TCG"

# Only process this extension (case-sensitive for the imagenette dataset)
IMAGE_EXT="JPEG"

verbose_usage() {
  cat <<'EOF'
Run the TCG executable over all .JPEG images under a folder (recursively).

Expects each image stem to have matching .edg and .cem (or _se_tcg.cem) files
in the same folder as the image (or in --edg-dir / --cem-dir if provided).

Usage:
  ./run_tcg_batch.sh <image_dir> [options]

Example:
  ./run_tcg_batch.sh /gpfs/data/bkimia/imagenette/images/

Options:
  -o, --output-dir DIR   Output directory          (default: <image_dir>/tcg_cpp_out)
                         Mirrors the image subdirectory layout.
  -f, --format FMT       cem | cemv | both         (default: both)
  -b, --binary PATH      Path to TCG executable    (default: <repo>/CPP/build/TCG)
      --edg-dir DIR      Directory for .edg files  (default: each image dir)
      --cem-dir DIR      Directory for .cem files  (default: each image dir)
  -n, --dry-run          Print commands without running them
  -h, --help             Show this help
EOF
}

die() {
  echo "error: $*" >&2
  exit 1
}

# Resolve input .cem for a stem: prefer <stem>.cem if present
resolve_cem() {
  local cem_dir="$1"
  local stem="$2"
  if [[ -f "${cem_dir}/${stem}.cem" ]]; then
    echo "${cem_dir}/${stem}.cem"
  else
    echo ""
  fi
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
if [[ -n "${EDG_DIR}" ]]; then
  [[ -d "${EDG_DIR}" ]] || die "edg directory not found: ${EDG_DIR}"
  EDG_DIR="$(cd "${EDG_DIR}" && pwd)"
fi
if [[ -n "${CEM_DIR}" ]]; then
  [[ -d "${CEM_DIR}" ]] || die "cem directory not found: ${CEM_DIR}"
  CEM_DIR="$(cd "${CEM_DIR}" && pwd)"
fi
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

# Recursively collect only .JPEG files under IMAGE_DIR (all sub-folders).
mapfile -d '' images < <(find "${IMAGE_DIR}" -type f -name "*.${IMAGE_EXT}" -print0 | sort -z)

if [[ ${#images[@]} -eq 0 ]]; then
  die "no .${IMAGE_EXT} images found under ${IMAGE_DIR}"
fi

echo "TCG binary : ${BINARY}"
echo "Image dir  : ${IMAGE_DIR}"
echo "Edg dir    : ${EDG_DIR:-<per-image directory>}"
echo "Cem dir    : ${CEM_DIR:-<per-image directory>}"
echo "Output dir : ${OUTPUT_DIR}"
echo "Format     : ${FORMAT}"
echo "Extension  : .${IMAGE_EXT} only"
echo "Images     : ${#images[@]}"
echo

ok=0
skip=0
fail=0

for img in "${images[@]}"; do
  # Drop trailing empty entry that mapfile can leave with -d ''.
  [[ -n "${img}" ]] || continue

  base="$(basename "${img}")"
  stem="${base%.*}"
  img_dir="$(cd "$(dirname "${img}")" && pwd)"
  edg_dir="${EDG_DIR:-${img_dir}}"
  cem_dir="${CEM_DIR:-${img_dir}}"
  edg="${edg_dir}/${stem}.edg"
  cem="$(resolve_cem "${cem_dir}" "${stem}")"

  # Mirror subdirectory layout under OUTPUT_DIR when images live in subfolders.
  rel="${img#${IMAGE_DIR}/}"
  rel_dir="$(dirname "${rel}")"
  if [[ "${rel_dir}" == "." ]]; then
    out_subdir="${OUTPUT_DIR}"
  else
    out_subdir="${OUTPUT_DIR}/${rel_dir}"
  fi
  mkdir -p "${out_subdir}"
  out_stem="${out_subdir}/${stem}_tcg_cpp"

  if [[ ! -f "${edg}" ]]; then
    echo "[SKIP] ${rel}: missing ${edg}"
    skip=$((skip + 1))
    continue
  fi
  if [[ -z "${cem}" ]]; then
    echo "[SKIP] ${rel}: missing ${cem_dir}/${stem}.cem (or ${stem}_se_tcg.cem)"
    skip=$((skip + 1))
    continue
  fi

  case "${FORMAT}" in
    cem)  out_arg="${out_stem}.cem" ;;
    cemv) out_arg="${out_stem}.cemv" ;;
    both) out_arg="${out_stem}.cem" ;;
  esac

  cmd=("${BINARY}" "${edg}" "${cem}" "${img}" "${FORMAT}" "${out_arg}")
  echo "[RUN] ${rel}"
  if [[ "${DRY_RUN}" -eq 1 ]]; then
    printf '  %q' "${cmd[@]}"
    echo
    ok=$((ok + 1))
    continue
  fi

  if "${cmd[@]}"; then
    ok=$((ok + 1))
  else
    echo "[FAIL] ${rel}" >&2
    fail=$((fail + 1))
  fi
  echo
done

echo "Done. ok=${ok} skip=${skip} fail=${fail}"
[[ "${fail}" -eq 0 ]]
