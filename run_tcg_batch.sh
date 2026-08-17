#!/usr/bin/env bash
# Run the TCG executable over all images under a folder (recursively).
#
# Image, .edg, and .cem roots are specified separately. Matching uses the
# image's relative subdirectory plus a filename pattern where '*' is the
# image stem (needed when multiple .edg/.cem variants exist per image).
#
# Usage:
#   ./run_tcg_batch.sh --image-dir DIR --edg-dir DIR --cem-dir DIR [options]
#   ./run_tcg_batch.sh <image_dir> <edg_dir> <cem_dir> [options]
#
# Example (imagenette class folders, TO edges + dborl contours):
#   ./run_tcg_batch.sh \
#     --image-dir .../imagenette/images \
#     --edg-dir   .../imagenette/edges \
#     --cem-dir   .../imagenette/contours \
#     --edg-name '*_to.edg' --cem-name '*_to_dborl.cem' \
#     -e JPEG -o .../out
#
# Options:
#   -i, --image-dir DIR    Root directory of input images (required)
#       --edg-dir DIR      Root directory of .edg files (required)
#       --cem-dir DIR      Root directory of .cem files (required)
#   -o, --output-dir DIR   Output directory (default: <image_dir>/tcg_cpp_out)
#                          Mirrors the image subdirectory layout.
#   -e, --ext EXT          Image extension to match, case-sensitive
#                          (default: JPEG; leading '.' optional)
#       --edg-name PAT     Edge filename pattern; '*' = image stem
#                          (default: '*.edg', e.g. '*_to.edg')
#       --cem-name PAT     Contour filename pattern; '*' = image stem
#                          (default: '*.cem', e.g. '*_to_dborl.cem')
#   -f, --format FMT       cem | cemv | both (default: both)
#   -b, --binary PATH      Path to TCG executable
#                          (default: <repo>/CPP/build/TCG)
#   -n, --dry-run          Print commands without running them
#   -h, --help             Show this help

set -euo pipefail

SCRIPT_DIR="$(cd "$(dirname "${BASH_SOURCE[0]}")" && pwd)"
DEFAULT_BINARY="${SCRIPT_DIR}/CPP/build/TCG"
DEFAULT_IMAGE_EXT="JPEG"
DEFAULT_EDG_NAME="*.edg"
DEFAULT_CEM_NAME="*.cem"

verbose_usage() {
  cat <<'EOF'
Run the TCG executable over all images under a folder (recursively).

Image, .edg, and .cem roots are specified separately. For an image at
  <image_dir>/<rel>/<stem>.EXT
the script looks for
  <edg_dir>/<rel>/<edg-name with * replaced by stem>
  <cem_dir>/<rel>/<cem-name with * replaced by stem>
and writes outputs under
  <output_dir>/<rel>/<stem>_tcg_cpp.{cem,cemv}

Usage:
  ./run_tcg_batch.sh --image-dir DIR --edg-dir DIR --cem-dir DIR [options]
  ./run_tcg_batch.sh <image_dir> <edg_dir> <cem_dir> [options]

Example:
  ./run_tcg_batch.sh \
    --image-dir /path/to/images \
    --edg-dir   /path/to/edges \
    --cem-dir   /path/to/contours \
    --edg-name '*_to.edg' \
    --cem-name '*_to_dborl.cem' \
    -e JPEG -o /path/to/out

Options:
  -i, --image-dir DIR    Root directory of input images   (required)
      --edg-dir DIR      Root directory of .edg files     (required)
      --cem-dir DIR      Root directory of .cem files     (required)
  -o, --output-dir DIR   Output directory                 (default: <cem-dir>/tcg_cpp_out)
                         Mirrors the image subdirectory layout.
  -e, --ext EXT          Image extension to match         (default: JPEG; case-sensitive;
                         leading '.' optional, e.g. JPEG or .jpg)
      --edg-name PAT     Edge filename pattern; '*' = stem
                         (default: '*.edg', e.g. '*_to.edg' or '*_se.edg')
      --cem-name PAT     Contour filename pattern; '*' = stem
                         (default: '*.cem', e.g. '*_to_dborl.cem' or '*_se_tcg.cem')
  -f, --format FMT       cem | cemv | both                (default: both)
  -b, --binary PATH      Path to TCG executable           (default: <repo>/CPP/build/TCG)
  -n, --dry-run          Print commands without running them
  -h, --help             Show this help
EOF
}

die() {
  echo "error: $*" >&2
  exit 1
}

require_dir() {
  local label="$1"
  local dir="$2"
  [[ -n "${dir}" ]] || die "${label} is required"
  [[ -d "${dir}" ]] || die "${label} not found: ${dir}"
  (cd "${dir}" && pwd)
}

# Validate a '*'-stem filename pattern and optionally check its extension.
validate_name_pattern() {
  local label="$1"
  local pat="$2"
  local expect_ext="$3" # e.g. .edg or .cem

  [[ -n "${pat}" ]] || die "${label} must not be empty"
  [[ "${pat}" == *"/"* ]] && die "${label} must be a filename pattern, not a path: ${pat}"
  [[ "${pat}" == *"*"* ]] || die "${label} must contain '*' for the image stem (got: ${pat})"
  local stars="${pat//[^*]/}"
  [[ "${#stars}" -eq 1 ]] || die "${label} must contain exactly one '*' (got: ${pat})"
  [[ "${pat}" == *"${expect_ext}" ]] || die "${label} must end with ${expect_ext} (got: ${pat})"
}

# Replace the single '*' in pattern with stem.
expand_name_pattern() {
  local pat="$1"
  local stem="$2"
  echo "${pat/\*/${stem}}"
}

IMAGE_DIR=""
EDG_DIR=""
CEM_DIR=""
OUTPUT_DIR=""
IMAGE_EXT="${DEFAULT_IMAGE_EXT}"
EDG_NAME="${DEFAULT_EDG_NAME}"
CEM_NAME="${DEFAULT_CEM_NAME}"
FORMAT="both"
BINARY="${DEFAULT_BINARY}"
DRY_RUN=0
POSITIONAL=()

while [[ $# -gt 0 ]]; do
  case "$1" in
    -i|--image-dir)
      [[ $# -ge 2 ]] || die "$1 requires an argument"
      IMAGE_DIR="$2"
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
    -o|--output-dir)
      [[ $# -ge 2 ]] || die "$1 requires an argument"
      OUTPUT_DIR="$2"
      shift 2
      ;;
    -e|--ext)
      [[ $# -ge 2 ]] || die "$1 requires an argument"
      IMAGE_EXT="$2"
      shift 2
      ;;
    --edg-name)
      [[ $# -ge 2 ]] || die "$1 requires an argument"
      EDG_NAME="$2"
      shift 2
      ;;
    --cem-name)
      [[ $# -ge 2 ]] || die "$1 requires an argument"
      CEM_NAME="$2"
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
      POSITIONAL+=("$1")
      shift
      ;;
  esac
done

# Allow: ./run_tcg_batch.sh <image_dir> <edg_dir> <cem_dir> ...
if [[ ${#POSITIONAL[@]} -gt 0 ]]; then
  if [[ -z "${IMAGE_DIR}" && ${#POSITIONAL[@]} -ge 1 ]]; then
    IMAGE_DIR="${POSITIONAL[0]}"
    POSITIONAL=("${POSITIONAL[@]:1}")
  fi
  if [[ -z "${EDG_DIR}" && ${#POSITIONAL[@]} -ge 1 ]]; then
    EDG_DIR="${POSITIONAL[0]}"
    POSITIONAL=("${POSITIONAL[@]:1}")
  fi
  if [[ -z "${CEM_DIR}" && ${#POSITIONAL[@]} -ge 1 ]]; then
    CEM_DIR="${POSITIONAL[0]}"
    POSITIONAL=("${POSITIONAL[@]:1}")
  fi
fi
[[ ${#POSITIONAL[@]} -eq 0 ]] || die "unexpected argument: ${POSITIONAL[0]}"

if [[ -z "${IMAGE_DIR}" || -z "${EDG_DIR}" || -z "${CEM_DIR}" ]]; then
  verbose_usage
  die "--image-dir, --edg-dir, and --cem-dir are all required (or pass them as three positional args)"
fi

# Allow "-e .JPEG" or "-e JPEG"; matching remains case-sensitive.
IMAGE_EXT="${IMAGE_EXT#.}"
[[ -n "${IMAGE_EXT}" ]] || die "image extension (-e/--ext) must not be empty"

validate_name_pattern "--edg-name" "${EDG_NAME}" ".edg"
validate_name_pattern "--cem-name" "${CEM_NAME}" ".cem"

IMAGE_DIR="$(require_dir "image directory" "${IMAGE_DIR}")"
EDG_DIR="$(require_dir "edg directory" "${EDG_DIR}")"
CEM_DIR="$(require_dir "cem directory" "${CEM_DIR}")"
OUTPUT_DIR="${OUTPUT_DIR:-${CEM_DIR}/tcg_cpp_out}"

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

# Recursively collect images with the requested extension under IMAGE_DIR.
mapfile -d '' images < <(find "${IMAGE_DIR}" -type f -name "*.${IMAGE_EXT}" -print0 | sort -z)

if [[ ${#images[@]} -eq 0 ]]; then
  die "no .${IMAGE_EXT} images found under ${IMAGE_DIR}"
fi

echo "TCG binary : ${BINARY}"
echo "Image dir  : ${IMAGE_DIR}"
echo "Edg dir    : ${EDG_DIR}"
echo "Cem dir    : ${CEM_DIR}"
echo "Output dir : ${OUTPUT_DIR}"
echo "Format     : ${FORMAT}"
echo "Extension  : .${IMAGE_EXT} only"
echo "Edg name   : ${EDG_NAME}  (* = image stem)"
echo "Cem name   : ${CEM_NAME}  (* = image stem)"
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

  # Mirror subdirectory layout across edg/cem/output roots.
  rel="${img#${IMAGE_DIR}/}"
  rel_dir="$(dirname "${rel}")"
  if [[ "${rel_dir}" == "." ]]; then
    edg_subdir="${EDG_DIR}"
    cem_subdir="${CEM_DIR}"
    out_subdir="${OUTPUT_DIR}"
  else
    edg_subdir="${EDG_DIR}/${rel_dir}"
    cem_subdir="${CEM_DIR}/${rel_dir}"
    out_subdir="${OUTPUT_DIR}/${rel_dir}"
  fi

  edg_base="$(expand_name_pattern "${EDG_NAME}" "${stem}")"
  cem_base="$(expand_name_pattern "${CEM_NAME}" "${stem}")"
  edg="${edg_subdir}/${edg_base}"
  cem="${cem_subdir}/${cem_base}"
  mkdir -p "${out_subdir}"
  out_stem="${out_subdir}/${stem}_tcg_cpp"

  if [[ ! -f "${edg}" ]]; then
    echo "[SKIP] ${rel}: missing ${edg}"
    skip=$((skip + 1))
    continue
  fi
  if [[ ! -f "${cem}" ]]; then
    echo "[SKIP] ${rel}: missing ${cem}"
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

printf 'Done. ok=%d skip=%d fail=%d\n' "${ok}" "${skip}" "${fail}"
[[ "${fail}" -eq 0 ]]
