#!/usr/bin/env bash
set -euo pipefail

ROOT_DIR="$(cd "$(dirname "${BASH_SOURCE[0]}")/.." && pwd)"
THIRD_PARTY_DIR="${ROOT_DIR}/third_party"
DEST_DIR="${THIRD_PARTY_DIR}/gvp-pocket_pred"
DEFAULT_REPO_URL="https://github.com/Mickdub/gvp.git"
DEFAULT_GH_REPO="Mickdub/gvp"

METHOD="auto"  # auto | git | gh
REPO_URL="${GVP_POCKET_PRED_REPO_URL:-${DEFAULT_REPO_URL}}"
GH_REPO="${GVP_POCKET_PRED_GH_REPO:-${DEFAULT_GH_REPO}}"

usage() {
  cat <<USAGE
Usage:
  bash scripts/setup_gvp_pocket_pred.sh
  bash scripts/setup_gvp_pocket_pred.sh --method git
  bash scripts/setup_gvp_pocket_pred.sh --method gh
  bash scripts/setup_gvp_pocket_pred.sh --repo-url https://github.com/Mickdub/gvp.git
  bash scripts/setup_gvp_pocket_pred.sh --gh-repo Mickdub/gvp

Default behavior:
  - clones ${DEFAULT_REPO_URL}
  - destination: third_party/gvp-pocket_pred

Environment overrides:
  GVP_POCKET_PRED_REPO_URL=${DEFAULT_REPO_URL}
  GVP_POCKET_PRED_GH_REPO=${DEFAULT_GH_REPO}
USAGE
}

while [[ $# -gt 0 ]]; do
  case "$1" in
    --method)
      METHOD="${2:-}"
      shift 2
      ;;
    --repo-url)
      REPO_URL="${2:-}"
      METHOD="git"
      shift 2
      ;;
    --gh-repo)
      GH_REPO="${2:-}"
      METHOD="gh"
      shift 2
      ;;
    -h|--help)
      usage
      exit 0
      ;;
    http*|git@*)
      REPO_URL="$1"
      METHOD="git"
      shift
      ;;
    */*)
      GH_REPO="$1"
      METHOD="gh"
      shift
      ;;
    *)
      echo "Unknown argument: $1" >&2
      usage >&2
      exit 2
      ;;
  esac
done

if [[ "${METHOD}" != "auto" && "${METHOD}" != "git" && "${METHOD}" != "gh" ]]; then
  echo "Invalid --method: ${METHOD}. Use auto, git, or gh." >&2
  exit 2
fi

mkdir -p "${THIRD_PARTY_DIR}"

if [[ -d "${DEST_DIR}/.git" ]]; then
  echo "gvp-pocket_pred already exists: ${DEST_DIR}"
  echo "To update it manually: cd ${DEST_DIR} && git pull"
elif [[ -e "${DEST_DIR}" && -n "$(ls -A "${DEST_DIR}" 2>/dev/null || true)" ]]; then
  echo "Destination exists and is not empty: ${DEST_DIR}" >&2
  echo "Move it away or remove it before cloning." >&2
  exit 1
else
  rm -rf "${DEST_DIR}"
  if [[ "${METHOD}" == "gh" ]]; then
    if ! command -v gh >/dev/null 2>&1; then
      echo "GitHub CLI 'gh' not found. Install gh or use --method git." >&2
      exit 1
    fi
    echo "Cloning PocketMiner with GitHub CLI: gh repo clone ${GH_REPO} ${DEST_DIR}"
    gh repo clone "${GH_REPO}" "${DEST_DIR}"
  else
    # In auto mode, use plain git clone by default because it works without GitHub CLI authentication.
    if ! command -v git >/dev/null 2>&1; then
      echo "git not found. Install git first." >&2
      exit 1
    fi
    echo "Cloning PocketMiner with git: git clone ${REPO_URL} ${DEST_DIR}"
    git clone "${REPO_URL}" "${DEST_DIR}"
  fi
fi

mkdir -p "${DEST_DIR}/models/pocketminer"

cat <<DONE

gvp-pocket_pred setup complete.
Expected config paths:
  pocketminer_root:       third_party/gvp-pocket_pred/src
  pocketminer_checkpoint: third_party/gvp-pocket_pred/models/pocketminer

If the PocketMiner checkpoint was not included in the cloned repo, put it under:
  ${DEST_DIR}/models/pocketminer

Common commands:
  bash scripts/setup_gvp_pocket_pred.sh
  bash scripts/setup_gvp_pocket_pred.sh --method gh
  gh repo clone Mickdub/gvp third_party/gvp-pocket_pred
DONE
