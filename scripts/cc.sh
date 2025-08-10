#!/usr/bin/env bash
set -euo pipefail

CLAUDE_CMD="${CLAUDE_CMD:-claude}"
PROMPT_FILE="${1:-}"
shift || true

[ -n "$PROMPT_FILE" ] || { echo "Usage: $0 .prompts/xyz.md"; exit 2; }

# Preflight before passing context to Claude
scripts/rs_preflight.sh fmt
scripts/rs_preflight.sh lint
scripts/rs_preflight.sh check
scripts/rs_preflight.sh test --smoke

TMP_DIR="$(mktemp -d)"
BUNDLE="$TMP_DIR/context.md"

{
  echo "# Prompt"
  cat "$PROMPT_FILE"
  echo
  echo "# Key project files"
  for f in Cargo.toml rust-toolchain.toml README.md; do
    [ -f "$f" ] && { echo "## FILE: $f"; echo '```'; sed -n '1,400p' "$f"; echo '```'; }
  done
} > "$BUNDLE"

"$CLAUDE_CMD" code --input "$BUNDLE" --apply=true

scripts/rs_preflight.sh test --full
