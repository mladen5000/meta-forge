#!/usr/bin/env bash
set -euo pipefail

CLAUDE_CMD="${CLAUDE_CMD:-claude}"
PROMPT_FILE="${1:-.prompts/refactor_surgical.md}"

TMP_DIR="$(mktemp -d)"
BUNDLE="$TMP_DIR/diff-request.md"

{
  echo "# Output unified diff only. No extra text."
  echo
  cat "$PROMPT_FILE"
  echo
  echo "## Project snapshot"
  for f in Cargo.toml rust-toolchain.toml README.md; do
    [ -f "$f" ] && { echo "### FILE: $f"; echo '```'; sed -n '1,400p' "$f"; echo '```'; }
  done
} > "$BUNDLE"

"$CLAUDE_CMD" code --input "$BUNDLE" --format=diff --apply=false
