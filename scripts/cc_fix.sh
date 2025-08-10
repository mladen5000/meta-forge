#!/usr/bin/env bash
set -euo pipefail

CLAUDE_CMD="${CLAUDE_CMD:-claude}"

# Run tests to capture failures
cargo test || true
cargo test -- --format json 2>/dev/null | tail -n 400 > .last-cargo-test.json || true

TMP_DIR="$(mktemp -d)"
BUNDLE="$TMP_DIR/fix-context.md"

{
  echo "# Minimal changes to make tests pass"
  echo "## PYTEST OUTPUT (last 400 lines)"
  echo '```'
  cat .last-cargo-test.json
  echo '```'
  echo
  echo "## Changed files"
  git diff --name-only | while read -r f; do
    [ -f "$f" ] && { echo "### FILE: $f"; echo '```'; sed -n '1,400p' "$f"; echo '```'; }
  done
} > "$BUNDLE"

"$CLAUDE_CMD" code --input "$BUNDLE" --apply=true

scripts/rs_preflight.sh test --full
