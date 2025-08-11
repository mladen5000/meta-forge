#!/usr/bin/env bash
set -euo pipefail

# Rust project preflight checks for Claude Code hooks
# Anthropic best-practices: keep this deterministic, local, and safe

cmd="${1:-}"
mode="${2:---smoke}"

run() {
  echo "▶ $*"
  local attempt=1
  until eval "$@"; do
    if [[ $attempt -ge 3 ]]; then
      echo "❌ Command failed after $attempt attempts"
      return 1
    fi
    echo "⏳ Cargo lock busy, retrying in 3s (attempt $attempt)..."
    sleep 3
    ((attempt++))
  done
}

have_nextest() {
  command -v cargo-nextest >/dev/null 2>&1
}

case "$cmd" in
  fmt)
    run "cargo fmt --all"
    ;;
  lint)
    run "cargo clippy --all-targets --all-features # -- -D warnings"
    ;;
  check)
    run "cargo check --all-targets --all-features"
    ;;
  test)
    if [ "$mode" = "--smoke" ]; then
      if have_nextest; then
        run "cargo nextest run --skip slow --failure-output=final"
      else
        run "cargo test -- --quiet --skip slow"
      fi
    else
      if have_nextest; then
        run "cargo nextest run --failure-output=final"
      else
        run "cargo test"
      fi
    fi
  ;;
  *)
    echo "Usage: $0 {fmt|lint|check|test [--smoke|--full]}"
    exit 2
    ;;
esac
