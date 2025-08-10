#!/bin/bash
cd /Users/mladenrasic/Projects/metagenomic_llm

echo "=== DIRECT CARGO CHECK OUTPUT ==="
cargo check --color=never 2>&1 | head -100

echo ""
echo "=== COUNTING ERROR TYPES ==="
cargo check --color=never 2>&1 | grep -c "error:" || echo "No 'error:' found"
cargo check --color=never 2>&1 | grep -c "warning:" || echo "No 'warning:' found"

echo ""
echo "=== FIRST 10 ERRORS ==="
cargo check --color=never 2>&1 | grep -A2 "error:" | head -30