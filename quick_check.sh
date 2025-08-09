#!/bin/bash
cd /Users/mladenrasic/Projects/metagenomic_llm

echo "=== QUICK CARGO CHECK ==="
cargo check --message-format=json 2>&1 | grep -E '"message"|"code"|"level":"error"' | head -50

echo ""
echo "=== QUICK TEST CHECK ==="
cargo test --no-run --message-format=json 2>&1 | grep -E '"message"|"code"|"level":"error"' | head -20