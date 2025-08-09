#!/bin/bash
set -e
cd /Users/mladenrasic/Projects/metagenomic_llm

echo "=== BATCH COMPILATION AND TEST CHECK ==="

# Check basic compilation
echo "1. Basic cargo check..."
cargo check 2>&1 | tee check_output.log

echo ""
echo "2. Test compilation check..."
cargo test --no-run 2>&1 | tee test_compile_output.log

echo ""
echo "3. Running tests..."
cargo test 2>&1 | tee test_run_output.log

echo ""
echo "4. Binary compilation check..."
cargo build --bins 2>&1 | tee binary_output.log

echo "=== BATCH CHECK COMPLETE ==="