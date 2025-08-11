#!/bin/bash
# Compilation check script

echo "🔍 Running compilation check across the entire codebase..."
echo "================================================"

# Run cargo check to identify compilation errors
echo "1. Running cargo check..."
cargo check --all-targets --all-features 2>&1 | tee cargo_check.log

echo ""
echo "2. Running cargo clippy for additional issues..."
cargo clippy --all-targets --all-features -- -D warnings 2>&1 | tee cargo_clippy.log

echo ""
echo "3. Running cargo test --no-run to check test compilation..."
cargo test --no-run 2>&1 | tee cargo_test_compile.log

echo ""
echo "4. Summary of issues found:"
echo "=========================="

# Extract compilation errors
echo "Compilation errors from cargo check:"
grep -E "error\[E[0-9]+\]|error:" cargo_check.log || echo "No compilation errors found in cargo check"

echo ""
echo "Compilation errors from cargo clippy:"
grep -E "error\[E[0-9]+\]|error:" cargo_clippy.log || echo "No compilation errors found in cargo clippy"

echo ""
echo "Test compilation errors:"
grep -E "error\[E[0-9]+\]|error:" cargo_test_compile.log || echo "No test compilation errors found"

echo ""
echo "Check complete. See logs for detailed output."