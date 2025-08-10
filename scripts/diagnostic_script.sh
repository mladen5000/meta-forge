#!/bin/bash

echo "=== COMPREHENSIVE RUST PROJECT DIAGNOSTICS ==="
echo "Project: MetaForge Metagenomic Pipeline"
echo "Date: $(date)"
echo "============================================="

echo ""
echo "1. CARGO CHECK - Compilation Errors"
echo "------------------------------------"
cargo check --all-targets --all-features 2>&1

echo ""
echo "2. CARGO CLIPPY - Code Quality Issues"
echo "------------------------------------"
cargo clippy --all-targets --all-features # -- -D warnings 2>&1

echo ""
echo "3. CARGO TEST - Test Failures"
echo "------------------------------"
# cargo test --all-targets --all-features 2>&1

echo ""
echo "4. CARGO BUILD - Full Build Attempt"
echo "-----------------------------------"
# cargo build --all-targets --all-features 2>&1

echo ""
echo "5. DEPENDENCY TREE CHECK"
echo "------------------------"
cargo tree 2>&1

echo ""
echo "=== DIAGNOSTICS COMPLETE ==="