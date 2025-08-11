#!/bin/bash

echo "🔍 Quick Compilation Analysis"
echo "============================"

# Simple cargo check to see immediate issues
cargo check --all-targets 2>&1 | head -50