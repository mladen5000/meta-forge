#!/bin/bash
cd /Users/mladenrasic/Projects/metagenomic_llm
echo "Checking Rust project compilation status..."
timeout 60 cargo check 2>&1 || echo "Check timed out or failed"