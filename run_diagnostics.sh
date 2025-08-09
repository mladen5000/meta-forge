#!/bin/bash
cd /Users/mladenrasic/Projects/metagenomic_llm
chmod +x diagnostic_script.sh
./diagnostic_script.sh > diagnostic_output.log 2>&1
echo "Diagnostics complete. Check diagnostic_output.log for results."
cat diagnostic_output.log