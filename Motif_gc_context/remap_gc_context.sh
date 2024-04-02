#!/bin/bash

# Path to the file containing TFs
FILE="/isdata/alab/people/pcr980/DeepCompare/TF_targets/tf_targets_k562.txt"

# Command to run in parallel
run_command() {
    TF=$1
    python3 remap_gc_context.py --tf "$TF"
}

export -f run_command

# Using parallel to process 80 TFs at a time
cat "$FILE" | tail -80 | parallel -j 80 run_command
