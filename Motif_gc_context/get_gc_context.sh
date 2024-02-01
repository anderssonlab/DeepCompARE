#!/bin/bash

# Path to the file containing TFs
FILE="/isdata/alab/people/pcr980/DeepCompare/Putative_binding_vs_actual_binding/tf_targets_k562.txt"

# Command to run in parallel
run_command() {
    TF=$1
    python3 get_gc_context.py --tf "$TF"
}

export -f run_command

# Using parallel to process 20 TFs at a time
cat "$FILE" | parallel -j 40 run_command
