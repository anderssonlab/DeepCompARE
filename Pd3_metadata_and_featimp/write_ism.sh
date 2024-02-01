#!/bin/bash

# Define the base paths for easier maintenance and readability
script_path="/isdata/alab/people/pcr980/DeepCompare/Scripts_python/ism_N.py"
seq_base_path="/isdata/alab/people/pcr980/DeepCompare/Pd2_sequences"
out_base_path="/isdata/alab/people/pcr980/DeepCompare/Pd3_metadata_and_featimp"

# Array of sequence types
seq_types=("CAGE" "DHS" "STARR" "SuRE")

# Loop through each sequence type for K562
for i in "${!seq_types[@]}"; do
    seq_file_k562="${seq_base_path}/seqs_${seq_types[$i]}_K562.csv"
    out_path_k562="${out_base_path}/ism_${seq_types[$i]}_K562.csv"
    device_k562="$i"  # GPUs 0, 1, 2, 3

    # Execute the script for K562
    "${script_path}" \
        --seq_file "${seq_file_k562}" \
        --seq_colname sequence \
        --out_path "${out_path_k562}" \
        --gpu "${device_k562}" &
done

# Loop through each sequence type for HepG2
for i in "${!seq_types[@]}"; do
    seq_file_hepg2="${seq_base_path}/seqs_${seq_types[$i]}_HepG2.csv"
    out_path_hepg2="${out_base_path}/ism_${seq_types[$i]}_HepG2.csv"
    device_hepg2=$((i + 4))  # GPUs 4, 5, 6, 7

    # Execute the script for HepG2
    "${script_path}" \
        --seq_file "${seq_file_hepg2}" \
        --seq_colname sequence \
        --out_path "${out_path_hepg2}" \
        --gpu "${device_hepg2}" &
done
