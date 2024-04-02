#!/bin/bash



# Define the base path more concisely
base_path="/isdata/alab/people/pcr980/DeepCompare/Pd4_promoters_enhancers_and_featimp"

# List of input files and their corresponding device IDs
files=("enhancers_hepg2")
devices=(0) # Assuming device IDs are still relevant and sequential

# Loop over files array
for i in "${!files[@]}"; do
    # Construct input and output file paths
    input_file="${files[$i]}.bed"
    output_file="gradxinp_${input_file}.csv"

    # Execute the command with the constructed file paths and device ID
    python3 /isdata/alab/people/pcr980/DeepCompare/Scripts_python/write_gradxinp_bed.py \
        --bed_file "$base_path/$input_file" \
        --out_path "$base_path/$output_file" \
        --device "${devices[$i]}" 
done

echo "Done"


# nohup bash write_gradxinp_from_bed.sh > write_gradxinp_from_bed.out &