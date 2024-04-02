#!/bin/bash

# conditions=("CAGE" "DHS" "STARR" "SuRE")
# cell_lines=("HepG2" "K562")

# for condition in "${conditions[@]}"; do
#     for cell_line in "${cell_lines[@]}"; do
#         /isdata/alab/people/pcr980/DeepCompare/Scripts_python/write_gradxinp_seq.py \
#             --seq_file /isdata/alab/people/pcr980/DeepCompare/Pd2_sequences/seqs_${condition}_${cell_line}.csv \
#             --seq_colname sequence \
#             --out_path /isdata/alab/people/pcr980/DeepCompare/Pd3_metadata_and_featimp/gradxinp_${condition}_${cell_line}.csv \
#             --device 1
#     done
# done

# echo Done




/isdata/alab/people/pcr980/DeepCompare/Scripts_python/write_gradxinp_seq.py \
  --seq_file /isdata/alab/people/pcr980/DeepCompare/Pd3_metadata_and_featimp/metadata_all_seqs.csv \
  --seq_colname seq \
  --out_path /isdata/alab/people/pcr980/DeepCompare/Pd3_metadata_and_featimp/gradxinp_all_seqs.csv 