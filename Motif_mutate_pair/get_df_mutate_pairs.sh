python3 get_df_mutate_pairs.py \
  --file_prefix "promoters_hepg2" \
  --remap_cell_type "Hep-G2" \
  --track_num 14 \
  --device 0 &

python3 get_df_mutate_pairs.py \
  --file_prefix "promoters_k562" \
  --remap_cell_type "K-562" \
  --track_num 15 \
  --device 1 &

python3 get_df_mutate_pairs.py \
  --file_prefix "enhancers_hepg2" \
  --remap_cell_type "Hep-G2" \
  --track_num 14 \
  --device 2 &

python3 get_df_mutate_pairs.py \
  --file_prefix "enhancers_k562" \
  --remap_cell_type "K-562" \
  --track_num 15 \
  --device 3 &




# nohup bash get_df_mutate_pairs.sh > get_df_mutate_pairs.out &