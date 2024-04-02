python3 get_df_mutate_pairs.py \
  --file_prefix "promoters_broad_common" \
  --remap_cell_type "Hep-G2" \
  --track_num 0 \
  --device 4 &

python3 get_df_mutate_pairs.py \
  --file_prefix "promoters_broad_common" \
  --remap_cell_type "K-562" \
  --track_num 1 \
  --device 5 &


python3 get_df_mutate_pairs.py \
  --file_prefix "promoters_sharp_common" \
  --remap_cell_type "Hep-G2" \
  --track_num 0 \
  --device 6 &

python3 get_df_mutate_pairs.py \
  --file_prefix "promoters_sharp_common" \
  --remap_cell_type "K-562" \
  --track_num 1 \
  --device 7 &


# nohup bash get_df_mutate_pairs.sh > get_df_mutate_pairs.out &