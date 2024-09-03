# python3 get_df_mutate_pairs.py \
#   --file_prefix "promoters_hepg2" \
#   --device 4 &

# python3 get_df_mutate_pairs.py \
#   --file_prefix "enhancers_hepg2" \
#   --device 5 &

# python3 get_df_mutate_pairs.py \
#   --file_prefix "promoters_k562" \
#   --device 6 &

# python3 get_df_mutate_pairs.py \
#   --file_prefix "enhancers_k562" \
#   --device 7 &



python3 get_df_mutate_pairs.py \
  --file_prefix "promoters_broad_common" \
  --device 4 &

python3 get_df_mutate_pairs.py \
  --file_prefix "promoters_broad_hepg2" \
  --device 5 &

python3 get_df_mutate_pairs.py \
  --file_prefix "promoters_broad_k562" \
  --device 6 &

python3 get_df_mutate_pairs.py \
  --file_prefix "promoters_sharp_common" \
  --device 4 &

python3 get_df_mutate_pairs.py \
  --file_prefix "promoters_sharp_hepg2" \
  --device 5 &

python3 get_df_mutate_pairs.py \
  --file_prefix "promoters_sharp_k562" \
  --device 6 &



# nohup bash get_df_mutate_pairs.sh > get_df_mutate_pairs.out &