python3 get_df_mutate_pairs.py \
  --file_name "promoters_hepg2" \
  --device 0 &

python3 get_df_mutate_pairs.py \
  --file_name "enhancers_hepg2" \
  --device 1 &

python3 get_df_mutate_pairs.py \
  --file_name "promoters_k562" \
  --device 2 &

python3 get_df_mutate_pairs.py \
  --file_name "enhancers_k562" \
  --device 3 &




# nohup bash get_df_mutate_pairs.sh > get_df_mutate_pairs.out &