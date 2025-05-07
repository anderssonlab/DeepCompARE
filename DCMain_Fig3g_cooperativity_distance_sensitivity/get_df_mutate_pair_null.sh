python3 get_df_mutate_pair_null.py \
  --file_prefix "promoters_hepg2" \
  --device 2 &

python3 get_df_mutate_pair_null.py \
  --file_prefix "enhancers_hepg2" \
  --device 3 &

python3 get_df_mutate_pair_null.py \
  --file_prefix "promoters_k562" \
  --device 4 &

python3 get_df_mutate_pair_null.py \
  --file_prefix "enhancers_k562" \
  --device 5 &






# nohup bash get_df_mutate_pair_null.sh > get_df_mutate_pair_null.out &