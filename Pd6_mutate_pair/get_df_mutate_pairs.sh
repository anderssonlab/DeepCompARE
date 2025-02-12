# python3 get_df_mutate_pairs.py \
#   --file_name "promoters_hepg2" \
#   --device 3 &

# python3 get_df_mutate_pairs.py \
#   --file_name "enhancers_hepg2" \
#   --device 4 &

# python3 get_df_mutate_pairs.py \
#   --file_name "promoters_k562" \
#   --device 6 &

# python3 get_df_mutate_pairs.py \
#   --file_name "enhancers_k562" \
#   --device 6 &





python3 get_df_mutate_pairs.py --file_name "dhs_distal_hepg2" --device 0 &
python3 get_df_mutate_pairs.py --file_name "dhs_distal_k562" --device 1 &
python3 get_df_mutate_pairs.py --file_name "dhs_proximal_hepg2" --device 4 &
python3 get_df_mutate_pairs.py --file_name "dhs_proximal_k562" --device 5 &



# nohup bash get_df_mutate_pairs.sh > get_df_mutate_pairs.out &