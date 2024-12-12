# python3 get_df_mutate_pairs.py \
#   --file_name dhs_constrained_distal_ti_hepg2 \
#   --device 6 &

# python3 get_df_mutate_pairs.py \
#   --file_name dhs_nonconstrained_distal_ti_hepg2 \
#   --device 7 &

python3 get_df_mutate_pairs.py \
  --file_name dhs_constrained_distal_ts_hepg2 \
  --device 0 &

# python3 get_df_mutate_pairs.py \
#   --file_name dhs_nonconstrained_distal_ts_hepg2 \
#   --device 7 



python3 get_df_mutate_pairs.py \
  --file_name dhs_constrained_proximal_ti_hepg2 \
  --device 0 &

python3 get_df_mutate_pairs.py \
  --file_name dhs_nonconstrained_proximal_ti_hepg2 \
  --device 1 &

python3 get_df_mutate_pairs.py \
  --file_name dhs_constrained_proximal_ts_hepg2 \
  --device 1 &

python3 get_df_mutate_pairs.py \
  --file_name dhs_nonconstrained_proximal_ts_hepg2 \
  --device 2 & 










# nohup bash get_df_mutate_pairs.sh > get_df_mutate_pairs.tits_hepg2.out 2>&1 &