# dhs_constrained_distal_ti_k562.tsv     dhs_nonconstrained_distal_ti_k562.tsv
# dhs_constrained_distal_ts_k562.tsv     dhs_nonconstrained_distal_ts_k562.tsv
# dhs_constrained_proximal_ti_k562.tsv   dhs_nonconstrained_proximal_ti_k562.tsv
# dhs_constrained_proximal_ts_k562.tsv   dhs_nonconstrained_proximal_ts_k562.tsv



python3 get_df_mutate_pairs.py \
  --file_name dhs_constrained_distal_ti_k562 \
  --device 0 &

python3 get_df_mutate_pairs.py \
  --file_name dhs_nonconstrained_distal_ti_k562 \
  --device 1 &

python3 get_df_mutate_pairs.py \
  --file_name dhs_constrained_distal_ts_k562 \
  --device 2 &

python3 get_df_mutate_pairs.py \
  --file_name dhs_nonconstrained_distal_ts_k562 \
  --device 4 &

python3 get_df_mutate_pairs.py \
  --file_name dhs_constrained_proximal_ti_k562 \
  --device 5 &

python3 get_df_mutate_pairs.py \
  --file_name dhs_nonconstrained_proximal_ti_k562 \
  --device 5 &

python3 get_df_mutate_pairs.py \
  --file_name dhs_constrained_proximal_ts_k562 \
  --device 7 &

python3 get_df_mutate_pairs.py \
  --file_name dhs_nonconstrained_proximal_ts_k562 \
  --device 7 &













# nohup bash get_df_mutate_pairs.sh > get_df_mutate_pairs.tits_k562.out 2>&1 &