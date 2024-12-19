

files=(dhs_constrained_distal_ti_hepg2.tsv dhs_constrained_distal_ti_k562.tsv dhs_constrained_distal_ts_hepg2.tsv dhs_constrained_distal_ts_k562.tsv dhs_constrained_proximal_ti_hepg2.tsv dhs_constrained_proximal_ti_k562.tsv dhs_constrained_proximal_ts_hepg2.tsv dhs_constrained_proximal_ts_k562.tsv dhs_nonconstrained_distal_ti_hepg2.tsv dhs_nonconstrained_distal_ti_k562.tsv dhs_nonconstrained_distal_ts_hepg2.tsv dhs_nonconstrained_distal_ts_k562.tsv dhs_nonconstrained_proximal_ti_hepg2.tsv dhs_nonconstrained_proximal_ti_k562.tsv dhs_nonconstrained_proximal_ts_hepg2.tsv dhs_nonconstrained_proximal_ts_k562.tsv)

for file in ${files[@]}
do
    echo $file
    python3 get_motif_df.py --file_name $file &
done


# nohup bash get_motif_df.sh > get_motif_df.out &