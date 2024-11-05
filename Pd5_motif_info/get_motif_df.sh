# python3 get_motif_df.py --file_name enhancers_hepg2 --device 2 &
# python3 get_motif_df.py --file_name promoters_k562 --device 2 &
# python3 get_motif_df.py --file_name enhancers_k562 --device 2 &


# python3 get_motif_df.py --file_name resize_600bp_CAGE_HepG2  --device 4 &

python3 get_motif_df.py --file_name DNase_E118_hepg2  --device 6 &
python3 get_motif_df.py --file_name DNase_E123_k562   --device 7 &


# nohup bash get_motif_df.sh > get_motif_df.out &