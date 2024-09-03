python3 get_motif_df.py --file_name promoters_hepg2 --device 2 &
python3 get_motif_df.py --file_name enhancers_hepg2 --device 2 &
python3 get_motif_df.py --file_name promoters_k562 --device 2 &
python3 get_motif_df.py --file_name enhancers_k562 --device 2 &


# python3 get_motif_df.py --file_name resize_600bp_CAGE_HepG2  --device 4 &



# nohup bash get_motif_df.sh > get_motif_df.out &