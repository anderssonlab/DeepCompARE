# python3 get_motif_df.py --file_name promoters_hepg2 --device 0 &
# python3 get_motif_df.py --file_name enhancers_hepg2 --device 1 &
# python3 get_motif_df.py --file_name promoters_k562 --device 2 &
# python3 get_motif_df.py --file_name enhancers_k562 --device 3 &





python3 get_motif_df.py --file_name dhs_distal_hepg2  --device 2 &
python3 get_motif_df.py --file_name dhs_distal_k562   --device 3 &
python3 get_motif_df.py --file_name dhs_proximal_hepg2  --device 6 &
python3 get_motif_df.py --file_name dhs_proximal_k562   --device 7 &















# nohup bash get_motif_df.sh > get_motif_df.out &