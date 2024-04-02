python3 get_motif_df.py --file_name promoters_hepg2 --track_num 0
python3 get_motif_df.py --file_name promoters_k562 --track_num 1
python3 get_motif_df.py --file_name enhancers_hepg2 --track_num 4
python3 get_motif_df.py --file_name enhancers_k562 --track_num 5
echo "Done"

#  nohup bash get_motif_df.sh > get_motif_df.out &