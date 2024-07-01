python3 get_motif_df.py --file_name promoters_hepg2 --track_num 8 --device 0 &
python3 get_motif_df.py --file_name promoters_k562 --track_num 9 --device 1 &
python3 get_motif_df.py --file_name promoters_hepg2 --track_num 10 --device 2 &
python3 get_motif_df.py --file_name promoters_k562 --track_num 11 --device 3 &

python3 get_motif_df.py --file_name promoters_hepg2 --track_num 12 --device 4 &
python3 get_motif_df.py --file_name promoters_k562 --track_num 13 --device 5 &
python3 get_motif_df.py --file_name promoters_hepg2 --track_num 14 --device 6 &
python3 get_motif_df.py --file_name promoters_k562 --track_num 15 --device 7 &


wait

python3 get_motif_df.py --file_name enhancers_hepg2 --track_num 8 --device 0 &
python3 get_motif_df.py --file_name enhancers_k562 --track_num 9 --device 1 &
python3 get_motif_df.py --file_name enhancers_hepg2 --track_num 10 --device 2 &
python3 get_motif_df.py --file_name enhancers_k562 --track_num 11 --device 3 &

python3 get_motif_df.py --file_name enhancers_hepg2 --track_num 12 --device 4 &
python3 get_motif_df.py --file_name enhancers_k562 --track_num 13 --device 5 &
python3 get_motif_df.py --file_name enhancers_hepg2 --track_num 14 --device 6 &
python3 get_motif_df.py --file_name enhancers_k562 --track_num 15 --device 7 &





#  nohup bash get_motif_df.sh > get_motif_df.out &