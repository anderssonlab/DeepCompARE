# python3 get_motif_df.py --file_name promoters_hepg2 --track_num 0 --device 0 &
python3 get_motif_df.py --file_name resize_600bp_CAGE_HepG2 --track_num 0 --device 4 &
python3 get_motif_df.py --file_name resize_600bp_CAGE_K562 --track_num 1 --device 5 &
python3 get_motif_df.py --file_name resize_600bp_DHS_HepG2 --track_num 2 --device 6 &
python3 get_motif_df.py --file_name resize_600bp_DHS_K562 --track_num 3 --device 7 &

wait

python3 get_motif_df.py --file_name resize_600bp_STARR_HepG2 --track_num 4 --device 4 &
python3 get_motif_df.py --file_name resize_600bp_STARR_K562 --track_num 5 --device 5 &
python3 get_motif_df.py --file_name resize_600bp_SuRE_HepG2 --track_num 6 --device 6 &
python3 get_motif_df.py --file_name resize_600bp_SuRE_K562 --track_num 7 --device 7 &

wait

python3 get_motif_df.py --file_name resize_600bp_CAGE_HepG2 --track_num 8 --device 4 &
python3 get_motif_df.py --file_name resize_600bp_CAGE_K562 --track_num 9 --device 5 &
python3 get_motif_df.py --file_name resize_600bp_DHS_HepG2 --track_num 10 --device 6 &
python3 get_motif_df.py --file_name resize_600bp_DHS_K562 --track_num 11 --device 7 &

wait

python3 get_motif_df.py --file_name resize_600bp_STARR_HepG2 --track_num 12 --device 4 &
python3 get_motif_df.py --file_name resize_600bp_STARR_K562 --track_num 13 --device 5 &
python3 get_motif_df.py --file_name resize_600bp_SuRE_HepG2 --track_num 14 --device 6 &
python3 get_motif_df.py --file_name resize_600bp_SuRE_K562 --track_num 15 --device 7 &

#  nohup bash get_motif_df.sh > get_motif_df.out &