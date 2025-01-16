python3 get_predictions.py --file_name promoters_hepg2 --device 0 &
python3 get_predictions.py --file_name enhancers_hepg2 --device 1 &
python3 get_predictions.py --file_name promoters_k562 --device 2 &
python3 get_predictions.py --file_name enhancers_k562 --device 3 &

# nohup bash get_predictions.sh > get_predictions.out &