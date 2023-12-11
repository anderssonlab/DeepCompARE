./train.py -m "BPNet" -d Dataset_final_rep -f "cage_hepg2" \
           -g "4" -r 1 -t "Class_ST" -n 100 -l 0.001 -b 4096 &

./train.py -m "BPNet" -d Dataset_final_rep -f "dhs_hepg2" \
           -g "5" -r 1 -t "Class_ST" -n 100 -l 0.001 -b 4096 &

./train.py -m "BPNet" -d Dataset_final_rep -f "starr_hepg2" \
           -g "6" -r 1 -t "Class_ST" -n 100 -l 0.001 -b 4096 &

./train.py -m "BPNet" -d Dataset_final_rep -f "sure_hepg2" \
           -g "7" -r 1 -t "Class_ST" -n 100 -l 0.001 -b 4096 &
