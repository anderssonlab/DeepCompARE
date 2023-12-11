#!/bin/bash
set -x -e


# Function to run a batch of tasks with specific parameters
run_batch() {
    local feature=$1

    ./train.py -m "Conv5D" -d Dataset_final_rep -f "${feature}_hepg2" -g "0" -r 1 -t "Reg_ST" -n 100 -l 0.001 -b 4096 &
    ./train.py -m "Conv5D" -d Dataset_final_rep -f "${feature}_hepg2" -g "1" -r 1 -t "Class_ST" -n 100 -l 0.001 -b 4096 &
    ./train.py -m "CRConv5D" -d Dataset_final_rep -f "${feature}_hepg2" -g "2" -r 1 -t "CR_ST" -n 100 -l 0.001 -b 4096 &

    ./train.py -m "Conv5D" -d Dataset_final_rep -f "${feature}_k562" -g "3" -r 1 -t "Reg_ST" -n 100 -l 0.001 -b 4096 &
    ./train.py -m "Conv5D" -d Dataset_final_rep -f "${feature}_k562" -g "4" -r 1 -t "Class_ST" -n 100 -l 0.001 -b 4096 &
    ./train.py -m "CRConv5D" -d Dataset_final_rep -f "${feature}_k562" -g "5" -r 1 -t "CR_ST" -n 100 -l 0.001 -b 4096 &

    # Wait for all background processes in this batch to finish
    wait
}

# Run all batches in sequence
for feature in "cage" "dhs" "starr" "sure"; do
    run_batch "$feature"
done
