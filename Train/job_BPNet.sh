#!/bin/bash
set -x -e


# Function to run a batch of tasks with specific parameters
run_batch() {
    local assay=$1
    local cell=$2
    ./train.py -m "BPNet21" -d Dataset_final_rep -f "${assay}_${cell}" -g "5" -r 1 -t "Reg_ST" -n 100 -l 0.001 -b 4096 &
    ./train.py -m "BPNet21" -d Dataset_final_rep -f "${assay}_${cell}" -g "5" -r 1 -t "Class_ST" -n 100 -l 0.001 -b 4096 &
    ./train.py -m "BPNet64" -d Dataset_final_rep -f "${assay}_${cell}" -g "6" -r 1 -t "Reg_ST" -n 100 -l 0.001 -b 4096 &
    ./train.py -m "BPNet64" -d Dataset_final_rep -f "${assay}_${cell}" -g "7" -r 1 -t "Class_ST" -n 100 -l 0.001 -b 4096 &

    wait
}

# Run all batches in sequence
for assay in "cage" "dhs" "starr" "sure"; do
    for cell in "hepg2" "k562"; do
        run_batch "$assay" "$cell"
    done
done