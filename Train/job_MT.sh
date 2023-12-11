#!/bin/bash
set -x -e


# Function to run a batch of tasks with specific parameters
run_dataset_with_multilabel() {
    local trainer=$1
    local model=$2
    ./train.py -m ${model} -d Dataset_5cv_with_multilabel_1 -f "dat" -g "0" -r 8 -t ${trainer} -n 100 -l 0.001 -b 4096 -w "[1,1,8,8,5,5,1,1]"&
    ./train.py -m ${model} -d Dataset_5cv_with_multilabel_2 -f "dat" -g "1" -r 8 -t ${trainer} -n 100 -l 0.001 -b 4096 -w "[1,1,8,8,5,5,1,1]"&
    ./train.py -m ${model} -d Dataset_5cv_with_multilabel_3 -f "dat" -g "2" -r 8 -t ${trainer} -n 100 -l 0.001 -b 4096 -w "[1,1,8,8,5,5,1,1]"&
    ./train.py -m ${model} -d Dataset_5cv_with_multilabel_4 -f "dat" -g "3" -r 8 -t ${trainer} -n 100 -l 0.001 -b 4096 -w "[1,1,8,8,5,5,1,1]"&
    ./train.py -m ${model} -d Dataset_5cv_with_multilabel_5 -f "dat" -g "4" -r 8 -t ${trainer} -n 100 -l 0.001 -b 4096 -w "[1,1,8,8,5,5,1,1]"&

    wait
}


run_dataset_without_multilabel() {
    local trainer=$1
    local model=$2
    ./train.py -m ${model} -d Dataset_5cv_without_multilabel_1 -f "dat" -g "0" -r 8 -t ${trainer} -n 100 -l 0.001 -b 4096 -w "[1,1,8,8,5,5,1,1]"&
    ./train.py -m ${model} -d Dataset_5cv_without_multilabel_2 -f "dat" -g "1" -r 8 -t ${trainer} -n 100 -l 0.001 -b 4096 -w "[1,1,8,8,5,5,1,1]"&
    ./train.py -m ${model} -d Dataset_5cv_without_multilabel_3 -f "dat" -g "2" -r 8 -t ${trainer} -n 100 -l 0.001 -b 4096 -w "[1,1,8,8,5,5,1,1]"&
    ./train.py -m ${model} -d Dataset_5cv_without_multilabel_4 -f "dat" -g "3" -r 8 -t ${trainer} -n 100 -l 0.001 -b 4096 -w "[1,1,8,8,5,5,1,1]"&
    ./train.py -m ${model} -d Dataset_5cv_without_multilabel_5 -f "dat" -g "4" -r 8 -t ${trainer} -n 100 -l 0.001 -b 4096 -w "[1,1,8,8,5,5,1,1]"&
    
    wait
}

run_regression() {
    ./train.py -m "AstigConv5D" -d Dataset_5cv_loose_1 -f "dat" -g "0" -r 8 -t "Reg_MT" -n 100 -l 0.001 -b 4096 -w "[1,1,8,8,5,5,1,1]"&
    ./train.py -m "AstigConv5D" -d Dataset_5cv_loose_2 -f "dat" -g "1" -r 8 -t "Reg_MT" -n 100 -l 0.001 -b 4096 -w "[1,1,8,8,5,5,1,1]"&
    ./train.py -m "AstigConv5D" -d Dataset_5cv_loose_3 -f "dat" -g "2" -r 8 -t "Reg_MT" -n 100 -l 0.001 -b 4096 -w "[1,1,8,8,5,5,1,1]"&
    ./train.py -m "AstigConv5D" -d Dataset_5cv_loose_4 -f "dat" -g "3" -r 8 -t "Reg_MT" -n 100 -l 0.001 -b 4096 -w "[1,1,8,8,5,5,1,1]"&
    ./train.py -m "AstigConv5D" -d Dataset_5cv_loose_5 -f "dat" -g "4" -r 8 -t "Reg_MT" -n 100 -l 0.001 -b 4096 -w "[1,1,8,8,5,5,1,1]"&
    
    wait
}



# Run all batches in sequence
for trainer in "Reg_MT" "CR_MT"; do
    if [ "$trainer" == "CR_MT" ]; then
        model="AstigCRConv5D"
    else
        model="AstigConv5D"
    fi
    run_dataset_with_multilabel "$trainer" "$model"
    run_dataset_without_multilabel "$trainer" "$model"
    run_regression
done



