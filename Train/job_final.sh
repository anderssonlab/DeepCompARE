#!/bin/bash
set -x -e

# ./train.py -m AstigCRConv5D -d Dataset_final_rep -f "dat" -g "3" -r 8 -t "CR_MT" -n 75 -l 0.001 -b 4096 -w "[1,1,8,8,5,5,1,1]"&
./train.py -m AstigCRConv6D -d Dataset_final_rep -f "dat" -g "7" -r 8 -t "CR_MT" -n 75 -l 0.001 -b 4096 -w "[1,1,8,8,5,5,1,1]"&


# nohup bash job_final.sh > job_final.log 2>&1 &