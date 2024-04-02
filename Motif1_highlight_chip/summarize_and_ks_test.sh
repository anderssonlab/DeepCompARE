#!/bin/bash

set -e

python3 b_summarize_and_ks_test.py --input "enhancers_hepg2" &
python3 b_summarize_and_ks_test.py --input "enhancers_k562" &
python3 b_summarize_and_ks_test.py --input "promoters_hepg2" &
python3 b_summarize_and_ks_test.py --input "promoters_k562" &

# nohup bash summarize_and_ks_test.sh > summarize_and_ks_test.out &