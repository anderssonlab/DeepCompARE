#!/bin/bash

set -e

python3 b_summary_ks_test.py --track cage &
python3 b_summary_ks_test.py --track dhs &
python3 b_summary_ks_test.py --track starr &
python3 b_summary_ks_test.py --track sure &

# nohup bash b_summary_ks_test.sh > b_summary_ks_test.out &