
# python3 /isdata/alab/people/pcr980/DeepCompare/Scripts_python/write_ism_N_bed.py \
#   --bed_file promoters_hepg2.bed \
#   --out_path ism_promoters_hepg2.csv \
#   --device 0 &

# python3 /isdata/alab/people/pcr980/DeepCompare/Scripts_python/write_ism_N_bed.py \
#   --bed_file promoters_k562.bed \
#   --out_path ism_promoters_k562.csv \
#   --device 1 &

# python3 /isdata/alab/people/pcr980/DeepCompare/Scripts_python/write_ism_N_bed.py \
#   --bed_file enhancers_k562.bed \
#   --out_path ism_enhancers_k562.csv \
#   --device 2 &



python3 /isdata/alab/people/pcr980/DeepCompare/Scripts_python/write_ism_N_bed.py \
  --bed_file enhancers_hepg2.bed \
  --out_path ism_enhancers_hepg2.csv \
  --device 2 

echo "Done"


# nohup bash write_ism.sh > write_ism.out &