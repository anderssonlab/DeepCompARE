
python3 /isdata/alab/people/pcr980/Scripts_python/write_base_importance.py \
  --bed_file promoters_hepg2.bed \
  --imp_type ism \
  --out_path ism_promoters_hepg2.csv \
  --device 1 &

python3 /isdata/alab/people/pcr980/Scripts_python/write_base_importance.py \
  --bed_file promoters_k562.bed \
  --imp_type ism \
  --out_path ism_promoters_k562.csv \
  --device 5 &

python3 /isdata/alab/people/pcr980/Scripts_python/write_base_importance.py \
  --bed_file enhancers_k562.bed \
  --imp_type ism \
  --out_path ism_enhancers_k562.csv \
  --device 6 &

python3 /isdata/alab/people/pcr980/Scripts_python/write_base_importance.py \
  --bed_file enhancers_hepg2.bed \
  --imp_type ism \
  --out_path ism_enhancers_hepg2.csv \
  --device 7 

echo "Done"


# nohup bash write_base_imp.sh > write_base_imp.out &