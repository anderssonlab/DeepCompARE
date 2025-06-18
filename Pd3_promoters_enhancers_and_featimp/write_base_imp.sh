
python3 /isdata/alab/people/pcr980/Scripts_python/write_base_importance.py \
  --file_name /isdata/alab/people/pcr980/DeepCompare/Pd1_bed_processed/CAGE_K562.bed
  --imp_type gradxinp \
  --out_path gradxinp_cage_hepg2.csv \
  --device 0 &

python3 /isdata/alab/people/pcr980/Scripts_python/write_base_importance.py \
  --file_name /isdata/alab/people/pcr980/DeepCompare/Pd1_bed_processed/CAGE_K562.bed
  --imp_type gradxinp \
  --out_path gradxinp_cage_k562.csv \
  --device 1 &







# nohup bash write_base_imp.sh > write_base_imp.out &