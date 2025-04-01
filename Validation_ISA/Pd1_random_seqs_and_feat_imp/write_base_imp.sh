
python3 /isdata/alab/people/pcr980/Scripts_python/write_base_importance.py \
  --file_name random_seqs.csv \
  --format seq \
  --imp_type ism \
  --out_path df_ism.csv \
  --device 1 &


python3 /isdata/alab/people/pcr980/Scripts_python/write_base_importance.py \
  --file_name random_seqs.csv \
  --format seq \
  --imp_type isa \
  --out_path df_isa.csv \
  --device 5 &



python3 /isdata/alab/people/pcr980/Scripts_python/write_base_importance.py \
  --file_name random_seqs.csv \
  --format seq \
  --imp_type gradxinp \
  --out_path df_gradxinp.csv \
  --device 6 

echo "Done!"


# nohup bash write_base_imp.sh > write_base_imp.out &