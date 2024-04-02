#python3 write_maf_and_effect_size.py --range_bed_file promoters_hepg2 
#python3 write_maf_and_effect_size.py --range_bed_file promoters_k562 
#python3 write_maf_and_effect_size.py --range_bed_file enhancers_k562 
python3 write_maf_and_effect_size.py --range_bed_file enhancers_hepg2
echo "Done"

# nohup bash write_maf_and_effect_size.sh > write_maf_and_effect_size.out &
