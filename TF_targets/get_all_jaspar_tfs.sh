# base dir: /isdata/alab/people/pcr980/DeepCompare/
# cell: hepg2 and k562
# promoter file: base dir + Motif1_highlight_chip/summary_and_ks_test_promoters_(cell).csv
# enhancer file: base dir + Motif1_highlight_chip/summary_and_ks_test_enhancers_(cell).csv
# remove header of promoter file and enhancer file, cat them, select first column, write unique TFs to file base_dir+TF_targets/jaspar_tfs_cell.txt

base_dir="/isdata/alab/people/pcr980/DeepCompare/"
cells=("hepg2" "k562")
for c in ${cells[@]}
do
    promoter_file="Motif1_highlight_chip/summary_and_ks_test_promoters_${c}.csv"
    enhancer_file="Motif1_highlight_chip/summary_and_ks_test_enhancers_${c}.csv"
    output_file="TF_targets/jaspar_tfs_${c}.txt"
    tail -n +2 $base_dir$promoter_file | cut -d ',' -f 1 > $base_dir$promoter_file.tmp
    tail -n +2 $base_dir$enhancer_file | cut -d ',' -f 1 > $base_dir$enhancer_file.tmp
    cat $base_dir$promoter_file.tmp $base_dir$enhancer_file.tmp | sort | uniq > $base_dir$output_file
    # remove tmp files
    rm $base_dir$promoter_file.tmp
    rm $base_dir$enhancer_file.tmp
done