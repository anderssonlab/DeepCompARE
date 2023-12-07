# run this in directory Quantify_regions/

dir_base="/isdata/alab/people/pcr980/DeepCompare/Datasets/Dataset_final_rep/"

# Define an array of directories
files=("cage_hepg2" "cage_k562"
       "dhs_hepg2" "dhs_k562"
       "starr_hepg2" "starr_k562"
       "sure_hepg2" "sure_k562"
)


# Loop through each directory
for my_file in "${files[@]}"; do
    echo "Processing file: $my_file"
    ./Python_code/csv2npy.py -d ${dir_base} -m CR -f "${my_file}_train.csv" -j no
    ./Python_code/csv2npy.py -d ${dir_base} -m CR -f "${my_file}_val.csv" -j no
done

 (cd "${dir_base}" && wc -l *.csv > file_info.txt)

