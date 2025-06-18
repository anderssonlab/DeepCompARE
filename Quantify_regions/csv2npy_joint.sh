# run this in directory Quantify_regions/

dir_base="/isdata/alab/people/pcr980/DeepCompare/Datasets/"

# Define an array of directories
dirs=(
    "Dataset_5cv_with_multilabel_1"
    "Dataset_5cv_with_multilabel_2"
    "Dataset_5cv_with_multilabel_3"
    "Dataset_5cv_with_multilabel_4"
    "Dataset_5cv_with_multilabel_5"
    "Dataset_5cv_without_multilabel_1"
    "Dataset_5cv_without_multilabel_2"
    "Dataset_5cv_without_multilabel_3"
    "Dataset_5cv_without_multilabel_4"
    "Dataset_5cv_without_multilabel_5"
    "Dataset_final_rep/"
)


# Loop through each directory
for my_dir in "${dirs[@]}"; do
    echo "Processing directory: $my_dir"
    ./Python_code/csv2npy.py -d "${dir_base}$my_dir/" -m CR -f "dat_train.csv" -j yes
    ./Python_code/csv2npy.py -d "${dir_base}$my_dir/" -m CR -f "dat_val.csv" -j yes
    (cd "${dir_base}$my_dir" && wc -l *.csv > file_info.txt)
done





dirs=("Dataset_5cv_loose_1"
      "Dataset_5cv_loose_2"
      "Dataset_5cv_loose_3"
      "Dataset_5cv_loose_4"
      "Dataset_5cv_loose_5"
)

# Loop through each directory
for my_dir in "${dirs[@]}"; do
    echo "Processing directory: $my_dir"
    ./Python_code/csv2npy.py -d "${dir_base}$my_dir/" -m regression -f "dat_train.csv" -j yes
    ./Python_code/csv2npy.py -d "${dir_base}$my_dir/" -m regression -f "dat_val.csv" -j yes
    (cd "${dir_base}$my_dir" && wc -l *.csv > file_info.txt)
done
