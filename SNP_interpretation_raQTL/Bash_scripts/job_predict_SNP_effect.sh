
#!/bin/bash

# Base directory
BASE_PATH="/isdata/alab/people/pcr980/"
TASK_PATH="SNP_interpretation_raQTL/"
PYTHON_SCRIPT="${BASE_PATH}Python_code/predict_SNP_effect.py"

# Function to call the Python script with the required parameters
run_prediction() {
    local data_file=$1
    local output_file=$2
    local mode=$3

    ${PYTHON_SCRIPT} \
        -d "${BASE_PATH}${TASK_PATH}${data_file}" \
        -r "ref.seq" \
        -a "alt.seq" \
        -o "${BASE_PATH}${TASK_PATH}${output_file}" \
        -m ${mode} \
        -p "both_ends"
}

# Arrays defining the data files and corresponding output files
declare -a data_files=("Raw_data/hepg2_sure_raQTL.csv" "Raw_data/k562_sure_raQTL.csv" "Pd1_reverse_complement_of_raw_data/rc_hepg2_sure_raQTL.csv" "Pd1_reverse_complement_of_raw_data/rc_k562_sure_raQTL.csv")
declare -a modes=("classification" "regression")

# Loop through the data files and call the prediction function
for data_file in "${data_files[@]}"; do
    if [[ "$data_file" =~ hepg2 ]]; then
        cell_type="hepg2"
    else
        cell_type="k562"
    fi

    for mode in "${modes[@]}"; do
        if [[ $mode == classification ]]; then
            output_prefix="z"
        else
            output_prefix="log1p_signal"
        fi

        # Determine if the file is forward or reverse
        if [[ $data_file == Raw_data/* ]]; then
            direction="forward"
        else
            direction="reverse"
        fi

        output_file="Pd2_DeepCompare_predictions/${output_prefix}_${cell_type}_${direction}.csv"
        run_prediction "${data_file}" "${output_file}" "${mode}"
    done
done
