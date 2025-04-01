
BASE_PATH="/isdata/alab/people/pcr980/"
TASK_PATH="SNP_interpretation_MPRA/"

${BASE_PATH}Python_code/predict_SNP_effect.py \
	-d ${BASE_PATH}${TASK_PATH}Pd1_MPRA_with_seq_info/MPRA_with_forward_seq.csv\
    -r "ref_seq" \
	-a "alt_seq" \
	-o ${BASE_PATH}${TASK_PATH}Pd2_DeepCompare_predictions/z_forward.csv\
	-m classification \
	-p "both_ends"


${BASE_PATH}Python_code/predict_SNP_effect.py \
	-d ${BASE_PATH}${TASK_PATH}Pd1_MPRA_with_seq_info/MPRA_with_reverse_seq.csv\
    -r "ref_seq" \
	-a "alt_seq" \
	-o ${BASE_PATH}${TASK_PATH}Pd2_DeepCompare_predictions/z_reverse.csv\
	-m classification \
	-p "both_ends"


${BASE_PATH}Python_code/predict_SNP_effect.py \
	-d ${BASE_PATH}${TASK_PATH}Pd1_MPRA_with_seq_info/MPRA_with_forward_seq.csv\
    -r "ref_seq" \
	-a "alt_seq" \
	-o ${BASE_PATH}${TASK_PATH}Pd2_DeepCompare_predictions/log1p_signal_forward.csv\
	-m regression \
	-p "both_ends"


${BASE_PATH}Python_code/predict_SNP_effect.py \
	-d ${BASE_PATH}${TASK_PATH}Pd1_MPRA_with_seq_info/MPRA_with_reverse_seq.csv\
    -r "ref_seq" \
	-a "alt_seq" \
	-o ${BASE_PATH}${TASK_PATH}Pd2_DeepCompare_predictions/log1p_signal_reverse.csv\
	-m regression \
	-p "both_ends"

