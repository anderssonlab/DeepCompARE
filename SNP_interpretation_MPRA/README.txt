Reference paper: https://www.nature.com/articles/s41467-019-11526-w
Raw_data/ is downloaded from: https://kircherlab.bihealth.org/satMutMPRA/


preprocess_MPRA.R: filter low quality positions, add reference and alternative sequences (both forward and reverse)

Pd1_MPRA_with_seq_info/: Raw MPRA data plus new column showing actual sequence. So models can predict

job_predict_SNP_effects.sh: Get DeepCompare predictions.

Pd2_DeepCompare_predictions/: DeepCompare predictions.

compare_DeepCompare_correlation.R: Compare various ways to predict effect size using DeepCompare output(forward, reverse, average). Results in Generated_plots/DeepCompare_forward_reverse_and_avg/

Pd3_Enformer_predictions/: Enformer predictions

compare_enformer_correlation.R: Compare various ways to predict effect size using Enformer output. Results in Generated_plots/Enformer_forward_reverse_and_avg/

deepCompare_vs_enformer.R: plot DeepCompare and Enformer results. Results in Generated_plots/Publication/

plot_effect_by_position.R: Plot experimental and predicted effect by position for each regulartory element.
Results in Generated_plots/Publication/






