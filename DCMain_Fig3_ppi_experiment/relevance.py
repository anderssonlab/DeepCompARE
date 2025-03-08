import pandas as pd
import numpy as np
import matplotlib.pyplot as plt

import sys
sys.path.insert(1,"/isdata/alab/people/pcr980/Scripts_python")
from utils import split_dimer
from stat_tests import fisher_exact_with_ci
from tf_cooperativity import assign_cooperativity


def match_dfs_nearest_neighbor(df_small, df_large, small_col='linear_distance', large_col='median_abs', remove_matched=True):
    # Make a copy of df_large so that we can remove rows if needed without modifying the original DataFrame.
    df_large_temp = df_large.copy()
    matched_indices = []
    # Loop through each value in the smaller DataFrame's matching column.
    for value in df_small[small_col]:
        # Calculate the absolute differences between the current value and all values in df_large_temp[large_col]
        diff = np.abs(df_large_temp[large_col] - value)
        # Identify the index of the smallest difference (nearest neighbor)
        best_idx = diff.idxmin()
        matched_indices.append(best_idx)
        # Optionally remove the matched row to avoid duplicate selections
        if remove_matched:
            df_large_temp = df_large_temp.drop(best_idx)
    return df_large.loc[matched_indices]






def fisher_exact_test(df_report,tfs_investigated,df_htfs, col):
    # reported: should be significant
    num_reported_investigated=sum(df_report[col])
    num_reported_not_investigated=sum(~df_report[col])
    num_not_reported_investigated=len(tfs_investigated)-num_reported_investigated
    num_not_reported=len(df_htfs)-len(df_report)
    num_not_reported_not_investigated=num_not_reported-num_not_reported_investigated
    table=np.array([[num_reported_investigated,num_reported_not_investigated],[num_not_reported_investigated,num_not_reported_not_investigated]])
    return fisher_exact_with_ci(table)



def relevance_of_investigation(bait, df):
    # Load all human TFs
    df_htfs=pd.read_csv(f"/isdata/alab/people/pcr980/Resource/Human_transcription_factors/DatabaseExtract_v_1.01.csv",index_col=0)
    df_htfs=df_htfs[df_htfs["Is TF?"]=="Yes"].reset_index(drop=True)
    # select K562 expressed TFs
    proteins_expressed=pd.read_csv(f"/isdata/alab/people/pcr980/DeepCompare/RNA_expression/expressed_tf_list_k562.tsv",sep='\t',header=None).iloc[:,0].values
    df_htfs=df_htfs[df_htfs["HGNC symbol"].isin(proteins_expressed)].reset_index(drop=True)
    # read df_coop
    df_coop=pd.read_csv(f"/isdata/alab/people/pcr980/DeepCompare/Pd7_TF_cooperativity/Temp/tf_pair_cooperativity_index_k562_pe.csv")
    df_coop=assign_cooperativity(df_coop,0.3,0.7)
    df_coop=df_coop[(df_coop["protein2"]==bait)].reset_index(drop=True)
    tfs_linear=df_coop[df_coop["cooperativity"]=="Linear"]["protein1"].values.tolist()
    tfs_nonlinear=df_coop[df_coop["cooperativity"]!="Linear"]["protein1"].values.tolist()
    # filter df to only include TFs
    df=df[df["gene"].isin(df_htfs["HGNC symbol"])].reset_index(drop=True)
    df["is_linear"]=df["gene"].isin(tfs_linear)
    df["is_nonlinear"]=df["gene"].isin(tfs_nonlinear)
    # enrichment test
    oddsratio_pos, pvalue_pos, ci_lower_pos, ci_higher_pos = fisher_exact_test(df, tfs_linear, df_htfs, "is_linear")
    # negative control: ChIP colocolization with similar distance
    df_chip=pd.read_csv("/isdata/alab/people/pcr980/Resource/ReMap2022/Colocolization/colocolization_k562_summary.csv")
    df_chip=df_chip[df_chip["tf2"]==bait].reset_index(drop=True)
    df_chip=df_chip[['tf1', 'tf2', 'tf_pair_count', 'median_abs']]
    # remove tf1 in tfs_investigated
    df_chip=df_chip[~df_chip["tf1"].isin(tfs_linear)].reset_index(drop=True)
    df_chip_distance_matched = match_dfs_nearest_neighbor(df_coop, df_chip)
    df["is_chip_matched"]=df["gene"].isin(df_chip_distance_matched["tf1"])
    oddsratio_neg, pvalue_neg, ci_lower_neg, ci_higher_neg = fisher_exact_test(df, df_chip_distance_matched["tf1"], df_htfs, "is_chip_matched")
    df_res=pd.DataFrame({"TF":bait,"OR_pos":[oddsratio_pos],"p_pos":[pvalue_pos],"ci_lower_pos":[ci_lower_pos],"ci_higher_pos":[ci_higher_pos],"OR_neg":[oddsratio_neg],"p_neg":[pvalue_neg],"ci_lower_neg":[ci_lower_neg],"ci_higher_neg":[ci_higher_neg]})
    return df_res




def plot_odds_ratios_with_ci(df, title, out_file, jitter=0.1):
    # Compute the log2-transformed odds ratios for both controls.
    df['log2_OR_pos'] = np.log2(df['OR_pos'])
    df['log2_OR_neg'] = np.log2(df['OR_neg'])
    # Calculate the error bars for positive control:
    # Lower error: difference between log2(OR) and log2(lower bound)
    # Upper error: difference between log2(upper bound) and log2(OR)
    df['log2_err_lower_pos'] = df['log2_OR_pos'] - np.log2(df['ci_lower_pos'])
    df['log2_err_upper_pos'] = np.log2(df['ci_higher_pos']) - df['log2_OR_pos']
    # Calculate the error bars for negative control:
    df['log2_err_lower_neg'] = df['log2_OR_neg'] - np.log2(df['ci_lower_neg'])
    df['log2_err_upper_neg'] = np.log2(df['ci_higher_neg']) - df['log2_OR_neg']
    
    # Create x positions for each TF
    x = np.arange(len(df))
    # Add small jitter for each group to avoid overlap
    x_pos = x + jitter   # jitter for positive control
    x_neg = x - jitter   # jitter for negative control
    
    plt.figure(figsize=(3,2.5))
    # Plot positive control (e.g. using black markers)
    plt.errorbar(x_pos, df['log2_OR_pos'], 
                 yerr=[df['log2_err_lower_pos'], df['log2_err_upper_pos']], 
                 fmt='o', capsize=0, markersize=2, label='Linear TFs', color='black', linewidth=0.5)
    # Plot negative control (e.g. using gray markers)
    plt.errorbar(x_neg, df['log2_OR_neg'], 
                 yerr=[df['log2_err_lower_neg'], df['log2_err_upper_neg']], 
                 fmt='o', capsize=0, markersize=2, label='Negative Control', color='gray', linewidth=0.5)
    
    # Draw a horizontal reference line at log2(OR) = 0 (which corresponds to OR = 1)
    plt.axhline(y=0, linestyle='--', color='gray', linewidth=0.5)
    
    # Set x-axis labels
    plt.xticks(x, df['TF'], fontsize=5)
    plt.yticks(fontsize=5)
    plt.xlabel('TF Bait', fontsize=7)
    plt.ylabel('log₂(Odds ratio)', fontsize=7)
    plt.title(title, fontsize=7)
    plt.legend(fontsize=5)
    plt.tight_layout()
    plt.savefig(out_file)
    plt.close()



df_res_passed=pd.DataFrame()
df_res_failed=pd.DataFrame()
for bait in ["BACH1","RFX5","IKZF1","MAFG","RREB1"]:
    # 1. get filtered proteins
    df_failed=pd.read_csv(f"Pd1_PPI_experiment/250107.{bait}.K562.Taplin.FilteredProteins.txt", sep='\t')
    df_failed["reference"]=df_failed["reference"].apply(lambda x: x.split("|")[-1])
    df_failed["gene"]=df_failed["reference"].apply(lambda x: x.split("_")[0])
    df_failed["species"]=df_failed["reference"].apply(lambda x: x.split("_")[-1])
    # remove nonhuman
    df_failed=df_failed[df_failed["species"]=="HUMAN"].reset_index(drop=True)
    # 2. get reported proteins
    df_passed=pd.read_csv(f"Pd1_PPI_experiment/250107.{bait}.K562.Taplin.GenoppiStats.txt", sep='\t')
    # 3. test relevance of investigation
    df_res_bait_failed=relevance_of_investigation(bait, df_failed)
    df_res_bait_passed=relevance_of_investigation(bait, df_passed)
    # 4. append results
    df_res_passed=pd.concat([df_res_passed,df_res_bait_passed],axis=0)
    df_res_failed=pd.concat([df_res_failed,df_res_bait_failed],axis=0)




plot_odds_ratios_with_ci(df_res_passed, "TFs passed filter", "relevance_linear_passed.pdf")
plot_odds_ratios_with_ci(df_res_failed, "TFs failed filter", "relevance_linear_failed.pdf")