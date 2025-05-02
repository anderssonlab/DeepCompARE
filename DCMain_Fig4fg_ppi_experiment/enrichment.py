import pandas as pd
import numpy as np
import seaborn as sns
import matplotlib.pyplot as plt

from utils_ppi import read_pooled_found_tf_pairs, mannwhitneyu_with_nan, read_pooled_subthresh_tf_pairs

import sys
sys.path.insert(1,"/isdata/alab/people/pcr980/Scripts_python")
from tf_cooperativity import assign_cooperativity
from stat_tests import fisher_exact_with_ci


import matplotlib
matplotlib.rcParams['pdf.fonttype']=42

# ci_suffix="dhs"
# redundancy_threshold=0.44
# codependent_threshold=0.81



ci_suffix="pe"
redundancy_threshold=0.3 
codependent_threshold=0.7




mode="cooperativity_index"






df_coop=pd.read_csv(f"/isdata/alab/people/pcr980/DeepCompare/Pd7_TF_cooperativity/tf_pair_cooperativity_index_k562_{ci_suffix}.csv")
df_coop=assign_cooperativity(df_coop,1,0.9,redundancy_threshold,codependent_threshold)
df_coop=df_coop[df_coop["protein2"].isin(["BACH1","MAFG","IKZF1","RREB1","RFX5"])].reset_index(drop=True)

# subset for proper background TFs (proteomics detectable), 
# because TFs detected in proteomics tend to have lower ci for the 5 baits
proteins_background=pd.read_csv('Pd1_PPI_experiment/Nusinow_CCLE_K562proteomics.txt', sep='\t')["Gene_Symbol"].tolist()
df_coop=df_coop[df_coop['protein1'].isin(proteins_background)].reset_index(drop=True)
# remove cases where protein1=protein2
df_coop=df_coop[df_coop["protein1"]!=df_coop["protein2"]].reset_index(drop=True)
df_coop["prediction_state"]="predicted"




# read pooled ppi pairs
df_found_pooled=read_pooled_found_tf_pairs()
df_subthreshold_pooled=read_pooled_subthresh_tf_pairs()
df_experiment=pd.concat([df_found_pooled,df_subthreshold_pooled],axis=0)


# merge with cooperativity index
df_coop=df_coop.merge(df_experiment,on=["protein1","protein2"],how="outer")
# if prediction_state is nan, set to "not_predicted"
# if experiment_state is nan, set to "not_found"
df_coop["prediction_state"]=df_coop["prediction_state"].fillna("not_predicted")
df_coop["experiment_state"]=df_coop["experiment_state"].fillna("not_found")





# plot distribution of cooperativity index
sns.kdeplot(df_coop[mode],color="blue",cut=0,label="All")
sns.kdeplot(df_coop[df_coop["experiment_state"]=="significant"][mode],color="red",cut=0,label="Significant TFs")
sns.kdeplot(df_coop[df_coop["experiment_state"]!="not_found"][mode],color="green",cut=0,label="Significant + Subthreshold TFs")
# add mann-whitney u test
stat,p=mannwhitneyu_with_nan(df_coop[mode],df_coop[df_coop["experiment_state"]=="significant"][mode])
plt.text(0.5,0.5,f"All v.s. sig p={p:.2e}",transform=plt.gca().transAxes)
stat,p=mannwhitneyu_with_nan(df_coop[mode],df_coop[df_coop["experiment_state"]=="subthreshold"][mode])
plt.text(0.5,0.4,f"All v.s. lenient p={p:.2e}",transform=plt.gca().transAxes)
plt.xlabel(mode)
plt.ylabel("Density")
plt.title(f"Distribution of {mode} ({ci_suffix})")
plt.legend()
plt.savefig(f"pooled_{mode}_distribution_{ci_suffix}.pdf")
plt.close()




# enrichment analysis: are codependent pairs more likely to be found?
def fisher_exact_test(df, class_label):
    # Count occurrences for the given class
    found_class = df[(df["cooperativity"] == class_label) & (df["experiment_state"] == "found")].shape[0]
    not_found_class = df[(df["cooperativity"] == class_label) & (df["experiment_state"] == "not_found")].shape[0]
    # Count occurrences for other classes
    found_other = df[(df["cooperativity"] != class_label) & (df["experiment_state"] ==  "found")].shape[0]
    not_found_other = df[(df["cooperativity"] != class_label) & (df["experiment_state"] == "not_found")].shape[0]
    # Construct contingency table
    contingency_table = np.array([[found_class, not_found_class], [found_other, not_found_other]])
    # Perform Fisher's exact test
    odds_ratio, p, ci_lower, ci_upper = fisher_exact_with_ci(contingency_table)
    df_res=pd.DataFrame({"class":[class_label],"odds_ratio":[odds_ratio],"p":[p],"ci_lower":[ci_lower],"ci_upper":[ci_upper]})
    return df_res




# are codependent TFs more likely found than other predicted class?
# remove the class of "not_predicted"
df_coop=df_coop[df_coop["prediction_state"]!="not_predicted"].reset_index(drop=True)

df_res=pd.DataFrame()
for class_label in ["Independent","Redundant","Intermediate","Synergistic"]:
    df_coop_copy=df_coop.copy()
    # stringent version: "experiment_state"=="found" only if "experiment_state"=="significant"
    df_coop_copy.loc[df_coop_copy["experiment_state"]=="significant","experiment_state"]="found"
    df_res_class=fisher_exact_test(df_coop_copy, class_label)
    df_res_class["method"]="stringent"
    df_res=pd.concat([df_res,df_res_class],axis=0).reset_index(drop=True)
    # lenient version: "experiment_state"=="found" if "experiment_state"=="significant" or "experiment_state"=="subthreshold"
    df_coop_copy.loc[df_coop_copy["experiment_state"]=="subthreshold","experiment_state"]="found"
    df_res_class=fisher_exact_test(df_coop_copy, class_label)
    df_res_class["method"]="lenient"
    df_res=pd.concat([df_res,df_res_class],axis=0).reset_index(drop=True)
    










plt.figure(figsize=(2.3,2.6))
# thin frame
plt.gca().spines['top'].set_linewidth(0.5)
plt.gca().spines['right'].set_linewidth(0.5)
plt.gca().spines['bottom'].set_linewidth(0.5)
plt.gca().spines['left'].set_linewidth(0.5)
# plot df_res

jitter = 0.1
x = np.arange(len(df_res["class"].unique()))
# Add small jitter for each group to avoid overlap
x_lenient = x + jitter   # jitter for positive control
x_stringent = x - jitter   # jitter for negative control

df_stringent=df_res[df_res["method"]=="stringent"].reset_index(drop=True)
df_lenient=df_res[df_res["method"]=="lenient"].reset_index(drop=True)
plt.errorbar(x_stringent, df_stringent['odds_ratio'],yerr=[df_stringent['odds_ratio']-df_stringent['ci_lower'], df_stringent['ci_upper']-df_stringent['odds_ratio']],fmt='o', capsize=0, markersize=2, label='Significant TFs', color='lightgray', linewidth=0.5)
plt.errorbar(x_lenient, df_lenient['odds_ratio'],yerr=[df_lenient['odds_ratio']-df_lenient['ci_lower'], df_lenient['ci_upper']-df_lenient['odds_ratio']],fmt='o', capsize=0, markersize=2, label='Significant + Subthreshold TFs', color='black', linewidth=0.5)
plt.xticks(x, df_stringent['class'], fontsize=5)
plt.yticks(fontsize=5)
plt.axhline(y=1, linestyle='--', color='gray', linewidth=0.5)
plt.xlabel('Type of TF partner', fontsize=7)
plt.ylabel('Odds ratio', fontsize=7)
plt.legend(fontsize=5)
plt.title("Pooled interactors of all 5 TFs", fontsize=7)
plt.tight_layout()
plt.savefig(f"pooled_enrichment_{ci_suffix}.pdf")
plt.close()


