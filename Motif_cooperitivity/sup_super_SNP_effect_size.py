import pandas as pd
import numpy as np
import seaborn as sns
import matplotlib.pyplot as plt
from scipy import stats

# read in additivity data
df=pd.read_csv("df_sup_sub.csv")
df=df[df.sig_sub_count>10].reset_index(drop=True)
df=df[df.sig_super_count>10].reset_index(drop=True)
df["track_num"]=df["track_num"].map({0: 'cage', 1: 'cage', 2: 'dhs', 3: 'dhs', 4: 'starr', 5: 'starr', 6: 'sure', 7: 'sure'})


# read in all SNPs in gnomad
tfbs_maf_enhancers_hepg2=pd.read_csv("/isdata/alab/people/pcr980/DeepCompare/Pd7_TFBS_of_maf/tfbs_maf_enhancers_hepg2.csv",header=None)
tfbs_maf_enhancers_k562=pd.read_csv("/isdata/alab/people/pcr980/DeepCompare/Pd7_TFBS_of_maf/tfbs_maf_enhancers_k562.csv",header=None)
tfbs_maf_promoters_hepg2=pd.read_csv("/isdata/alab/people/pcr980/DeepCompare/Pd7_TFBS_of_maf/tfbs_maf_promoters_hepg2.csv",header=None)
tfbs_maf_promoters_k562=pd.read_csv("/isdata/alab/people/pcr980/DeepCompare/Pd7_TFBS_of_maf/tfbs_maf_promoters_k562.csv",header=None)

tfbs_maf_enhancers_hepg2.columns=["Chromosome","Start","End","ID","REF","ALT","AF","motif","score","chip_evidence"]
tfbs_maf_enhancers_k562.columns=["Chromosome","Start","End","ID","REF","ALT","AF","motif","score","chip_evidence"]
tfbs_maf_promoters_hepg2.columns=["Chromosome","Start","End","ID","REF","ALT","AF","motif","score","chip_evidence"]
tfbs_maf_promoters_k562.columns=["Chromosome","Start","End","ID","REF","ALT","AF","motif","score","chip_evidence"]

tfbs_maf_enhancers_hepg2["dataset"]="enhancers_hepg2"
tfbs_maf_enhancers_k562["dataset"]="enhancers_k562"
tfbs_maf_promoters_hepg2["dataset"]="promoters_hepg2"
tfbs_maf_promoters_k562["dataset"]="promoters_k562"

df_tfbs_maf=pd.concat([tfbs_maf_enhancers_hepg2,tfbs_maf_enhancers_k562,tfbs_maf_promoters_hepg2,tfbs_maf_promoters_k562],axis=0).reset_index(drop=True)
df_tfbs_maf["log_AF"]=np.log10(df_tfbs_maf["AF"])




#--------------------------------------------------------------------------------------------
# global sub additive TFs and allele frequency
#--------------------------------------------------------------------------------------------

# get global sub and super TFs
series_tf_counts=df.protein.value_counts()
df_sub=df[df["super_sub_ratio"]< -0.5]
series_sub_counts=df_sub.protein.value_counts()
df_super=df[df["super_sub_ratio"]> 0.5]
series_super_counts=df_super.protein.value_counts()
df_counts=pd.concat([series_tf_counts,series_sub_counts,series_super_counts],axis=1)
df_counts.columns=["total","sub","super"]
super_tfs=df_counts[df_counts["super"]>10].index.tolist() # None
sub_tfs=df_counts[df_counts["sub"]>10].index.tolist()

# plot distribution of AF for global sub and super TFs
df_tfbs_maf["additivity"]=np.where(df_tfbs_maf["motif"].isin(sub_tfs),"sub","rest")
sns.kdeplot(data=df_tfbs_maf,x="log_AF",hue="additivity",common_norm=False)
plt.title("All putative TFBS")
plt.savefig("Plots_maf/AF_distribution_all.pdf")
plt.close()
# is there significant difference in AF between sub and rest?
u_stat, p_value = stats.mannwhitneyu(df_tfbs_maf[df_tfbs_maf["additivity"]=="sub"]["AF"],df_tfbs_maf[df_tfbs_maf["additivity"]=="rest"]["AF"]) # p<1e-183
df_tfbs_maf[df_tfbs_maf["additivity"]=="sub"]["AF"].median() # 6.97e-6
df_tfbs_maf[df_tfbs_maf["additivity"]=="rest"]["AF"].median() # 6.75e-6


# subset by ChIP_evidence=True
df_chip_true=df_tfbs_maf[df_tfbs_maf["chip_evidence"]==True].reset_index(drop=True)
sns.kdeplot(data=df_chip_true,x="log_AF",hue="additivity",common_norm=False)
plt.title("TFBS with ChIP evidence")
plt.savefig("Plots_maf/AF_distribution_chip_true.pdf")
plt.close()
# is there significant difference in AF between sub and rest?
u_stat, p_value=stats.mannwhitneyu(df_chip_true[df_chip_true["additivity"]=="sub"]["AF"],df_chip_true[df_chip_true["additivity"]=="rest"]["AF"]) # p<e-124
df_chip_true[df_chip_true["additivity"]=="sub"]["AF"].median() # 6.89e-6
df_chip_true[df_chip_true["additivity"]=="rest"]["AF"].median() # 6.76e-6


# are global sub tfs cell-type agnostic?
df_tf_dispersion_hepg2=pd.read_csv("/isdata/alab/people/pcr980/DeepCompare/TF_targets/hepG2.TFs.dispersionEstimates.tab",sep="\t")
df_tf_dispersion_k562=pd.read_csv("/isdata/alab/people/pcr980/DeepCompare/TF_targets/k562.TFs.dispersionEstimates.tab",sep="\t")
df_tf_dispersion=pd.concat([df_tf_dispersion_hepg2,df_tf_dispersion_k562],axis=0).reset_index(drop=True)
df_tf_dispersion=df_tf_dispersion.drop_duplicates().reset_index(drop=True)
# get gini rank
df_tf_dispersion["gini_rank"]=df_tf_dispersion["gini"].rank(ascending=False)
df_tf_dispersion["adjusted_dispersion_rank"]=df_tf_dispersion["adjusted_dispersion"].rank(ascending=False)
# subset by sub_tfs
df_sub=df_tf_dispersion[df_tf_dispersion["gene"].isin(sub_tfs)].reset_index(drop=True)





#---------------------------------------
# contextual sup and super TFs and allele frequency
#---------------------------------------

df["super_add"]=df["super_sub_ratio"] > 0
df["sub_add"]=df["super_sub_ratio"] < 0
df=df.groupby(["dataset","protein"]).agg({"super_add":"sum","sub_add":"sum"}).reset_index()
df.rename(columns={"super_add":"super_profile_count","sub_add":"sub_profile_count"},inplace=True)
df["tf_property"]=np.where(df["super_profile_count"]>=3,"super",np.where(df["sub_profile_count"]>=3,"sub","unknown"))
df.pivot(index="protein",columns="dataset",values="tf_property").reset_index().to_csv("tf_property.csv",index=False)
df_tfbs_maf.rename(columns={"motif":"protein"},inplace=True)
df_merged=df_tfbs_maf.merge(df,on=["dataset","protein"],how="inner")
df_merged=df_merged[df_merged["tf_property"].isin(["super","sub"])].reset_index(drop=True)
df_merged=df_merged[df_merged["chip_evidence"]==True].reset_index(drop=True)

for dataset in df_merged["dataset"].unique():
    df_subset=df_merged[df_merged["dataset"]==dataset].reset_index(drop=True).copy()
    u_stat, p_value=stats.mannwhitneyu(df_subset[df_subset["tf_property"]=="sub"]["AF"],
                                       df_subset[df_subset["tf_property"]=="super"]["AF"])
    if df_subset[df_subset["tf_property"]=="sub"]["AF"].median() > df_subset[df_subset["tf_property"]=="super"]["AF"].median():
        larger_median_group = "sub"
    else:
        larger_median_group = "super"
    df_subset["tf_property"]=pd.Categorical(df_subset["tf_property"],categories=["sub","super"])
    sns.kdeplot(data=df_subset,x="log_AF",hue="tf_property",common_norm=False)
    plt.title("Contextual TFBS for "+dataset)
    plt.annotate(f"p={p_value:.1e}",xy=(0.5,0.5),xycoords="axes fraction")
    plt.annotate(f"larger median: {larger_median_group}",xy=(0.5,0.4),xycoords="axes fraction")
    plt.savefig(f"Plots_maf/AF_distribution_contextual_{dataset}.pdf")
    plt.close()


#--------------------------------------------------------------------------------------------
# global sub additive TFs and effect size
#--------------------------------------------------------------------------------------------
maf_effect_size_enhancers_hepg2=pd.read_csv(f"maf_with_effect_size_enhancers_hepg2.csv",header=None,index_col=0)
maf_effect_size_enhancers_k562=pd.read_csv(f"maf_with_effect_size_enhancers_k562.csv",header=None,index_col=0)
maf_effect_size_promoters_hepg2=pd.read_csv(f"maf_with_effect_size_promoters_hepg2.csv",header=None,index_col=0)
maf_effect_size_promoters_k562=pd.read_csv(f"maf_with_effect_size_promoters_k562.csv",header=None,index_col=0)

maf_effect_size_enhancers_hepg2.reset_index(drop=True,inplace=True)
maf_effect_size_enhancers_k562.reset_index(drop=True,inplace=True)
maf_effect_size_promoters_hepg2.reset_index(drop=True,inplace=True)
maf_effect_size_promoters_k562.reset_index(drop=True,inplace=True)

maf_effect_size_enhancers_hepg2.columns=["chromosome","start","end","ID","REF","ALT","AF",'Name','Score','Strand']+ ["track_"+str(i) for i in range(16)]
maf_effect_size_enhancers_k562.columns=["chromosome","start","end","ID","REF","ALT","AF",'Name','Score','Strand']+ ["track_"+str(i) for i in range(16)]
maf_effect_size_promoters_hepg2.columns=["chromosome","start","end","ID","REF","ALT","AF",'Name','Score','Strand']+ ["track_"+str(i) for i in range(16)]
maf_effect_size_promoters_k562.columns=["chromosome","start","end","ID","REF","ALT","AF",'Name','Score','Strand']+ ["track_"+str(i) for i in range(16)]

