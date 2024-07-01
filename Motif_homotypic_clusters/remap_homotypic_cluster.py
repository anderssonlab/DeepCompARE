from scipy.stats import chi2_contingency,ks_2samp
import matplotlib.pyplot as plt
import seaborn as sns
import pandas as pd
import numpy as np
import pyranges as pr
from loguru import logger  
from statsmodels.stats.multitest import multipletests
from adjustText import adjust_text


import os
import sys
sys.path.append("/isdata/alab/people/pcr980/DeepCompare/Scripts_python")
from stat_tests import match_by_decile

#-------------
# Remap evidence
#-------------
def olap_with_homotypic_clusters(df_region,hc):
    gr_df_region = pr.PyRanges(df_region.rename(columns={'start': 'Start', 'end': 'End', 'chromosome': 'Chromosome', 'strand': 'Strand'}))
    gr_hc = pr.PyRanges(hc.rename(columns={'start': 'Start', 'end': 'End', 'chromosome': 'Chromosome', 'strand': 'element_strand'}))
    overlaps = gr_df_region.join(gr_hc)
    overlaps_df = overlaps.as_df()
    overlaps_df_unique = overlaps_df[['Start', 'End', 'Chromosome', 'Strand']].rename(columns={'Start': 'start', 'End': 'end', 'Chromosome': 'chromosome', 'Strand': 'strand'})
    result_df = pd.merge(df_region, overlaps_df_unique, on=['start', 'end', 'chromosome', 'strand'], how='left', indicator=True)
    result_df['overlap'] = result_df['_merge'] == 'both'
    result_df.drop(columns=['_merge'], inplace=True)
    return result_df['overlap']



tf_files=os.listdir("/isdata/alab/people/pcr980/Resource/JASPAR2022_tracks/Pd1_homotypic_clusters")
# tfs are the names of the files, strip the _distance_0.csv
tfs=[tf.split("_distance_0.csv")[0] for tf in tf_files]




# olaps_chip_true=[]
# olaps_chip_false=[]
# pvals=[]
# tf_list=[]

# for tf in tfs:
#     logger.info(f"Processing {tf}")
#     df_region=pd.read_csv(f"/isdata/alab/people/pcr980/DeepCompare/Motif_gc_context/Pd1_motif_loc_and_context_gc/{tf}.csv",
#                           skiprows=1,
#                           header=None,
#                           index_col=False,
#                           names=['start','end','protein','score','strand','chromosome','chip_evidence','motif_seq','motif_gc'
#                                  'context_gc_2bp','context_gc_10bp','context_gc_50bp','context_gc_100bp','context_gc_300bp'])
#     logger.info(f"Before filter: {df_region.shape}")
#     # if chip_evidence contain only False, no True, skip
#     if df_region['chip_evidence'].sum()==0:
#         continue
#     hc=pd.read_csv(f"/isdata/alab/people/pcr980/Resource/JASPAR2022_tracks/Pd1_homotypic_clusters/{tf}_distance_0.csv",header=None)
#     hc.columns=['chromosome', 'start', 'end', 'element_start', 'element_end', 'element_strand', 'element_score', 'count', 'protein']
#     # since hc is partially done, we need to filter out the clusters that are not in the same chromosome
#     logger.info(hc['chromosome'].unique())
#     df_region=df_region[df_region['chromosome'].isin(hc['chromosome'].unique())]
#     df_region=match_by_decile(df_region,'score','chip_evidence')
#     logger.info(f"After filter: {df_region.shape}")
#     df_region['overlap']=olap_with_homotypic_clusters(df_region,hc)
#     # percentage of overlap=True and chip_evidence=True  in chip_evidence=True 
#     olaps_chip_true.append(df_region[df_region['chip_evidence']==True]['overlap'].sum()/df_region['chip_evidence'].sum())
#     olaps_chip_false.append(df_region[df_region['chip_evidence']==False]['overlap'].sum()/(df_region['chip_evidence']==False).sum())
#     # chi-square test
#     obs = pd.crosstab(df_region['chip_evidence'], df_region['overlap'])
#     chi2, p, dof, ex = chi2_contingency(obs)
#     pvals.append(p)
#     tf_list.append(tf)

# df_res=pd.DataFrame({'tf':tf_list,'olaps_chip_true':olaps_chip_true,'olaps_chip_false':olaps_chip_false,'pvals':pvals})
# df_res.to_csv("remap_evidence.csv",index=False)



#---------------------
# DeepCompare evidence
#---------------------

df_promoters_k562=pd.read_csv("/isdata/alab/people/pcr980/DeepCompare/Pd5_motif_info/motif_info_promoters_k562.csv")
df_promoters_k562["file_name"]="promoters_k562"

df_promoters_hepg2=pd.read_csv("/isdata/alab/people/pcr980/DeepCompare/Pd5_motif_info/motif_info_promoters_hepg2.csv")
df_promoters_hepg2["file_name"]="promoters_hepg2"

df_enhancers_k562=pd.read_csv("/isdata/alab/people/pcr980/DeepCompare/Pd5_motif_info/motif_info_enhancers_k562.csv")
df_enhancers_k562["file_name"]="enhancers_k562"

df=pd.concat([df_promoters_k562,df_promoters_hepg2,df_enhancers_k562],axis=0)

true_medians=[]
false_medians=[]
dstats=[]
pvals=[]
tf_list=[]


for tf in df.protein.unique():
    logger.info(f"Processing {tf}")
    motif_df_sub=df[df['protein']==tf]
    motif_df_sub=match_by_decile(motif_df_sub,'score','chip_evidence')
    if motif_df_sub.shape[0]==0:
        continue
    # if chip_evidence contain only False, or only True, skip
    if len(motif_df_sub['chip_evidence'].unique())==1:
        continue
    #overlap_true=motif_df_sub[motif_df_sub['homotypic_environment'].notna()]['feat_imp_orig']
    #overlap_false=motif_df_sub[pd.isna(motif_df_sub['homotypic_environment'])]['feat_imp_orig']
    overlap_true=motif_df_sub[motif_df_sub['homotypic_clusters']==True]['feat_imp_orig']
    overlap_false=motif_df_sub[motif_df_sub['homotypic_clusters']==False]['feat_imp_orig']
    
    logger.info(f"Overlap True: {len(overlap_true)}")
    logger.info(f"Overlap False: {len(overlap_false)}")
    if len(overlap_true)==0 or len(overlap_false)==0:
        continue
    dstat,pval=ks_2samp(overlap_true,overlap_false)
    true_medians.append(overlap_true.median())
    false_medians.append(overlap_false.median())
    dstats.append(dstat)
    pvals.append(pval)
    tf_list.append(tf)

df_res=pd.DataFrame({'tf':tf_list,'true_medians':true_medians,'false_medians':false_medians,'dstats':dstats,'pvals':pvals})
df_res["fdr"]=multipletests(df_res.pvals, method="fdr_bh")[1]

df_res.to_csv("deepCompare_evidence.csv",index=False)

# plot df_res
# scatter plot. x=true_medians, y=false_medians
# color by fdr: if fdr<0.05, color=red, else color=blue, add legend
# add text to the points with tf names
df_res=pd.read_csv("deepCompare_evidence.csv")
df_res["sig"]=np.where(df_res["fdr"]<0.05,"Sig","Nonsig")
df_res["log_true_medians"]=np.log(df_res["true_medians"])
df_res["log_false_medians"]=np.log(df_res["false_medians"])


texts = []
categories = df_res['sig'].unique()

# Generate a color palette with seaborn
palette = sns.color_palette("husl", len(categories))
palette_dict = dict(zip(categories, palette))

for category in categories:
    # Filter data for the current category
    category_data = df_res[df_res['sig'] == category]
    # Plot each category with its corresponding color
    scatter = plt.scatter(category_data['log_true_medians'], category_data['log_false_medians'],
                          color=palette_dict[category], label=category,s=1)
    # Add text for each point in the current category
    for line in range(0, category_data.shape[0]):
        text = plt.text(category_data['log_true_medians'].iloc[line], category_data['log_false_medians'].iloc[line], 
                        category_data['tf'].iloc[line], horizontalalignment='left', 
                        size='x-small', color=palette_dict[category])
        texts.append(text)

# Adjust texts to avoid overlap
adjust_text(texts)
plt.legend()
# add diagonal line
plt.plot([-8,-0.5],[-8,-0.5],color="black")
plt.xlabel("log_true_medians")
plt.ylabel("log_false_medians")
plt.title("Medians of feat_imp_orig grouped by homotypic clusters")
plt.savefig("DeepCompare_evidence.pdf",dpi=300)
plt.close()