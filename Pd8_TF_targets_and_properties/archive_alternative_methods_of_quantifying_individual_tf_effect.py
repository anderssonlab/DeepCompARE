import pandas as pd
import matplotlib.pyplot as plt
import seaborn as sns
from scipy.stats import pearsonr

import sys
sys.path.insert(1,"/isdata/alab/people/pcr980/DeepCompare/Scripts_python")
from prediction import compute_predictions


df_enhancers_hepg2=pd.read_csv(f"/isdata/alab/people/pcr980/DeepCompare/Pd5_motif_info/motif_info_thresh_500_enhancers_hepg2_track4.csv")
df_promoters_hepg2=pd.read_csv(f"/isdata/alab/people/pcr980/DeepCompare/Pd5_motif_info/motif_info_thresh_500_promoters_hepg2_track0.csv")
df_enhancers_k562=pd.read_csv(f"/isdata/alab/people/pcr980/DeepCompare/Pd5_motif_info/motif_info_thresh_500_enhancers_k562_track5.csv")
df_promoters_k562=pd.read_csv(f"/isdata/alab/people/pcr980/DeepCompare/Pd5_motif_info/motif_info_thresh_500_promoters_k562_track1.csv")

df_enhancers_hepg2["dataset"]="enhancers_hepg2"
df_promoters_hepg2["dataset"]="promoters_hepg2"
df_enhancers_k562["dataset"]="enhancers_k562"
df_promoters_k562["dataset"]="promoters_k562"

df=pd.concat([df_enhancers_hepg2,df_promoters_hepg2,df_enhancers_k562,df_promoters_k562])
# # TODO: whether to choose ChIP_evidence=True
# df=df[df["chip_evidence"]==True] 


#------------------
# Analysis 1: pred(N*600)=?
#------------------
pred_null=compute_predictions("N"*600).flatten()
# round to 2 decimal places: [-0.03, -0.08, 0.02, -0.0, 0.82, 0.83, 0.9, 0.69, -5.65, -6.07, -2.66, -1.6, -0.04, -0.84, -1.62, -7.2]
pred_null=[round(k,2) for k in pred_null] 




#------------------------------------------------------------------
# Analysis 2: Can ablating/dishuffling the context sequence quantify e_tf? No
#------------------------------------------------------------------



# def plot(df,title,out_name):
#     df_pred_ablate_context=df.loc[:,["pred_ablate_context","dataset"],].rename(columns={"pred_ablate_context":"value"})
#     df_pred_ablate_context["data"]="pred_ablate_context"

#     df_pred_dishuffle_context=df.loc[:,["pred_dishuffle_context","dataset"],].rename(columns={"pred_dishuffle_context":"value"})
#     df_pred_dishuffle_context["data"]="pred_dishuffle_context"

#     df_feat_imp=df.loc[:,["feat_imp_orig","dataset"]].rename(columns={"feat_imp_orig":"value"})
#     df_feat_imp["data"]="feat_imp"

#     df_ism_motif=df.loc[:,["ism_motif","dataset"]].rename(columns={"ism_motif":"value"})
#     df_ism_motif["data"]="ism_motif"

#     df_plot=pd.concat([
#         df_pred_ablate_context,
#         df_pred_dishuffle_context,
#         df_feat_imp,
#         df_ism_motif
#     ])
#     plt.figure(figsize=(10,8))
#     sns.boxplot(x="dataset",y="value",data=df_plot,hue="data",fliersize=0.5)
#     plt.xticks(rotation=45)
#     plt.subplots_adjust(bottom=0.3)
#     plt.title(title)
#     plt.legend(loc='upper left')
#     plt.savefig(out_name)
#     plt.close()
    
    

# for tf in df["protein"].unique():
#     df_sub=df[df["protein"]==tf]
#     plot(df_sub,f"{tf}",f"ism_distribution_chip_true_{tf}.png")
    
# plot(df,"All","ism_distribution_all_chip_true.png")

    
    
    

#------------------------------------------------------------------
# Analysis 3: Do ism add up? No
#------------------------------------------------------------------

def plot_sum_ism_vs_prediction(df, title, out_name):
    sum_ism=df.loc[:,["seq_idx", "ism_motif"]].groupby("seq_idx").sum().reset_index()    
    df_pred=df.loc[:,["seq_idx", "pred_orig"]].groupby("seq_idx").mean().reset_index()
    df_corr=pd.merge(sum_ism,df_pred,on="seq_idx")
    corr,pval=pearsonr(df_corr["ism_motif"],df_corr["pred_orig"]) # (0.52, 0.0)
    sns.scatterplot(x="ism_motif",y="pred_orig",data=df_corr)
    plt.text(0.05, 0.95, f"corr={round(corr,2)}\npval={round(pval,2)}", horizontalalignment='left', verticalalignment='top', transform=plt.gca().transAxes)
    # min_val=min(min(df_corr["ism_motif"]),min(df_corr["pred_orig"]))
    # max_val=max(max(df_corr["ism_motif"]),max(df_corr["pred_orig"]))
    # plt.plot([min_val,max_val],[min_val,max_val],color="red")
    plt.xlabel("Sum of ISM")
    plt.ylabel("Prediction")
    plt.title(title)
    plt.savefig(out_name)
    plt.close()



plot_sum_ism_vs_prediction(df,"All","ism_vs_prediction_all.png")
plot_sum_ism_vs_prediction(df[df["dataset"]=="enhancers_hepg2"],"Enhancers HepG2","ism_vs_prediction_enhancers_hepg2.png")
plot_sum_ism_vs_prediction(df[df["dataset"]=="promoters_hepg2"],"Promoters HepG2","ism_vs_prediction_promoters_hepg2.png")
plot_sum_ism_vs_prediction(df[df["dataset"]=="enhancers_k562"],"Enhancers K562","ism_vs_prediction_enhancers_k562.png")
plot_sum_ism_vs_prediction(df[df["dataset"]=="promoters_k562"],"Promoters K562","ism_vs_prediction_promoters_k562.png")


# Final conclusion: Use avg(ism)


    
# nohup python3 ism_distribution.py > ism_distribution.out &