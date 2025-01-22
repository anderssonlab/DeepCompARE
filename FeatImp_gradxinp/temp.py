import pandas as pd
from scipy.stats import pearsonr
from adjustText import adjust_text
import seaborn as sns
from matplotlib import pyplot as plt


def read_file(file_suffix):
    df=pd.read_csv(f'/isdata/alab/people/pcr980/DeepCompare/Pd5_motif_info/motif_info_thresh_500_{file_suffix}.csv',index_col=0)
    df["max_af"] = df["gnomad_af"].apply(lambda x: max([float(i) for i in str(x).split(":")]))
    df["max_241way"] = df["phylop_241way"].apply(lambda x: max([float(i) for i in str(x).split(":")]))
    for track_num in range(8):
        df[f"max_gradxinp_{track_num}"] = df[f"gradxinp_{track_num}"].apply(lambda x: max([float(i) for i in str(x).split(":")]))
        df[f"min_gradxinp_{track_num}"] = df[f"gradxinp_{track_num}"].apply(lambda x: min([float(i) for i in str(x).split(":")]))
        df[f"max_ism_{track_num}"] = df[f"ism_{track_num}"].apply(lambda x: max([float(i) for i in str(x).split(":")]))
        df[f"max_isa_{track_num}"] = df[f"isa_{track_num}"].apply(lambda x: max([float(i) for i in str(x).split(":")]))
    # remove columns starting with "ism" and "gradxinp" and "isa"
    cols_to_remove = [col for col in df.columns if col.startswith("ism") or col.startswith("gradxinp")]
    df.drop(cols_to_remove, axis=1, inplace=True)
    df["dataset"]=file_suffix
    return df




def scatter(df,xcol,ycol):
    sns.scatterplot(x=xcol,y=ycol,data=df)
    # annotate pearson r
    r,p=pearsonr(df[xcol],df[ycol])
    plt.annotate(f"r={r:.2f}, p={p:.2e}",(0.8,0.8),xycoords="axes fraction")
    plt.savefig(f"scatter_{xcol}_{ycol}.png")
    plt.close()


df=read_file("promoters_hepg2")
df_pred=pd.read_csv("/isdata/alab/people/pcr980/DeepCompare/Pd5_motif_info/predictions_promoters_hepg2.csv")
df_pred.columns=["region"]+[f"pred_{col}" for col in df_pred.columns[1:]]
df=pd.merge(df,df_pred,left_on="region",right_on="region",how="inner")



#---------------------
# Analysis 1: scatter plot of each importance
#---------------------


# motif importance comparison
for ycol in ["min_gradxinp_0","max_ism_0","max_isa_0","isa_track0"]:
    scatter(df,"max_gradxinp_0",ycol)


scatter(df,"isa_track0","min_gradxinp_0")
scatter(df,"max_isa_0","max_ism_0")





#---------------------
# does TF effect change according to prediction? (ELK1 example)
#---------------------
tf="ELF1"
df_sub=df[df["protein"]==tf].reset_index(drop=True)


for ycol in ["max_gradxinp_0","min_gradxinp_0","max_ism_0","max_isa_0","isa_track0"]:
    scatter(df_sub,"pred_0",ycol)




#--------------------------------------
# does TF effect change according to prediction? (all TFs)
#----------------------------------

corr_dict={"tfs":[],
           "motif_count":[],
           "max_gradxinp_0":[],
           "max_ism_0":[],
           "isa_track0":[],
           "min_gradxinp_0":[]}


for tf in df["protein"].unique():
    df_sub=df[df["protein"]==tf].reset_index(drop=True)
    corr_dict["tfs"].append(tf)
    corr_dict["motif_count"].append(len(df_sub))
    for ycol in corr_dict.keys():
        if ycol in ["tfs","motif_count"]:
            continue
        r,p=pearsonr(df_sub["pred_0"],df_sub[ycol])
        corr_dict[ycol].append(r)


df_result=pd.DataFrame(corr_dict,index=corr_dict["tfs"])


# histogram
for col in df_result.columns[2:]:
    sns.histplot(df_result[col],bins=50,fill=False)
    bottom=df_result.sort_values(col).head(int(len(df_result)*0.03))["tfs"].to_list()
    top=df_result.sort_values(col).tail(int(len(df_result)*0.03))["tfs"].to_list()
    tfs_to_write=bottom+top
    texts=[]
    for i,tf in enumerate(tfs_to_write):
        texts.append(plt.text(df_result[col][tf],0,f"{tf} (n={df_result['motif_count'][tf]})",rotation=90,color="red"))
    #
    adjust_text(texts)
    plt.title(f"correlarion between sequence activity prediction and {col}")
    plt.savefig(f"hist_corr_pred_vs_imp_{col}.png")
    plt.close()




df_result.max_gradxinp_0.median()
df_result.min_gradxinp_0.median()
df_result.max_ism_0.median()
df_result.isa_track0.median()




#-----------------------------------------------------
# Which tfs deviate most between ism and isa
#-----------------------------------------------------
# is this deviation related to prediction?

df["diff_gradxinp_ism"]=df["max_gradxinp_0"]-df["max_ism_0"]
df["diff_isa_ism"]=df["isa_track0"]-df["max_ism_0"]
# cut df_pred to 10 levels
df["pred_class"]=pd.cut(df["pred_0"],10,labels=False)
pearsonr(df["pred_0"],df["diff_gradxinp_ism"])
pearsonr(df["pred_0"],df["diff_isa_ism"])
# remove class 0 and recalculate pearson r
pearsonr(df[df["pred_class"]!=0]["pred_0"],df[df["pred_class"]!=0]["diff_gradxinp_ism"])
pearsonr(df[df["pred_class"]!=0]["pred_0"],df[df["pred_class"]!=0]["diff_isa_ism"])
# boxplot
sns.boxplot(x="pred_class",y="diff_gradxinp_ism",data=df)
plt.xticks(rotation=90)
plt.savefig("box_pred_vs_diff_gradxinp_ism.png")
plt.close()

sns.boxplot(x="pred_class",y="diff_isa_ism",data=df)
plt.xticks(rotation=90)
plt.savefig("box_pred_vs_diff_isa_ism.png")
plt.close()


# Conclusion: have to match prediction to identify other determinants of deviation between ism and isa


#-----------------------------------------------------------------------------------
# Within each class, does redundant tfs deviate more between ism and gradxinp?
#-----------------------------------------------------------------------------------

tfs_redundant=pd.read_csv("/isdata/alab/people/pcr980/DeepCompare/Pd7_TF_cooperativity/tfs_redundant_hepg2.txt",header=None)[0].to_list()
tfs_codependent=pd.read_csv("/isdata/alab/people/pcr980/DeepCompare/Pd7_TF_cooperativity/tfs_codependent_hepg2.txt",header=None)[0].to_list()
df["tf_class"]="unknown"
df.loc[df["protein"].isin(tfs_redundant),"tf_class"]="redundant"
df.loc[df["protein"].isin(tfs_codependent),"tf_class"]="codependent"


sns.boxplot(x="pred_class",y="diff_gradxinp_ism",hue="tf_class",data=df,showfliers=False)
plt.xticks(rotation=90)
plt.savefig("box_diff_gradxinp_ism_split_by_tf_class.png")
plt.close()

sns.boxplot(x="pred_class",y="diff_isa_ism",hue="tf_class",data=df,showfliers=False)
plt.xticks(rotation=90)
plt.savefig("box_diff_isa_ism_split_by_tf_class.png")
plt.close()

# calculate p values within each class, are redundant TFs more likely to deviate between ism and gradxinp?
from scipy.stats import ttest_ind
p_values=[]
for i in range(1,10):
    p_values.append(ttest_ind(df[(df["pred_class"]==i)&(df["tf_class"]=="redundant")]["diff_gradxinp_ism"],
                              df[(df["pred_class"]==i)&(df["tf_class"]=="unknown")]["diff_gradxinp_ism"]).pvalue)




#-----------------------------------------------------------------------------------
# Within each class, which tfs deviate most between ism and gradxinp?
#-----------------------------------------------------------------------------------
