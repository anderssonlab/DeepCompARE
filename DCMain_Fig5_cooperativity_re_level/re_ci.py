import pandas as pd
import numpy as np
from loguru import logger
import seaborn as sns
import matplotlib.pyplot as plt
from scipy.stats import mannwhitneyu
import sys
sys.path.insert(1,"/isdata/alab/people/pcr980/Scripts_python")
from tf_cooperativity import read_cooperativity


# ----------------------------------------------------
# Helper functions
# ----------------------------------------------------

def get_re_ci(file_name):
    # get cooperative pairs
    df_coop=read_cooperativity(f"Pd2_mutate_pairs/mutate_pairs_{file_name}.csv",threshold=0.01, track_nums=[1,3,5,7])
    # group by region_idx and codependency, sum column c within group
    df_coop=df_coop.groupby(["region_idx","codependency"]).agg({"c":"sum"}).unstack(fill_value=0)
    # reduce to single index
    df_coop.columns = df_coop.columns.droplevel(0)
    df_coop.reset_index(inplace=True)
    df_coop.columns=["region_idx","redundant","codependent"]
    df_coop["redundant"]=df_coop["redundant"].abs()
    df_coop=df_coop[(df_coop["redundant"]+df_coop["codependent"])>0].reset_index(drop=True)
    df_coop["ci_re"]=df_coop["codependent"]/(df_coop["codependent"]+df_coop["redundant"])
    # add region information
    return df_coop["ci_re"].tolist()






# ----------------------------------------------------
# E1E2P1P2
# ----------------------------------------------------
ci_e1=get_re_ci("E1_k562")
logger.info(f"ci_e1: {len(ci_e1)}")
ci_e2=get_re_ci("E2_k562")
logger.info(f"ci_e2: {len(ci_e2)}")
ci_p1=get_re_ci("P1_k562")
logger.info(f"ci_p1: {len(ci_p1)}")
ci_p2=get_re_ci("P2_k562")
logger.info(f"ci_p2: {len(ci_p2)}")


df_res=pd.DataFrame({"ci":ci_e1+ci_e2+ci_p1+ci_p2,
                     "region":["E1"]*len(ci_e1)+["E2"]*len(ci_e2)+["P1"]*len(ci_p1)+["P2"]*len(ci_p2)})

df_res["region"]=pd.Categorical(df_res["region"],categories=["E1","E2","P1","P2"],ordered=True)

df_res.to_csv("ci_E1E2P1P2.csv",index=False)


# violin plot
plt.figure(figsize=(6,8))
sns.boxplot(x="region", y="ci", data=df_res)
plt.title("Cooperativity index of regulatory elements")
plt.ylabel("Cooperativity index")
plt.xlabel("Region")
plt.savefig("ci_E1E2P1P2.pdf")
plt.close()





# ----------------------------------------------------
# constrained v.s. non-constrained
# ----------------------------------------------------
constrained_distal_ti=get_re_ci("dhs_constrained_distal_ti_k562")
logger.info(f"ci_proximal_constrained: {len(constrained_distal_ti)}")   
constrained_distal_ts=get_re_ci("dhs_constrained_distal_ts_k562")
logger.info(f"ci_proximal_non_constrained: {len(constrained_distal_ts)}")
nonconstrained_distal_ti=get_re_ci("dhs_nonconstrained_distal_ti_k562")
logger.info(f"ci_distal_constrained: {len(nonconstrained_distal_ti)}")
nonconstrained_distal_ts=get_re_ci("dhs_nonconstrained_distal_ts_k562")
logger.info(f"ci_distal_non_constrained: {len(nonconstrained_distal_ts)}")
constrained_proximal_ti=get_re_ci("dhs_constrained_proximal_ti_k562")
logger.info(f"ci_proximal_constrained: {len(constrained_proximal_ti)}")
constrained_proximal_ts=get_re_ci("dhs_constrained_proximal_ts_k562")
logger.info(f"ci_proximal_non_constrained: {len(constrained_proximal_ts)}")
nonconstrained_proximal_ti=get_re_ci("dhs_nonconstrained_proximal_ti_k562")
logger.info(f"ci_proximal_constrained: {len(nonconstrained_proximal_ti)}")
nonconstrained_proximal_ts=get_re_ci("dhs_nonconstrained_proximal_ts_k562")
logger.info(f"ci_proximal_non_constrained: {len(nonconstrained_proximal_ts)}")


df_res=pd.DataFrame({"ci":constrained_distal_ti+constrained_distal_ts+nonconstrained_distal_ti+nonconstrained_distal_ts+constrained_proximal_ti+constrained_proximal_ts+nonconstrained_proximal_ti+nonconstrained_proximal_ts,
                    "region":["constrained_distal_ti"]*len(constrained_distal_ti)+["constrained_distal_ts"]*len(constrained_distal_ts)+["nonconstrained_distal_ti"]*len(nonconstrained_distal_ti)+["nonconstrained_distal_ts"]*len(nonconstrained_distal_ts)+["constrained_proximal_ti"]*len(constrained_proximal_ti)+["constrained_proximal_ts"]*len(constrained_proximal_ts)+["nonconstrained_proximal_ti"]*len(nonconstrained_proximal_ti)+["nonconstrained_proximal_ts"]*len(nonconstrained_proximal_ts)})
df_res.to_csv("re_ci.csv",index=False)



df_res["constraint"]=df_res["region"].apply(lambda x: x.split("_")[0])
df_res["region"]=df_res["region"].apply(lambda x: x.split("_")[1]+"_"+x.split("_")[2])



plt.figure(figsize=(6,8))
sns.boxplot(x="region", y="ci", data=df_res, hue="constraint")
plt.title("CI of regulatory elements")
plt.ylabel("Cooperativity index")
plt.xlabel("Region")
plt.xticks(rotation=90)
plt.tight_layout()
plt.savefig("re_ci.pdf")
plt.close()




# ----------------------------------------------------
# ti v.s. ts
# ----------------------------------------------------

ci_proximal_ti=get_re_ci("proximal_ti_k562")
logger.info(f"ci_proximal_ti: {len(ci_proximal_ti)}")
ci_proximal_ts=get_re_ci("proximal_ts_k562")
logger.info(f"ci_proximal_ts: {len(ci_proximal_ts)}")
ci_distal_ti=get_re_ci("distal_ti_k562")
logger.info(f"ci_distal_ti: {len(ci_distal_ti)}")
ci_distal_ts=get_re_ci("distal_ts_k562")
logger.info(f"ci_distal_ts: {len(ci_distal_ts)}")

mannwhitneyu(ci_proximal_ti,ci_proximal_ts)
np.median(ci_proximal_ti)
np.median(ci_proximal_ts)

mannwhitneyu(ci_distal_ti,ci_distal_ts)
np.median(ci_distal_ti)
np.median(ci_distal_ts)

df_res=pd.DataFrame({"ci":ci_proximal_ti+ci_proximal_ts+ci_distal_ti+ci_distal_ts,
                    "region":["proximal_ti"]*len(ci_proximal_ti)+["proximal_ts"]*len(ci_proximal_ts)+["distal_ti"]*len(ci_distal_ti)+["distal_ts"]*len(ci_distal_ts)})

df_res.to_csv("ci_ti_ts.csv",index=False)


plt.figure(figsize=(6,8))
sns.boxplot(x="region", y="ci", data=df_res)
plt.title("CI of regulatory elements")
plt.ylabel("Cooperativity index")
plt.xlabel("Region")
plt.xticks(rotation=90)
plt.tight_layout()
plt.savefig("ci_ti_ts.pdf")
plt.close()




# nohup python3 re_ci.py > re_ci.out &












# ----------------------------------------------------
# 3. Why do promoter and enhancer differ in additivity profile?
# H1: same TF pair has different additivity profile 
# H2: different TF pair distribution
# ----------------------------------------------------

# cell_type="hepg2"
# remap="Hep-G2"

# cell_type="k562"
# remap="K-562"

# df_promoter=aggregate_tracks("promoters",cell_type,remap)
# # group by protein1, protein2, sum sub_additive and super_additive
# df_promoter=df_promoter.groupby(["protein1","protein2"])[["sub_additive","super_additive"]].sum().reset_index()
# df_promoter.rename(columns={"sub_additive":"sub_additive_promoter","super_additive":"super_additive_promoter"},inplace=True)
# df_enhancer=aggregate_tracks("enhancers",cell_type,remap)
# df_enhancer=df_enhancer.groupby(["protein1","protein2"])[["sub_additive","super_additive"]].sum().reset_index()
# df_enhancer.rename(columns={"sub_additive":"sub_additive_enhancer","super_additive":"super_additive_enhancer"},inplace=True)

# # merge promoter and enhancer by 'region_idx', 'protein1', 'protein2', 'chromosome1', 'start_rel1','end_rel1', 'strand', 'score', 'chromosome2', 'start_rel2', 'end_rel2', 'score2'
# df=pd.merge(df_promoter,df_enhancer,on=['protein1','protein2'],how="outer")
# # replace NaN with 0
# df.fillna(0,inplace=True)


# df["additivity_promoter"]="unknown"
# df.loc[df["sub_additive_promoter"]>df["super_additive_promoter"],"additivity_promoter"]="sub"
# df.loc[df["sub_additive_promoter"]<df["super_additive_promoter"],"additivity_promoter"]="super"

# df["additivity_enhancer"]="unknown"
# df.loc[df["sub_additive_enhancer"]>df["super_additive_enhancer"],"additivity_enhancer"]="sub"
# df.loc[df["sub_additive_enhancer"]<df["super_additive_enhancer"],"additivity_enhancer"]="super"

# df_known=df[(df["additivity_promoter"]!="unknown") & (df["additivity_enhancer"]!="unknown")].reset_index(drop=True)
# # how many TF pairs have different additivity profile
# df_diff=df_known[df_known["additivity_promoter"]!=df_known["additivity_enhancer"]].reset_index(drop=True)



# Conclusions: 
# For HepG2 
# union(tf_pair_promoter,tf_pair_enhancer)=7068
# intersection(tf_pair_promoter,tf_pair_enhancer)=2264
# inconsistent TF_pair behavior between promoter and enhancer: 802

# For K562:
# union(tf_pair_promoter,tf_pair_enhancer)=10804
# intersection(tf_pair_promoter,tf_pair_enhancer)=4050
# inconsistent TF_pair behavior between promoter and enhancer: 1320

