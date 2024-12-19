import pandas as pd
import numpy as np
from loguru import logger
import seaborn as sns
import matplotlib.pyplot as plt
from scipy.stats import mannwhitneyu



import sys
sys.path.insert(1,"/isdata/alab/people/pcr980/Scripts_python")
from stat_tests import bin_and_label

# ----------------------------------------------------
# Helper functions
# ----------------------------------------------------



def get_tf_count(file_name):
    df=pd.read_csv(f"Pd2_motif_info/motif_info_thresh_500_{file_name}.tsv.csv")
    if "hepg2" in file_name:
        cell_line="hepg2"
    elif "k562" in file_name:
        cell_line="k562"
    else:
        raise ValueError("cell line not found")
    tfs_codependent=pd.read_csv(f"/isdata/alab/people/pcr980/DeepCompare/Pd7_TF_cooperativity/tfs_codependent_{cell_line}.txt",header=None).iloc[:,0].tolist()
    tfs_redundant=pd.read_csv(f"/isdata/alab/people/pcr980/DeepCompare/Pd7_TF_cooperativity/tfs_redundant_{cell_line}.txt",header=None).iloc[:,0].tolist()
    df["tf_type"]="other"
    df.loc[df["protein"].isin(tfs_codependent),"tf_type"]="codependent"
    df.loc[df["protein"].isin(tfs_redundant),"tf_type"]="redundant"
    # count number of redundant and codependent TFs in each region
    df=df.groupby(["region","tf_type"]).size().unstack(fill_value=0).reset_index()
    df.rename(columns={"codependent":"codependent_tf_count","redundant":"redundant_tf_count","other":"other_tf_count"},inplace=True)
    # remove rows with nan
    df["region_type"]=file_name
    return df



# ----------------------------------------------------
# get data
# ----------------------------------------------------
for cell_line in ["hepg2","k562"]:
    constrained_distal_ti=get_tf_count(f"dhs_constrained_distal_ti_{cell_line}")
    constrained_distal_ts=get_tf_count(f"dhs_constrained_distal_ts_{cell_line}")
    nonconstrained_distal_ti=get_tf_count(f"dhs_nonconstrained_distal_ti_{cell_line}")
    nonconstrained_distal_ts=get_tf_count(f"dhs_nonconstrained_distal_ts_{cell_line}")
    constrained_proximal_ti=get_tf_count(f"dhs_constrained_proximal_ti_{cell_line}")
    constrained_proximal_ts=get_tf_count(f"dhs_constrained_proximal_ts_{cell_line}")
    nonconstrained_proximal_ti=get_tf_count(f"dhs_nonconstrained_proximal_ti_{cell_line}")
    nonconstrained_proximal_ts=get_tf_count(f"dhs_nonconstrained_proximal_ts_{cell_line}")
    df=pd.concat([constrained_distal_ti,constrained_distal_ts,nonconstrained_distal_ti,nonconstrained_distal_ts,constrained_proximal_ti,constrained_proximal_ts,nonconstrained_proximal_ti,nonconstrained_proximal_ts],ignore_index=True)
    df.to_csv(f"tf_count_{cell_line}.csv",index=False)



# ----------------------------------------------------
# boxplot for distribution of TF count
# ----------------------------------------------------


# Neighbor pairs for Mann-Whitney U tests
neighbor_pairs = [
    ("distal_ts", "distal_ti"),
    ("distal_ti", "proximal_ts"),
    ("proximal_ts", "proximal_ti")
]

for cell_line in ["k562", "hepg2"]:
    df = pd.read_csv(f"tf_count_{cell_line}.csv")
    df["region_type"] = df["region_type"].apply(lambda x: x.split("_")[2] + "_" + x.split("_")[3])
    df["region_type"] = pd.Categorical(
        df["region_type"],
        categories=["distal_ts", "distal_ti", "proximal_ts", "proximal_ti"],
        ordered=True
    )
    # Add pseudocount
    df["redundant_tf_count"] += 1
    df["codependent_tf_count"] += 1
    df["other_tf_count"] += 1
    #
    for col in ["redundant_tf_count", "codependent_tf_count", "other_tf_count"]:
        logger.info(f"Plotting {col}")
        plt.figure(figsize=(6, 6))
        sns.boxplot(x="region_type", y=col, data=df, boxprops={'facecolor': 'none'})
        plt.xticks(rotation=90)
        plt.yscale("log")
        #
        # Calculate and annotate Mann-Whitney U test p-values
        for pair in neighbor_pairs:
            group1 = df[df["region_type"] == pair[0]][col]
            group2 = df[df["region_type"] == pair[1]][col]
            stat, p_value = mannwhitneyu(group1, group2, alternative='two-sided')
            # Add p-value annotation
            y_max = max(group1.max(), group2.max())
            y_position = y_max * 1.5  # Adjust the multiplier as needed for spacing
            #
            plt.plot([
                df["region_type"].cat.categories.get_loc(pair[0]),
                df["region_type"].cat.categories.get_loc(pair[1])
            ], [y_position, y_position], lw=1.5, color='black')
            #
            plt.text(
                (df["region_type"].cat.categories.get_loc(pair[0]) + df["region_type"].cat.categories.get_loc(pair[1])) / 2,
                y_position * 1.1,
                f"p={p_value:.2e}",
                ha='center',
                va='bottom',
                fontsize=8
            )
            #
        plt.tight_layout()
        plt.subplots_adjust(bottom=0.3)  # Larger lower margin
        plt.savefig(f"Plots/{col}_{cell_line}.pdf")
        plt.close()


# ----------------------------------------------------
# statistical test
# ----------------------------------------------------

cell_line="k562"
df=pd.read_csv(f"tf_count_{cell_line}.csv")
df["region_type"]=df["region_type"].replace(f"_{cell_line}","",regex=True)
df["distal_proximal"]=df["region_type"].apply(lambda x: x.split("_")[2])
df["ti_ts"]=df["region_type"].apply(lambda x: x.split("_")[3])


# mannwhitneyu test on TF count for distal and proximal
df_class1=df[df["distal_proximal"]=="distal"]
df_class2=df[df["distal_proximal"]=="proximal"]

mannwhitneyu(df_class1["redundant_tf_count"],df_class2["redundant_tf_count"],alternative="two-sided")
df_class1["redundant_tf_count"].median()
df_class2["redundant_tf_count"].median()

mannwhitneyu(df_class1["codependent_tf_count"],df_class2["codependent_tf_count"],alternative="two-sided")
df_class1["codependent_tf_count"].median()
df_class2["codependent_tf_count"].median()

mannwhitneyu(df_class1["other_tf_count"],df_class2["other_tf_count"],alternative="two-sided")
df_class1["other_tf_count"].median()
df_class2["other_tf_count"].median()



# mannwhitneyu test on TF count for ti and ts
df_class1=df[df["ti_ts"]=="ti"]
df_class2=df[df["ti_ts"]=="ts"]

mannwhitneyu(df_class1["redundant_tf_count"],df_class2["redundant_tf_count"],alternative="two-sided")
df_class1["redundant_tf_count"].median()
df_class2["redundant_tf_count"].median()

mannwhitneyu(df_class1["codependent_tf_count"],df_class2["codependent_tf_count"],alternative="two-sided")
df_class1["codependent_tf_count"].median()
df_class2["codependent_tf_count"].median()











# ----------------------------------------------------
# classifier
# ----------------------------------------------------

cell_line="k562"
df=pd.read_csv(f"tf_count_{cell_line}.csv")
df["region_type"]=df["region_type"].replace(f"_{cell_line}","",regex=True)
df["distal_proximal"]=df["region_type"].apply(lambda x: x.split("_")[2])
df["ti_ts"]=df["region_type"].apply(lambda x: x.split("_")[3])
df["region_type"]=df["region_type"].apply(lambda x: x.split("_")[2]+"_"+x.split("_")[3])


pred_col = "region_type"
x = df[['codependent_tf_count', 'redundant_tf_count']]
y = df[pred_col]

# Encode target labels
le = LabelEncoder()
y = le.fit_transform(y)

# Get class counts for weights
class_counts = df[pred_col].value_counts()
class_weights = {le.transform([cls])[0]: 1000.0 / count for cls, count in class_counts.items()}

# Initialize classifier with class weights
clf = RandomForestClassifier(random_state=0)
scorers = {
    'f1': make_scorer(f1_score, average='weighted'),  # F1 score
    'accuracy': make_scorer(accuracy_score)          # Accuracy
}

cv_results = cross_validate(clf, x, y, cv=5, scoring=scorers)
f1_scores = cv_results['test_f1']
accuracy_scores = cv_results['test_accuracy']




# train test split
from sklearn.model_selection import train_test_split
x_train, x_test, y_train, y_test = train_test_split(x, y, test_size=0.2, random_state=0)
clf.fit(x_train, y_train)
y_pred = clf.predict(x_test)
# count each type of y_pred
pd.Series(y_pred).value_counts()
pd.Series(y).value_counts()
le.inverse_transform([3])


# confusion matrix
from sklearn.metrics import confusion_matrix
confusion_matrix(y_test, y_pred)








# predict on another cell line
df_pred=pd.read_csv(f"re_ci_hepg2.csv")
# remove "_hepg2" from region
df_pred["region_type"]=df_pred["region_type"].replace("_hepg2","",regex=True)
df_pred["constraint"]=df_pred["region_type"].apply(lambda x: x.split("_")[1])
df_pred["distal_proximal"]=df_pred["region_type"].apply(lambda x: x.split("_")[2])
df_pred["ti_ts"]=df_pred["region_type"].apply(lambda x: x.split("_")[3])
x_pred = df_pred[['codependent_tf_count', 'redundant_tf_count']]
y_truth = df_pred[pred_col]
y_truth = le.transform(y_truth)

y_pred=clf.predict(x_pred)
# how many are correctly predicted
sum(y_truth==y_pred)/len(y_truth)






























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

