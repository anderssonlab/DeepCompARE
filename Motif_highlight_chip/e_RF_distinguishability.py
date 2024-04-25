import pandas as pd
from sklearn.ensemble import RandomForestClassifier
from sklearn.model_selection import KFold
from sklearn.metrics import accuracy_score
from scipy.stats import pearsonr
import matplotlib.pyplot as plt
import numpy as np
import pyranges as pr
import seaborn as sns
from loguru import logger

import sys
sys.path.append("/isdata/alab/people/pcr980/DeepCompare/Scripts_python")
from seq_ops import SeqExtractor, dinucleotide_frequencies
seq_extractor=SeqExtractor()

# TF level annotation
def add_colocolization_info(df,cell_type):
    df_stripe=pd.read_csv(f'/isdata/alab/people/pcr980/Resource/ReMap2022/Colocolization/colocolization_potential_{cell_type}.csv')
    df=pd.merge(df, df_stripe, left_on="protein", right_on="tf2", how="inner")
    df.drop(columns=['tf2.1'], inplace=True)
    return df

def add_occupancy_info(df, cell_type):
    df_occupancy=pd.read_csv(f"/isdata/alab/people/pcr980/DeepCompare/Motif_highlight_chip/summary_and_ks_test_promoters_{cell_type}.csv")
    df_occupancy=df_occupancy[['protein', 'conditional_occupancy']]
    df=pd.merge(df, df_occupancy, on="protein", how="inner")
    return df

# Sequence level annotation
def add_cpg_info(df,bed_file):
    df_cpg=pd.read_csv('/isdata/alab/people/pcr980/Resource/hg38_cpgIslandExtUnmasked.txt', sep='\t', header=None)
    df_cpg.columns=['bin','chrom', 'chromStart', 'chromEnd', 'name', 'length', 'cpgNum', 'gcNum', 'perCpg', 'perGc', 'obsExp']
    gr_cpg = pr.PyRanges(chromosomes=df_cpg["chrom"],starts=df_cpg["chromStart"],ends=df_cpg["chromEnd"])
    # read sequence bed file to gr
    df_seq=pd.read_csv(f"/isdata/alab/people/pcr980/DeepCompare/Pd4_promoters_enhancers_and_featimp/{bed_file}.bed", header=None,sep='\t')
    df_seq=df_seq.iloc[:,[0,1,2]]
    df_seq.columns=['Chromosome', 'Start', 'End']
    df_seq["seq_idx"] = "Seq"+df_seq.index.astype(str)
    gr_seq=pr.PyRanges(df_seq)
    # overlap
    overlaps = gr_seq.join(gr_cpg).df
    overlaps['overlap_length'] = overlaps.apply(lambda row: min(row['End'], row['End_b']) - max(row['Start'], row['Start_b']),axis=1)
    overlaps['cpg_percentage']=overlaps['overlap_length']/600
    overlaps=overlaps[["seq_idx","cpg_percentage"]]
    # merge with df
    df=pd.merge(df, overlaps, left_on="seq_idx", right_on="seq_idx", how="left")
    df['cpg_percentage'].fillna(0, inplace=True)
    return df

def add_dinucleotide(df,bed_file):
    # convert bed file to sequence
    df_seq=pd.read_csv(f"/isdata/alab/people/pcr980/DeepCompare/Pd4_promoters_enhancers_and_featimp/{bed_file}.bed", header=None,sep='\t')
    df_seq=df_seq.iloc[:,[0,1,2]]
    df_seq["seq"]=df_seq.apply(lambda row: seq_extractor.get_seq(row[0],row[1],row[2]), axis=1)
    df_difreq=dinucleotide_frequencies(df_seq["seq"])
    df_difreq["seq_idx"]="Seq"+df_difreq.index.astype(str)
    df=pd.merge(df, df_difreq, left_on="seq_idx", right_on="seq_idx", how="left")
    return df

def rf_summary(X,y,features):
    features_dinucleotide=[a+b for a in 'ACGT' for b in 'ACGT']
    features_trinucleotide=[a+b+c for a in 'ACGT' for b in 'ACGT' for c in 'ACGT']
    kf = KFold(n_splits=5, shuffle=True, random_state=42)
    all_feature_importances = []
    accuracies = []
    
    for train_index, test_index in kf.split(X):
        X_train, X_test = X[train_index], X[test_index]
        y_train, y_test = y[train_index], y[test_index]
        clf = RandomForestClassifier()
        clf.fit(X_train, y_train)
        all_feature_importances.append(clf.feature_importances_)
        accuracies.append(accuracy_score(y_test, clf.predict(X_test)))
    all_feature_importances = np.array(all_feature_importances)
    df=pd.DataFrame(all_feature_importances)
    df.columns=features
    df_unsummarized=df.copy()
    # summarize features
    df.index=["fold1","fold2","fold3","fold4","fold5"]
    max_di=df[features_dinucleotide].max(axis=0).sort_values(ascending=False).index[0]
    min_di=df[features_dinucleotide].min(axis=0).sort_values(ascending=True).index[0]
    df["dinucleotide_mean"]=df[features_dinucleotide].mean(axis=1)
    df[f"dinucleotide_max({max_di})"]=df[features_dinucleotide].max(axis=1)
    df[f"dinucleotide_min({min_di})"]=df[features_dinucleotide].min(axis=1)
    df.drop(columns=features_dinucleotide, inplace=True)
    
    df=pd.melt(df, var_name="feature", value_name="importance", ignore_index=False).reset_index()
    df_unsummarized=pd.melt(df_unsummarized, var_name="feature", value_name="importance")
    df.rename(columns={"index":"fold"}, inplace=True)
    return df, df_unsummarized, accuracies


def whole_analysis(features,file_suffix):    
    # read in data
    df=pd.read_csv(f'/isdata/alab/people/pcr980/DeepCompare/Pd5_motif_info/motif_info_thresh_500_{file_suffix}.csv')
    df['feat_imp_orig_quartile']=pd.qcut(df['feat_imp_orig'], q=4, labels=["Q1", "Q2", "Q3", "Q4"])
    df=add_cpg_info(df,file_suffix)
    df=add_dinucleotide(df,file_suffix)
    # df=add_trinucleotide(df,file_suffix)
    
    X=df[features].values
    y_featimp=df['feat_imp_orig_quartile'].values
    y_chip=df['chip_evidence'].values

    # train random forest
    logger.info("Training to predict feature importance quartile")
    featimp_featimp, featimp_featimp_unsummarized, acc_featimp=rf_summary(X,y_featimp,features)
    logger.info("Training to predict chip evidence")
    featimp_chip, featimp_chip_unsummarized,acc_chip=rf_summary(X,y_chip,features)
    df_rf_featimp=pd.concat([featimp_featimp.assign(type="feat_imp"), featimp_chip.assign(type="chip")])
    df_rf_featimp['feature']=df_rf_featimp['feature'].replace("score", "motif_score")
    df_rf_featimp.to_csv(f"feature_importance_{file_suffix}.csv")
    
    # correlation between feature importance of feat_imp and chip
    mean_featimp_featimp=featimp_featimp_unsummarized.groupby("feature").mean().reset_index()
    mean_featimp_chip=featimp_chip_unsummarized.groupby("feature").mean().reset_index()
    mean_featimp=pd.merge(mean_featimp_featimp, mean_featimp_chip, on="feature", suffixes=("_featimp", "_chip"))
    corr,p=pearsonr(mean_featimp["importance_featimp"], mean_featimp["importance_chip"])
    
    # plot
    df_rf_featimp['feature'] = pd.Categorical(df_rf_featimp['feature'], categories=df_rf_featimp[df_rf_featimp.type=="feat_imp"].sort_values("importance", ascending=False).feature.unique())
    sns.barplot(data=df_rf_featimp, x="feature", y="importance", hue="type",errorbar=('ci', 95),capsize=0,errwidth=1)
    plt.text(0.3, 0.95, f"Mean acc (feat_imp): {round(np.mean(acc_featimp),2)}", transform=plt.gca().transAxes)
    plt.text(0.3, 0.88, f"Mean acc (chip): {round(np.mean(acc_chip),2)}", transform=plt.gca().transAxes)
    plt.text(0.3, 0.81, f"r={round(corr,2)}, p={round(p,2)}", transform=plt.gca().transAxes)
    plt.xticks(rotation=90)
    plt.subplots_adjust(bottom=0.4)
    # put legend at top right
    plt.legend(loc='upper right')
    plt.title(file_suffix)
    plt.savefig(f"feature_importance_{file_suffix}.pdf")
    plt.close()



def main():
    features=['score','count_all_TFs_no_thresh', 'count_TF_no_thresh','cpg_percentage',
        'seq_gc', 'context_gc_2bp','context_gc_10bp', 'context_gc_50bp', 'context_gc_100bp',
        'count_all_TFs_thresh_500', 'count_TF_thresh_500']
    features_dinucleotide=[a+b for a in 'ACGT' for b in 'ACGT']
    #features_trinucleotide=[a+b+c for a in 'ACGT' for b in 'ACGT' for c in 'ACGT']
    features=features+features_dinucleotide   #+features_trinucleotide
    # for file_suffix in ["promoters_hepg2", "enhancers_hepg2", "promoters_k562", "enhancers_k562"]:
    for file_suffix in ["enhancers_hepg2"]:
        logger.info(f"Running {file_suffix}")
        whole_analysis(features,file_suffix)
    

if __name__ == "__main__":
    main()
    
    
# nohup python3 e_RF_distinguishability.py > e_RF_distinguishability.out &