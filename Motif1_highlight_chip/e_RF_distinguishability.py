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


def add_colocolization_info(df,cell_type):
    df_stripe=pd.read_csv(f'/isdata/alab/people/pcr980/Resource/ReMap2022/Colocolization/colocolization_potential_{cell_type}.csv')
    df=pd.merge(df, df_stripe, left_on="protein", right_on="tf2", how="inner")
    df.drop(columns=['tf2.1'], inplace=True)
    return df

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

def add_occupancy_info(df, cell_type):
    df_occupancy=pd.read_csv(f"/isdata/alab/people/pcr980/DeepCompare/Motif1_highlight_chip/summary_and_ks_test_promoters_{cell_type}.csv")
    df_occupancy=df_occupancy[['protein', 'conditional_occupancy']]
    df=pd.merge(df, df_occupancy, on="protein", how="inner")
    return df

def rf_summary(X,y,features):
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
    df=pd.melt(df, var_name="feature", value_name="importance")
    return df, accuracies


def whole_analysis(features,file_suffix,cell_type):
    
    # read in data
    df=pd.read_csv(f'/isdata/alab/people/pcr980/DeepCompare/Pd5_motif_info/motif_info_thresh_500_{file_suffix}.csv')
    df['feat_imp_orig_quartile']=pd.qcut(df['feat_imp_orig'], q=4, labels=["Q1", "Q2", "Q3", "Q4"])
    df=add_colocolization_info(df, cell_type)
    df=add_cpg_info(df,file_suffix)
    df=add_occupancy_info(df, cell_type.lower())

    X=df[features].values
    y_featimp=df['feat_imp_orig_quartile'].values
    y_chip=df['chip_evidence'].values

    # train random forest
    featimp_featimp, acc_featimp=rf_summary(X,y_featimp,features)
    featimp_chip, acc_chip=rf_summary(X,y_chip,features)
    df_rf_featimp=pd.concat([featimp_featimp.assign(type="feat_imp"), featimp_chip.assign(type="chip")])
    df_rf_featimp['feature']=df_rf_featimp['feature'].replace("score", "motif_score")

    # correlation between feature importance of feat_imp and chip
    mean_featimp_featimp=featimp_featimp.groupby("feature").mean().reset_index()
    mean_featimp_chip=featimp_chip.groupby("feature").mean().reset_index()
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
    plt.title(file_suffix)
    # TODO: change with/without_TF_annotation
    plt.savefig(f"feature_importance_{file_suffix}_no_TF_annotation.pdf")
    plt.close()



def main():
    features1=['score','count_all_TFs_no_thresh', 'count_TF_no_thresh','cpg_percentage',
        'seq_gc', 'context_gc_2bp','context_gc_10bp', 'context_gc_50bp', 'context_gc_100bp',
        'count_all_TFs_thresh_500', 'count_TF_thresh_500']
    features2=features1+['conditional_occupancy','colocolized_by_t2']
    dict_cell_type={"promoters_hepg2":"HepG2",
                    "enhancers_hepg2":"HepG2",
                    "promoters_k562":"K562",
                    "enhancers_k562":"K562"}
    for file_suffix in ["promoters_hepg2", "enhancers_hepg2", "promoters_k562", "enhancers_k562"]:
        logger.info(f"Running {file_suffix}")
        # TODO: change features1/2
        whole_analysis(features1,file_suffix,dict_cell_type[file_suffix])
    

if __name__ == "__main__":
    main()
    
    
# nohup python3 e_RF_distinguishability.py > e_RF_distinguishability.out &