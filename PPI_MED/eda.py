import pandas as pd
import matplotlib.pyplot as plt
import seaborn as sns

from scipy.stats import mannwhitneyu
from sklearn.decomposition import PCA


mode="dhs"



df_med=pd.read_csv('2024-08-07_MED-TF_interactions.txt', sep='\t')
# select only significant interactions
df_med=df_med[df_med['significant']].reset_index(drop=True)

df_tf=pd.read_csv(f'/isdata/alab/people/pcr980/DeepCompare/Pd8_TF_effect_and_constraint/tf_effect_and_constraints_k562_{mode}.csv')
# subset for proper background TFs
proteins_background=pd.read_csv('MED_experiment_K562proteomics.txt', sep='\t',header=None)[0].tolist()
df_tf=df_tf[df_tf['protein'].isin(proteins_background)].reset_index(drop=True)
# add a column to indicate if the TF is in df_med
df_tf['in_MED']=df_tf['protein'].isin(df_med['gene']).astype(int)



for col in ['avg_isa_cage_activity', 'avg_isa_dhs_activity',
       'avg_isa_starr_activity', 'avg_isa_sure_activity',
       'avg_isa_cage_probability', 'avg_isa_dhs_probability',
       'avg_isa_starr_probability', 'avg_isa_sure_probability',
       'dstat_isa_cage_activity', 'dstat_isa_dhs_activity',
       'dstat_isa_starr_activity', 'dstat_isa_sure_activity',
       'dstat_isa_cage_probability', 'dstat_isa_dhs_probability',
       'dstat_isa_starr_probability', 'dstat_isa_sure_probability',
       'num_expected_variants', 'num_rare_variants', 'oe', 'chisq', 'z']:
    plt.figure()
    sns.kdeplot(df_tf, x=col, hue='in_MED', fill=True, common_norm=False, common_grid=True)
    # mannwhitneyu test
    u, p=mannwhitneyu(df_tf[df_tf['in_MED']==0][col], df_tf[df_tf['in_MED']==1][col])
    plt.text(0.5, 0.5, f'U={u}\np={p}', horizontalalignment='center', verticalalignment='center', transform=plt.gca().transAxes)
    plt.title(col)
    plt.savefig(f'Plots/{col}_{mode}.png')   
    plt.close()
    







# pca on 'avg_isa_cage_activity', 'avg_isa_dhs_activity','avg_isa_starr_activity', 'avg_isa_sure_activity','avg_isa_cage_probability', 'avg_isa_dhs_probability','avg_isa_starr_probability', 'avg_isa_sure_probability',
# color by in_MED

# for aggregation_method in ["avg","dstat"]:
#     pca=PCA(n_components=2)
#     X=df_tf[[f'{aggregation_method}_isa_cage_activity', f'{aggregation_method}_isa_dhs_activity',f'{aggregation_method}_isa_starr_activity', f'{aggregation_method}_isa_sure_activity',f'{aggregation_method}_isa_cage_probability', f'{aggregation_method}_isa_dhs_probability',f'{aggregation_method}_isa_starr_probability', f'{aggregation_method}_isa_sure_probability']]
#     X_pca=pca.fit_transform(X)
#     plt.figure()
#     plt.scatter(X_pca[:,0], X_pca[:,1], c=df_tf['in_MED'],s=10)
#     plt.title('PCA')
#     plt.legend()
#     plt.savefig(f'Plots/PCA_{aggregation_method}.png')
#     plt.close()