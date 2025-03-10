import pandas as pd
import matplotlib.pyplot as plt
import seaborn as sns

from scipy.stats import mannwhitneyu
from sklearn.decomposition import PCA


df_med=pd.read_csv('2024-08-07_MED-TF_interactions.txt', sep='\t')
# select only significant interactions
df_med=df_med[df_med['significant']].reset_index(drop=True)

df_tf=pd.read_csv('/isdata/alab/people/pcr980/DeepCompare/Pd8_TF_effect_and_constraint/tf_effect_and_constraints_distal_k562.csv')
# add a column to indicate if the TF is in df_med
df_tf['in_MED']=df_tf['protein'].isin(df_med['gene']).astype(int)


# calculate activity_probability difference between d_stat
df_tf['activity_probability_diff_cage']=df_tf['dstat_isa_cage_activity']-df_tf['avg_isa_cage_activity']
df_tf['activity_probability_diff_dhs']=df_tf['dstat_isa_dhs_activity']-df_tf['avg_isa_dhs_activity']
df_tf['activity_probability_diff_starr']=df_tf['dstat_isa_starr_activity']-df_tf['avg_isa_starr_activity']
df_tf['activity_probability_diff_sure']=df_tf['dstat_isa_sure_activity']-df_tf['avg_isa_sure_activity']

# calculate cage_sure diff
df_tf['cage_sure_diff']=df_tf['dstat_isa_cage_activity']-df_tf['dstat_isa_sure_activity']


# plot distribution of 
# 'avg_isa_cage_activity', 'avg_isa_dhs_activity','avg_isa_starr_activity', 'avg_isa_sure_activity','avg_isa_cage_probability', 'avg_isa_dhs_probability',
# 'avg_isa_starr_probability', 'avg_isa_sure_probability', 'dstat_isa_cage_activity', 'dstat_isa_dhs_activity', 'dstat_isa_starr_activity', 'dstat_isa_sure_activity',
# 'dstat_isa_cage_probability', 'dstat_isa_dhs_probability',
# 'dstat_isa_starr_probability', 'dstat_isa_sure_probability', 'z',
# 'activity_probability_diff_cage', 'activity_probability_diff_dhs', 'activity_probability_diff_starr', 'activity_probability_diff_sure'
# split by "in_MED"




for col in ['avg_isa_cage_activity', 'avg_isa_dhs_activity',
       'avg_isa_starr_activity', 'avg_isa_sure_activity',
       'avg_isa_cage_probability', 'avg_isa_dhs_probability',
       'avg_isa_starr_probability', 'avg_isa_sure_probability',
       'dstat_isa_cage_activity', 'dstat_isa_dhs_activity',
       'dstat_isa_starr_activity', 'dstat_isa_sure_activity',
       'dstat_isa_cage_probability', 'dstat_isa_dhs_probability',
       'dstat_isa_starr_probability', 'dstat_isa_sure_probability',
       'num_expected_variants', 'num_rare_variants', 'oe', 'chisq', 'z', 'activity_probability_diff_cage',
       'activity_probability_diff_dhs', 'activity_probability_diff_starr',
       'activity_probability_diff_sure', 'cage_sure_diff']:
    plt.figure()
    sns.kdeplot(df_tf, x=col, hue='in_MED', fill=True, common_norm=False, common_grid=True)
    # mannwhitneyu test
    u, p=mannwhitneyu(df_tf[df_tf['in_MED']==0][col], df_tf[df_tf['in_MED']==1][col])
    plt.text(0.5, 0.5, f'U={u}\np={p}', horizontalalignment='center', verticalalignment='center', transform=plt.gca().transAxes)
    plt.title(col)
    plt.savefig(f'Plots/{col}.png')   
    plt.close()
    

# pca on 'avg_isa_cage_activity', 'avg_isa_dhs_activity','avg_isa_starr_activity', 'avg_isa_sure_activity','avg_isa_cage_probability', 'avg_isa_dhs_probability','avg_isa_starr_probability', 'avg_isa_sure_probability',
# color by in_MED

for aggregation_method in ["avg","dstat"]:
    pca=PCA(n_components=2)
    X=df_tf[[f'{aggregation_method}_isa_cage_activity', f'{aggregation_method}_isa_dhs_activity',f'{aggregation_method}_isa_starr_activity', f'{aggregation_method}_isa_sure_activity',f'{aggregation_method}_isa_cage_probability', f'{aggregation_method}_isa_dhs_probability',f'{aggregation_method}_isa_starr_probability', f'{aggregation_method}_isa_sure_probability']]
    X_pca=pca.fit_transform(X)
    plt.figure()
    plt.scatter(X_pca[:,0], X_pca[:,1], c=df_tf['in_MED'],s=10)
    plt.title('PCA')
    plt.legend()
    plt.savefig(f'Plots/PCA_{aggregation_method}.png')
    plt.close()