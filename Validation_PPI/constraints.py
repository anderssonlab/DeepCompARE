import pandas as pd

df_med=pd.read_csv('2024-08-07_MED-TF_interactions.txt', sep='\t')
# select only significant interactions
# df_med=df_med[df_med['significant']].reset_index(drop=True)

df_tf=pd.read_csv('/isdata/alab/people/pcr980/DeepCompare/Pd8_TF_effect_and_constraint/tf_effect_and_constraints_distal_k562.csv')
# add a column to indicate if the TF is in df_med
df_tf['in_MED']=df_tf['protein'].isin(df_med['gene']).astype(int)




# plot distribution of 
# 'avg_isa_cage_activity', 'avg_isa_dhs_activity','avg_isa_starr_activity', 'avg_isa_sure_activity','avg_isa_cage_probability', 'avg_isa_dhs_probability',
# 'avg_isa_starr_probability', 'avg_isa_sure_probability', 'dstat_isa_cage_activity', 'dstat_isa_dhs_activity', 'dstat_isa_starr_activity', 'dstat_isa_sure_activity',
# 'dstat_isa_cage_probability', 'dstat_isa_dhs_probability',
# 'dstat_isa_starr_probability', 'dstat_isa_sure_probability', 'z'
# split by "in_MED"

import matplotlib.pyplot as plt
import seaborn as sns


for col in df_tf.columns[1:]:
    plt.figure()
    sns.kdeplot(df_tf, x=col, hue='in_MED', fill=True, common_norm=False, common_grid=True)
    plt.title(col)
    plt.savefig(f'{col}.png')   
    plt.close()