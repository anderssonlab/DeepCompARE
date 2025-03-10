import pandas as pd
import seaborn as sns
from matplotlib import pyplot as plt
from adjustText import adjust_text

df=pd.read_csv('/isdata/alab/people/pcr980/DeepCompare/Pd9_TF_effect_and_constraint/tf_effect_and_constraints_promoters_hepg2.csv')
plt.figure(figsize=(10, 10))
sns.scatterplot(data=df, x='dstat_ism_cage_activity', y='dstat_ism_cage_probability', hue='z')
texts =[]
for i in range(len(df)):
    if df['dstat_ism_cage_activity'][i]*df['dstat_ism_cage_probability'][i]<0:
        texts.append(plt.text(df['dstat_ism_cage_activity'][i], df['dstat_ism_cage_probability'][i], df['protein'][i]))
    if abs(df['dstat_ism_cage_activity'][i]-df['dstat_ism_cage_probability'][i])>0.15:
       texts.append(plt.text(df['dstat_ism_cage_activity'][i], df['dstat_ism_cage_probability'][i], df['protein'][i]))

adjust_text(texts, arrowprops=dict(arrowstyle='-', color='black'))
plt.savefig('activity_vs_probability_promoters_hepg2.pdf')
plt.close()