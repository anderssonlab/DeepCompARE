import pandas as pd


# read /isdata/alab/people/pcr980/DeepCompare/Pd8_TF_effect_and_constraint/tf_effect_and_constraints_enhancers_k562.csv
df=pd.read_csv('/isdata/alab/people/pcr980/DeepCompare/Pd8_TF_effect_and_constraint/tf_effect_and_constraints_enhancers_k562.csv')
# plot isa v.s z

import matplotlib.pyplot as plt
import seaborn as sns
sns.set(style="whitegrid")
# Draw a nested barplot to show survival
sns.scatterplot(x="avg_isa_sure_activity", y="z", data=df)
plt.savefig('avg_isa_sure_activity_v_z.png')
plt.close()