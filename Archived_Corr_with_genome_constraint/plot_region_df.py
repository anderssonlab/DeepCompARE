import pandas as pd
import seaborn as sns
import matplotlib.pyplot as plt

# Load data
df = pd.read_csv("regions_df.csv")
df=df.dropna().reset_index(drop=True)

sns.scatterplot(data=df, x="avg_zscore",y="avg_gradxinp",s=2,alpha=0.1)
plt.savefig("avg_zscore_vs_avg_gradxinp.png")
plt.close()

sns.scatterplot(data=df, x="avg_zscore",y="avg_ism",s=2,alpha=0.1)
plt.savefig("avg_zscore_vs_avg_ism.png")
plt.close()