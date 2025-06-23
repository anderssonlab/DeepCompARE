import pandas as pd
from scipy.stats import pearsonr
import matplotlib.pyplot as plt
import seaborn as sns

df = pd.read_csv('motif_ism_isa.csv')
# shuffle
df = df.sample(frac=1).reset_index(drop=True)


# select 100 rows every time until all rows are selected
# calculate pearsom r between 'ism_score' and 'isa_track0'
corrs= []
for i in range(df.shape[0] // 100):
    df_subset = df.iloc[i * 100:(i + 1) * 100]
    r, _ = pearsonr(df_subset['ism_score'], df_subset['isa_track0'])
    corrs.append(r)

# for corrs, calculate mean and 95% confidence interval
mean_corr = sum(corrs) / len(corrs)
conf_interval = 1.96 * (pd.Series(corrs).std() / (len(corrs) ** 0.5))
print(f"Mean correlation: {mean_corr:.4f}, 95% CI: {conf_interval:.4f}")


plt.figure(figsize=(2,2))
# thin frame
plt.gca().spines['top'].set_visible(False)
plt.gca().spines['bottom'].set_linewidth(0.5)
plt.gca().spines['left'].set_linewidth(0.5)
plt.gca().spines['right'].set_visible(False)
sns.kdeplot(corrs, linewidth=1)
plt.title("Motif level ISM vs ISA", fontsize=7)
plt.xlabel("Pearson correlation",fontsize=7)
plt.ylabel("Density",fontsize=7)
plt.xticks(fontsize=7)
plt.yticks(fontsize=7)
plt.xlim(0,1)
plt.tight_layout()
plt.savefig(f"corr_motif_level_ism_isa.pdf")
plt.close()
