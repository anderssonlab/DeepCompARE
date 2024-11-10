import numpy as np
import pandas as pd
from loguru import logger
from scipy.stats import pearsonr
import seaborn as sns
import matplotlib.pyplot as plt
import sys
sys.path.insert(1, "/isdata/alab/people/pcr980/Scripts_python/")

from seq_ops import generate_random_seqs
from in_silico_mutagenesis import compute_mutagenesis_score
from gradxinp import compute_gradxinp


#-------------------------------------
# ism vs isa
#-------------------------------------

# corr_list=[]
# for i in range(100):
#     seqs=generate_random_seqs(1000, 600)
#     # compute ism
#     df_isa=compute_mutagenesis_score(seqs,"isa")
#     df_ism=compute_mutagenesis_score(seqs,"ism")
#     # flatten the values and calculate correlation
#     isa = df_isa.values.flatten()
#     ism = df_ism.values.flatten()
#     corr, _ = pearsonr(isa, ism)
#     corr_list.append(corr)
#     logger.info(f"iteration {i}, correlation: {corr}")

# # output the correlation
# pd.DataFrame(corr_list).to_csv("corr_base_ism_vs_isa.csv", index=False)


df=pd.read_csv("corr_base_ism_vs_isa.csv")
df.columns=["corr"]
df["corr"].mean()
df["corr"].std()
# histogram
sns.histplot(data=df, x="corr", bins=20)
plt.xlabel("Pearson correlation")
plt.ylabel("Frequency")
plt.title("Pearson correlation between ISA and ISM")
plt.savefig("hist_corr_isa_ism.pdf")
plt.show()
plt.close()




#-------------------------------------
# |isa| vs |gradxinp|
#-------------------------------------


corr_list=[]
corr_abs_list=[]
for i in range(100):
    seqs=generate_random_seqs(1000, 600)
    # compute isa
    df_isa=compute_mutagenesis_score(seqs,"isa")
    df_gradxinp=compute_gradxinp(seqs)
    # flatten the values and calculate correlation
    isa = df_isa.values.flatten()
    gradxinp = df_gradxinp.values.flatten()
    corr, _ = pearsonr(isa, gradxinp)
    corr_list.append(corr)
    corr_abs, _ = pearsonr(np.abs(isa), np.abs(gradxinp))
    corr_abs_list.append(corr_abs)
    logger.info(f"iteration {i}, correlation: {corr}, correlation_abs: {corr_abs}")

# output the correlation
pd.DataFrame(corr_list).to_csv("corr_base_isa_vs_gradxinp.csv", index=False)
pd.DataFrame(corr_abs_list).to_csv("corr_abs_base_isa_vs_gradxinp.csv", index=False)





df=pd.read_csv("corr_base_isa_vs_gradxinp.csv")
df.columns=["corr"]
df["corr"].mean()
df["corr"].std()
# histogram
sns.histplot(data=df, x="corr", bins=20)
plt.xlabel("Pearson correlation")
plt.ylabel("Frequency")
plt.title("Pearson correlation between ISA and gradxinp")
plt.savefig("hist_corr_isa_gradxinp.pdf")
plt.close()




df=pd.read_csv("corr_abs_base_isa_vs_gradxinp.csv")
df.columns=["corr"]
df["corr"].mean()
df["corr"].std()
# histogram
sns.histplot(data=df, x="corr", bins=20)
plt.xlabel("Pearson correlation")
plt.ylabel("Frequency")
plt.title("Pearson correlation between ISA and |gradxinp|")
plt.savefig("hist_corr_isa_abs_gradxinp.pdf")
plt.close()






# nohup python3 per_base_importance.py > per_base_importance.log 2>&1 &