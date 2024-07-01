import pandas as pd
import matplotlib.pyplot as plt
import seaborn as sns

df=pd.read_csv("counts2.csv")
df.drop(columns=['count_neg_ism', 'count_pos_ism'], inplace=True)

df["phylop_missing_rate"]=1-df["count_base_remain"]/df["count_base_total"]

df["cons_acc_ratio"]=df["count_conserved"]/df["count_accelerated"]
df["sig_cons_acc_ratio"]=df["count_sig_conserved"]/df["count_sig_accelerated"]
df["sig_ratio_cons"]=df["count_sig_conserved"]/df["count_conserved"]
df["sig_ratio_acc"]=df["count_sig_accelerated"]/df["count_accelerated"]


# plot phylop_missing_rate
plt.figure(figsize=(8,6))
sns.scatterplot(x="file_prefix", y="phylop_missing_rate", hue="phylop", data=df)
plt.xticks(rotation=30)
plt.title("#bases without phyloP score/#total bases")
plt.subplots_adjust(bottom=0.2)
plt.savefig("phylop_missing_rate_by_file.png")
plt.close()

plt.figure(figsize=(8,6))
sns.scatterplot(x="phylop", y="phylop_missing_rate", hue="file_prefix", data=df)
plt.xticks(rotation=30)
plt.title("#bases without phyloP score/#total bases")
plt.subplots_adjust(bottom=0.2)
plt.savefig("phylop_missing_rate_by_phylop.png")
plt.close()



# plot cons_acc_ratio
plt.figure(figsize=(8,6))
sns.scatterplot(x="file_prefix", y="cons_acc_ratio", hue="phylop", data=df)
plt.xticks(rotation=30)
plt.title("#conserved/#accelerated")
plt.subplots_adjust(bottom=0.2)
plt.savefig("cons_acc_ratio_by_file.png")
plt.close()

plt.figure(figsize=(8,6))
sns.scatterplot(x="phylop", y="cons_acc_ratio", hue="file_prefix", data=df)
plt.xticks(rotation=30)
plt.title("#conserved/#accelerated")
plt.subplots_adjust(bottom=0.2)
plt.savefig("cons_acc_ratio_by_phylop.png")
plt.close()



# plot sig_cons_acc_ratio
plt.figure(figsize=(8,6))
sns.scatterplot(x="file_prefix", y="sig_cons_acc_ratio", hue="phylop", data=df)
plt.xticks(rotation=30)
plt.title("#significantly_conserved/#significantly_accelerated")
plt.subplots_adjust(bottom=0.2)
plt.savefig("sig_cons_acc_ratio_by_file.png")
plt.close()

plt.figure(figsize=(8,6))
sns.scatterplot(x="phylop", y="sig_cons_acc_ratio", hue="file_prefix", data=df)
plt.xticks(rotation=30)
plt.title("#significantly_conserved/#significantly_accelerated")
plt.subplots_adjust(bottom=0.2)
# put legend at bottom left corner
plt.legend(loc='lower left')
plt.savefig("sig_cons_acc_ratio_by_phylop.png")
plt.close()




# plot sig_ratio_cons
plt.figure(figsize=(8,6))
sns.scatterplot(x="file_prefix", y="sig_ratio_cons", hue="phylop", data=df)
plt.xticks(rotation=30)
plt.title("#significantly_conserved/#conserved")
plt.subplots_adjust(bottom=0.2)
plt.savefig("sig_ratio_cons_by_file.png")
plt.close()

plt.figure(figsize=(8,6))
sns.scatterplot(x="phylop", y="sig_ratio_cons", hue="file_prefix", data=df)
plt.xticks(rotation=30)
plt.title("#significantly_conserved/#conserved")
plt.subplots_adjust(bottom=0.2)
plt.savefig("sig_ratio_cons_by_phylop.png")
plt.close()





# plot sig_ratio_acc
plt.figure(figsize=(8,6))
sns.scatterplot(x="file_prefix", y="sig_ratio_acc", hue="phylop", data=df)
plt.xticks(rotation=30)
plt.title("#significantly_accelerated/#accelerated")
plt.subplots_adjust(bottom=0.2)
plt.savefig("sig_ratio_acc_by_file.png")
plt.close()

plt.figure(figsize=(8,6))
sns.scatterplot(x="phylop", y="sig_ratio_acc", hue="file_prefix", data=df)
plt.xticks(rotation=30)
plt.title("#significantly_accelerated/#accelerated")
plt.subplots_adjust(bottom=0.2)
plt.savefig("sig_ratio_acc_by_phylop.png")
plt.close()
