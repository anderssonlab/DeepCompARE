import matplotlib.pyplot as plt
import pandas as pd
import seaborn as sns


df = pd.read_csv("rf_summary.csv", index_col=0)

motif_imp_df = df[["file", "track_num", "mode", "acc_motif_imp"]].rename(columns={"acc_motif_imp": "acc"})
motif_imp_df["model"] = "motif_imp"
chip_df = df[["file", "track_num", "mode", "acc_chip"]].rename(columns={"acc_chip": "acc"})
chip_df["model"] = "chip"
transformed_df = pd.concat([motif_imp_df, chip_df])


di_df = transformed_df[transformed_df["mode"] == "di"].rename(columns={"acc": "acc_di"}).drop(columns=["mode"])
tri_df = transformed_df[transformed_df["mode"] == "tri"].rename(columns={"acc": "acc_tri"}).drop(columns=["mode"])
merged_df = pd.merge(di_df, tri_df, on=["file", "track_num", "model"])


# Creating the scatter plot
plt.figure(figsize=(6, 6))
sns.scatterplot(data=merged_df, x="acc_di", y="acc_tri", hue="model", style="model", s=100, palette="viridis")
min_val = min(merged_df['acc_di'].min(), merged_df['acc_tri'].min())
max_val = max(merged_df['acc_di'].max(), merged_df['acc_tri'].max())
plt.plot([min_val, max_val], [min_val, max_val], 'r--')
plt.title("Comparison of accuracy")
plt.xlabel("Accuracy of dinucleotide+#motif model")
plt.ylabel("Accuracy of trinucleotide model")
plt.savefig("rf_summary.pdf")
plt.close()