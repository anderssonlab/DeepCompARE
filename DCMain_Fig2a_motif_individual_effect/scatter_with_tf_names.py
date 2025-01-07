import pandas as pd
import matplotlib.pyplot as plt
from adjustText import adjust_text


TRACK = "starr"

# Load datasets
hepg2_file = '/isdata/alab/people/pcr980/DeepCompare/Pd8_TF_effect_and_constraint/tf_effect_and_constraints_hepg2.csv'
k562_file = '/isdata/alab/people/pcr980/DeepCompare/Pd8_TF_effect_and_constraint/tf_effect_and_constraints_k562.csv'

hepg2_data = pd.read_csv(hepg2_file)
k562_data = pd.read_csv(k562_file)

# Merge the files by "protein" using an inner join
merged_data = pd.merge(hepg2_data, k562_data, on="protein", suffixes=("_hepg2", "_k562"), how="inner")

# Create a new column to classify points as HepG2-biased, K562-biased, or neutral
merged_data["bias"] = merged_data.apply(
    lambda row: "HepG2-biased" if row[f"dstat_isa_{TRACK}_activity_hepg2"] > row[f"dstat_isa_{TRACK}_activity_k562"] 
                else ("K562-biased" if row[f"dstat_isa_{TRACK}_activity_hepg2"] < row[f"dstat_isa_{TRACK}_activity_k562"] 
                      else "Neutral"),
    axis=1
)

# Calculate transparency based on the difference between x and y coordinates
merged_data["alpha"] = merged_data.apply(
    lambda row: min(1.0, max(0.2, abs(row[f"dstat_isa_{TRACK}_activity_hepg2"] - row[f"dstat_isa_{TRACK}_activity_k562"]))), axis=1
)

# Scatter plot with bias coloring and dynamic transparency
plt.figure(figsize=(6,6))
scatter = plt.scatter(
    x=merged_data[f"dstat_isa_{TRACK}_activity_hepg2"],
    y=merged_data[f"dstat_isa_{TRACK}_activity_k562"],
    c=merged_data["bias"].map({"HepG2-biased": "red", "K562-biased": "blue", "Neutral": "gray"}),
    alpha=merged_data["alpha"],
    s=50,  # Dot size
    edgecolor="none"
)

# Prepare text annotations for significant points
texts = []
for _, row in merged_data.iterrows():
    x = row[f"dstat_isa_{TRACK}_activity_hepg2"]
    y = row[f"dstat_isa_{TRACK}_activity_k562"]
    if x * y < 0 or abs(x + y) > 1:  # Opposite signs or large sum of coordinates
        texts.append(plt.text(x, y, row["protein"], fontsize=8, alpha=0.7))

# Adjust text to minimize overlaps
adjust_text(texts)

# Label the axes and title
plt.xlabel(f"Effect on HepG2 (dstat_isa_{TRACK}_activity)", fontsize=14)
plt.ylabel(f"Effect on K562 (dstat_isa_{TRACK}_activity)", fontsize=14)
plt.title(f"{TRACK} Effect on HepG2 vs. K562", fontsize=16)
plt.grid(False)  # Remove grid
plt.savefig(f"scatter_with_tf_names_{TRACK}.pdf")
plt.close()
