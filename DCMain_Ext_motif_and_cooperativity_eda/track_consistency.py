import pandas as pd
import seaborn as sns
import matplotlib.pyplot as plt
import sys
sys.path.insert(1,"/isdata/alab/people/pcr980/Scripts_python")
from tf_cooperativity import read_cooperativity

#------------------
# Step 1: get log
#------------------

thresholds=[0,0.0001,0.001,0.01,0.02,0.05,0.1,0.2,0.5]


def read_various_threholds(file_suffix, track_nums):
    df_shapes=[]
    for threshold in thresholds:
        df=read_cooperativity(f"/isdata/alab/people/pcr980/DeepCompare/Pd6_mutate_pair/mutate_pairs_{file_suffix}.csv",track_nums=track_nums,threshold=threshold)
        df_shapes.append(df.shape[0])
    return df_shapes

shapes_promoter_hepg2=read_various_threholds("promoters_hepg2",[0,2,4,6])
shapes_promoter_k562=read_various_threholds("promoters_k562",[1,3,5,len(thresholds)])
shapes_enhancer_hepg2=read_various_threholds("enhancers_hepg2",[0,2,4,6])
shapes_enhancer_k562=read_various_threholds("enhancers_k562",[1,3,5,len(thresholds)])


# output to csv
df=pd.DataFrame({"threshold":thresholds,
                 "promoter_hepg2":shapes_promoter_hepg2,
                 "promoter_k562":shapes_promoter_k562,
                 "enhancer_hepg2":shapes_enhancer_hepg2,
                 "enhancer_k562":shapes_enhancer_k562})
df.to_csv("tf_pairs_remain.csv",index=False)


#------------------
# Step 1: ananlyze log
#------------------

import re
log_data = open("track_consistency.out").read()
total_tf_pairs = re.findall(r"# Total tf pairs: (\d+)", log_data)
confusing_tf_pairs = re.findall(r"# Confusing tf pairs: (\d+)", log_data)

total_tf_pairs = list(map(int, total_tf_pairs))
confusing_tf_pairs = list(map(int, confusing_tf_pairs))

orig_shape_promoter_hepg2 = pd.read_csv("/isdata/alab/people/pcr980/DeepCompare/Pd6_mutate_pair/mutate_pairs_promoters_hepg2.csv").shape[0]
orig_shape_promoter_k562 = pd.read_csv("/isdata/alab/people/pcr980/DeepCompare/Pd6_mutate_pair/mutate_pairs_promoters_k562.csv").shape[0]
orig_shape_enhancer_hepg2 = pd.read_csv("/isdata/alab/people/pcr980/DeepCompare/Pd6_mutate_pair/mutate_pairs_enhancers_hepg2.csv").shape[0]
orig_shape_enhancer_k562 = pd.read_csv("/isdata/alab/people/pcr980/DeepCompare/Pd6_mutate_pair/mutate_pairs_enhancers_k562.csv").shape[0]

datasets=[[i]*len(thresholds) for i in ["promoter_hepg2", "promoter_k562", "enhancer_hepg2", "enhancer_k562"]]
datasets = [item for sublist in datasets for item in sublist]
df = pd.DataFrame({"tf_pairs_passing_threshold": total_tf_pairs,
                   "confusing_tf_pairs": confusing_tf_pairs,
                   "dataset": datasets,
                   "threshold": thresholds*4,
                   "orig_shape": [orig_shape_promoter_hepg2]*len(thresholds) + [orig_shape_promoter_k562]*len(thresholds) + [orig_shape_enhancer_hepg2]*len(thresholds) + [orig_shape_enhancer_k562]*len(thresholds)})

df["tf_pairs_passing_threshold"]=df["tf_pairs_passing_threshold"]/df["orig_shape"]
df["confusing_tf_pairs"]=df["confusing_tf_pairs"]/df["orig_shape"]

# seaborn line plot, hue=dataset
plt.figure(figsize=(8,6))
sns.lineplot(x="threshold", y="tf_pairs_passing_threshold", hue="dataset", data=df, linestyle="--",legend=False)
sns.lineplot(x="threshold", y="confusing_tf_pairs", hue="dataset", data=df, linestyle="-")
# Customize the title, labels, and scale
plt.title("Cooperativity consistency across tracks")
plt.xscale("log")
plt.xlabel("Threshold")
plt.ylabel("Fraction of TF pairs")
# Handle the legend
# Get the hue legend entries
handles, labels = plt.gca().get_legend_handles_labels()
datasets = df['dataset'].unique()
# Create custom legend entries for line styles
custom_lines = [
    plt.Line2D([0], [0], color="black", linestyle="--", label="% TF pairs passing threshold"),
    plt.Line2D([0], [0], color="black", linestyle="-", label="% TF pairs with contradiction across tracks"),
]
# Combine dataset (hue) and line style legends
plt.legend(
    handles + custom_lines,
    labels + ["% TF pairs passing threshold", "% TF pairs with contradiction across tracks"],
    loc="best"
)
# Save and close the plot
plt.savefig("track_consistency.pdf")
plt.close()


# nohup python3 track_consistency.py > track_consistency.out &
