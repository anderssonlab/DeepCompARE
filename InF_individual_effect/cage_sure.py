import pandas as pd
import matplotlib.pyplot as plt
import seaborn as sns


import pandas as pd
import seaborn as sns
import matplotlib.pyplot as plt
from scipy.cluster.hierarchy import linkage


for file_suffix in ["promoters_hepg2", "promoters_k562", "enhancers_hepg2", "enhancers_k562"]:
    df = pd.read_csv(f"/isdata/alab/people/pcr980/DeepCompare/Pd9_TF_effect_and_constraint/tf_effect_and_constraints_{file_suffix}.csv")
    df_sub = df[["protein", f"dstat_ism_cage_activity", f"dstat_ism_sure_activity",
                    f"dstat_ism_cage_probability", f"dstat_ism_sure_probability"]].copy()
    # set index to protein
    df_sub.set_index("protein", inplace=True)
    # replace NaN with 0
    df_sub.fillna(0, inplace=True)
    # Perform hierarchical clustering on the standardized data
    row_clusters = linkage(df_sub, method='average', metric='euclidean')
    col_clusters = linkage(df_sub.T, method='average', metric='euclidean')
    # Create a clustermap based on the row linkage and column linkage
    sns.clustermap(df_sub, row_linkage=row_clusters, col_linkage=col_clusters, figsize=(3, 10), annot=False)
    # Save the clustermap
    plt.title(f"CAGE vs. SuRE for {file_suffix}")
    plt.tight_layout()
    plt.savefig(f"cage_sure_clustering_{file_suffix}.pdf")
    plt.close()


