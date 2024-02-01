import matplotlib.pyplot as plt
import seaborn as sns
from adjustText import adjust_text
from scipy.cluster.hierarchy import linkage, dendrogram


def scatter_plot_with_annotation(df, x, y, text_col, file_name, xlab, ylab):
    texts = []
    plt.figure(figsize=(10, 8))  
    scatter_plot = sns.scatterplot(data=df, x=x, y=y)
    for line in range(0, df.shape[0]):
        text=scatter_plot.text(df[x][line], df[y][line], 
                                df[text_col][line], horizontalalignment='left', 
                                size='xx-small', color='black')
        texts.append(text)
    adjust_text(texts)
    plt.xlabel(xlab)
    plt.ylabel(ylab)
    plt.savefig(file_name)
    plt.close()
    
    
    
def heatmap_plot(df,x,y,color,out_path):
    pivot_df = df.pivot(x, y, color)
    pivot_df = pivot_df.fillna(0)  
    # Clustering
    row_clusters = linkage(pivot_df, method='ward', optimal_ordering=True)
    col_clusters = linkage(pivot_df.T, method='ward', optimal_ordering=True)
    # Creating a clustered heatmap
    sns.clustermap(pivot_df, row_linkage=row_clusters, col_linkage=col_clusters, cmap="viridis")
    plt.savefig(out_path)
    plt.close()