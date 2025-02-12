import pandas as pd
from matplotlib import pyplot as plt
import seaborn as sns
from scipy.stats import mannwhitneyu
from loguru import logger

from adjustText import adjust_text




import matplotlib
matplotlib.rcParams['pdf.fonttype']=42




universal_stripe_factors=pd.read_csv("/isdata/alab/people/pcr980/Resource/universal_stripe_factors.txt",sep='\t').iloc[:,0].tolist()

def create_dimer_list(tf_list):
    dimer_list=[]
    for tf1 in tf_list:
        for tf2 in tf_list:
            dimer_list.append(f"{tf1}::{tf2}")
    return dimer_list


for cell_line in ["hepg2","k562"]:
    for suffix in ["pe","dhs"]:
        
        logger.info(f"Processing {cell_line}")
        df_tf=pd.read_csv(f"/isdata/alab/people/pcr980/DeepCompare/Pd7_TF_cooperativity/tf_cooperativity_index_{cell_line}_{suffix}.csv")
        df_tf=df_tf[df_tf["c_sum"]>1].reset_index(drop=True)
        # sort by cooperativity index
        df_tf["rank"]=df_tf["cooperativity_index"].rank(ascending=True)

        universal_stripe_factors=pd.read_csv("/isdata/alab/people/pcr980/Resource/universal_stripe_factors.txt",sep='\t').iloc[:,0].tolist()
        usf_dimers=create_dimer_list(universal_stripe_factors)
        universal_stripe_factors=universal_stripe_factors+usf_dimers
        pioneer_factors=pd.read_csv("/isdata/alab/people/pcr980/Resource/pioneer_factor_list.txt",sep='\t').iloc[:,0].tolist()
        pf_dimers=create_dimer_list(pioneer_factors)
        pioneer_factors=pioneer_factors+pf_dimers
        
        df_usf=df_tf[df_tf["protein2"].isin(universal_stripe_factors)]
        df_pf=df_tf[df_tf["protein2"].isin(pioneer_factors)]

        # mann-whitney u test
        df_non_usf=df_tf[~df_tf["protein2"].isin(universal_stripe_factors)]
        stat,p=mannwhitneyu(df_usf["cooperativity_index"],df_non_usf["cooperativity_index"],alternative='less')
        logger.info(f"USF vs non-USF: {p}")
        df_non_pf=df_tf[~df_tf["protein2"].isin(pioneer_factors)]
        stat,p=mannwhitneyu(df_pf["cooperativity_index"],df_non_pf["cooperativity_index"],alternative='greater')
        logger.info(f"PF vs non-PF: {p}")

        plt.figure(figsize=(3.2,3.2))
        # thin frame
        plt.gca().spines['top'].set_linewidth(0.5)
        plt.gca().spines['right'].set_linewidth(0.5)
        plt.gca().spines['bottom'].set_linewidth(0.5)
        plt.gca().spines['left'].set_linewidth(0.5)
        sns.scatterplot(x="rank",y="cooperativity_index",data=df_tf,s=1,color='black')
        plt.xlim(-70,df_tf.shape[0]+70)
        plt.ylim(-0.1,1.1)
        plt.xlabel("Rank")
        plt.ylabel("TF cooperativity index")

        texts=[]
        for i in range(df_usf.shape[0]):
            texts.append(plt.text(df_usf.iloc[i]["rank"],df_usf.iloc[i]["cooperativity_index"],df_usf.iloc[i]["protein2"],color='#4169E1',fontsize=5))

        for i in range(df_pf.shape[0]):
            texts.append(plt.text(df_pf.iloc[i]["rank"],df_pf.iloc[i]["cooperativity_index"],df_pf.iloc[i]["protein2"],color='darkorange',fontsize=5))

        adjust_text(texts)

        # add legend for text color
        plt.scatter([],[],color='#4169E1',label='Universal stripe factors')
        plt.scatter([],[],color='darkorange',label='Pioneer factors')
        plt.legend(fontsize=5)
        # ticks 7
        plt.xticks(fontsize=5)
        plt.yticks(fontsize=5)
        # label 7
        plt.xlabel("Rank",fontsize=7)
        plt.ylabel("TF cooperativity index",fontsize=7)
        plt.tight_layout()
        plt.savefig(f"usf_pf_{cell_line}_{suffix}.pdf",dpi=300)
        plt.close()

