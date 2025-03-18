import pandas as pd
from matplotlib import pyplot as plt
import seaborn as sns
from scipy.stats import mannwhitneyu
from loguru import logger

from adjustText import adjust_text


import sys
sys.path.insert(1,"/isdata/alab/people/pcr980/Scripts_python")
from tf_cooperativity import assign_cooperativity


import matplotlib
matplotlib.rcParams['pdf.fonttype']=42

# are linear TFs depleted in USF and PF?



def create_dimer_list(tf_list):
    dimer_list=[]
    for tf1 in tf_list:
        for tf2 in tf_list:
            dimer_list.append(f"{tf1}::{tf2}")
    return dimer_list


threshold_dict={"hepg2":{"pe":[0.3,0.7],"dhs":[0.48,0.78]},
                "k562":{"pe":[0.3,0.7],"dhs":[0.44,0.81]}}


for cell_line in ["hepg2","k562"]:
    for suffix in ["pe","dhs"]:
        for col in ["cooperativity_index","linearity_index"]:
            logger.info(f"Processing {cell_line}, {suffix}, {col}")
            df_tf=pd.read_csv(f"/isdata/alab/people/pcr980/DeepCompare/Pd7_TF_cooperativity/tf_cooperativity_index_{cell_line}_{suffix}.csv")
            df_tf=assign_cooperativity(df_tf,
                                       5,0.95,
                                       threshold_dict[cell_line][suffix][0],
                                       threshold_dict[cell_line][suffix][1])
            # make sure col does not contain nan
            df_tf=df_tf.dropna(subset=[col])
            # sort by cooperativity index
            df_tf["rank"]=df_tf[col].rank(ascending=True)
            #
            universal_stripe_factors=pd.read_csv("/isdata/alab/people/pcr980/Resource/universal_stripe_factors.txt",sep='\t').iloc[:,0].tolist()
            usf_dimers=create_dimer_list(universal_stripe_factors)
            universal_stripe_factors=universal_stripe_factors+usf_dimers
            pioneer_factors=pd.read_csv("/isdata/alab/people/pcr980/Resource/pioneer_factor_list.txt",sep='\t').iloc[:,0].tolist()
            pf_dimers=create_dimer_list(pioneer_factors)
            pioneer_factors=pioneer_factors+pf_dimers
            #
            df_usf=df_tf[df_tf["protein2"].isin(universal_stripe_factors)]
            df_pf=df_tf[df_tf["protein2"].isin(pioneer_factors)]
            #
            # mann-whitney u test
            df_non_usf=df_tf[~df_tf["protein2"].isin(universal_stripe_factors)]
            stat,p=mannwhitneyu(df_usf[col],df_non_usf[col])
            logger.info(f"USF vs non-USF: {p}")
            df_non_pf=df_tf[~df_tf["protein2"].isin(pioneer_factors)]
            stat,p=mannwhitneyu(df_pf[col],df_non_pf[col])
            logger.info(f"PF vs non-PF: {p}")

            plt.figure(figsize=(3.2,3.2))
            # thin frame
            plt.gca().spines['top'].set_linewidth(0.5)
            plt.gca().spines['right'].set_linewidth(0.5)
            plt.gca().spines['bottom'].set_linewidth(0.5)
            plt.gca().spines['left'].set_linewidth(0.5)
            sns.scatterplot(x="rank",y=col,data=df_tf,s=1,color='black')
            plt.xticks(fontsize=5)
            plt.yticks(fontsize=5)
            plt.xlabel("Rank")
            plt.ylabel("TF cooperativity index")
            plt.xlim(-70,df_tf.shape[0]+70)
            plt.ylim(df_tf[col].min()-0.1,df_tf[col].max()+0.1)
            texts=[]
            for i in range(df_usf.shape[0]):
                texts.append(plt.text(df_usf.iloc[i]["rank"],df_usf.iloc[i][col],df_usf.iloc[i]["protein2"],color='#4169E1',fontsize=5))

            for i in range(df_pf.shape[0]):
                texts.append(plt.text(df_pf.iloc[i]["rank"],df_pf.iloc[i][col],df_pf.iloc[i]["protein2"],color='darkorange',fontsize=5))

            adjust_text(texts)

            # add legend for text color
            plt.scatter([],[],color='#4169E1',label='Universal stripe factors')
            plt.scatter([],[],color='darkorange',label='Pioneer factors')
            plt.legend(fontsize=5,loc='lower right')
            # label 7
            plt.xlabel("Rank",fontsize=7)
            plt.ylabel(col,fontsize=7)
            plt.tight_layout()
            plt.savefig(f"usf_pf_{col}_{cell_line}_{suffix}.pdf",dpi=300)
            plt.close()

