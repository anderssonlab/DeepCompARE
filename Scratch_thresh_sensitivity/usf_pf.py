import pandas as pd
from matplotlib import pyplot as plt
import seaborn as sns
from scipy.stats import mannwhitneyu
from loguru import logger

from adjustText import adjust_text


for cell_line in ["hepg2","k562"]:
    for tpm_thresh in [0.1,0.5,1.0]:
        for nonlinearity_thresh in [0.01,0.05,0.1]:
            logger.info(f"Processing {cell_line} {tpm_thresh} {nonlinearity_thresh}")
            df_tf=pd.read_csv(f"Pd1_cis/tf_cooperativity_index_{tpm_thresh}_{nonlinearity_thresh}_{cell_line}.csv")
            df_tf=df_tf[df_tf["c_sum"]>5].reset_index(drop=True)
            # sort by cooperativity index
            df_tf["rank"]=df_tf["cooperativity_index"].rank(ascending=True)

            universal_stripe_factors=pd.read_csv("/isdata/alab/people/pcr980/Resource/universal_stripe_factors.txt",sep='\t').iloc[:,0].tolist()
            pioneer_factors=pd.read_csv("/isdata/alab/people/pcr980/Resource/pioneer_factor_list.txt",sep='\t').iloc[:,0].tolist()
            df_usf=df_tf[df_tf["protein2"].isin(universal_stripe_factors)]
            df_pf=df_tf[df_tf["protein2"].isin(pioneer_factors)]

            sns.scatterplot(x="rank",y="cooperativity_index",data=df_tf,s=5,color='black')
            plt.xlim(-40,df_tf.shape[0]+40)
            plt.ylim(-0.1,1.1)
            plt.xlabel("Rank")
            plt.ylabel("TF cooperativity index")

            texts=[]
            for i in range(df_usf.shape[0]):
                texts.append(plt.text(df_usf.iloc[i]["rank"],df_usf.iloc[i]["cooperativity_index"],df_usf.iloc[i]["protein2"],color='#4169E1'))

            for i in range(df_pf.shape[0]):
                texts.append(plt.text(df_pf.iloc[i]["rank"],df_pf.iloc[i]["cooperativity_index"],df_pf.iloc[i]["protein2"],color='darkorange'))

            adjust_text(texts)

            # add legend for text color
            plt.scatter([],[],color='#4169E1',label='Universal stripe factors')
            plt.scatter([],[],color='darkorange',label='Pioneer factors')
            plt.legend()
            
            # mann-whitney u test
            df_non_usf=df_tf[~df_tf["protein2"].isin(universal_stripe_factors)]
            _,p_usf=mannwhitneyu(df_usf["cooperativity_index"],df_non_usf["cooperativity_index"],alternative='less')
            df_non_pf=df_tf[~df_tf["protein2"].isin(pioneer_factors)]
            _,p_pf=mannwhitneyu(df_pf["cooperativity_index"],df_non_pf["cooperativity_index"],alternative='greater')
            plt.title(f"tpm_thresh={tpm_thresh}, nonlinearity_thresh={nonlinearity_thresh}, cell_line={cell_line}\nUSF vs non-USF p={p_usf:.2e}, PF vs non-PF p={p_pf:.2e}")

            plt.savefig(f"Plots/usf_pf_{tpm_thresh}_{nonlinearity_thresh}_{cell_line}.png")
            plt.close()

