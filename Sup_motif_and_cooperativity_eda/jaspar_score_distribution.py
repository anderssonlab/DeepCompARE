import pyranges as pr
import sys
import argparse
from loguru import logger

sys.path.insert(1,"/isdata/alab/people/pcr980/Scripts_python")
from seq_annotators import JasparAnnotator

#--------------------------------------------------------------------------------
# Step 1: get scores
#--------------------------------------------------------------------------------

def main(file_name):
    df_regions=pr.read_bed(f"/isdata/alab/people/pcr980/DeepCompare/Pd4_promoters_enhancers_and_featimp/{file_name}.bed").df.iloc[:,0:3]
    df_regions.columns=["chromosome","start","end"]
    df_regions["end"]=df_regions["end"]-1
    if "hepg2" in file_name:
        jaspar_annotator=JasparAnnotator("/isdata/alab/people/pcr980/Resource/JASPAR2022_tracks/JASPAR2022_hg38.bb",
                                    score_thresh=0,
                                    chip_file="/isdata/alab/people/pcr980/Resource/ReMap2022/ReMap2022_hg38_hepg2.bed",
                                    rna_file="/isdata/alab/people/pcr980/DeepCompare/RNA_expression/expressed_tf_list_hepg2.tsv")                                    
    if "k562" in file_name:
        jaspar_annotator=JasparAnnotator("/isdata/alab/people/pcr980/Resource/JASPAR2022_tracks/JASPAR2022_hg38.bb",
                                    score_thresh=0,
                                    chip_file="/isdata/alab/people/pcr980/Resource/ReMap2022/ReMap2022_hg38_k562.bed",
                                    rna_file="/isdata/alab/people/pcr980/DeepCompare/RNA_expression/expressed_tf_list_k562.tsv")
    jaspar_annotator.annotate(df_regions,outpath=f"{file_name}_tfbs.csv")


if __name__ == "__main__":
    parser = argparse.ArgumentParser(description='Annotate motif info')
    parser.add_argument('--file_name', type=str)
    
    args=parser.parse_args()
    main(args.file_name)
    logger.info(f"Finished annotating {args.file_name}.")




#--------------------------------------------------------------------------------
# Step 2: plot distribution
#--------------------------------------------------------------------------------
import pandas as pd
import seaborn as sns
import matplotlib.pyplot as plt

df_promoters_hepg2=pd.read_csv("promoters_hepg2_tfbs.csv")
df_promoters_hepg2["dataset"]="promoters_hepg2"
df_promoters_k562=pd.read_csv("promoters_k562_tfbs.csv")
df_promoters_k562["dataset"]="promoters_k562"
df_enhancers_hepg2=pd.read_csv("enhancers_hepg2_tfbs.csv")
df_enhancers_hepg2["dataset"]="enhancers_hepg2"
df_enhancers_k562=pd.read_csv("enhancers_k562_tfbs.csv")
df_enhancers_k562["dataset"]="enhancers_k562"


df=pd.concat([df_promoters_hepg2,df_promoters_k562,df_enhancers_hepg2,df_enhancers_k562])

# index are the datasets
df_count_all=pd.DataFrame()
for threshold in [0,100,200,300,400,500,600]:
    df["above_threshold"]=df["score"].apply(lambda x: x>threshold)
    df_count=df.groupby(["dataset","above_threshold"]).size().unstack().reset_index()
    df_count["threshold"]=threshold
    df_count_all=pd.concat([df_count_all,df_count])


# bar plot on the count:
sns.barplot(data=df_count_all,x="threshold",y=1,hue="dataset")
plt.savefig("jaspar_score_distribution.pdf")
plt.close()