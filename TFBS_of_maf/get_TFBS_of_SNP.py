import sys
import pyranges as pr
import pandas as pd
import os
from loguru import logger

sys.path.insert(1,"/isdata/alab/people/pcr980/DeepCompare/Scripts_python")
from seq_annotators import JasparAnnotator, ReMapAnnotator, gnomADSNPAnnotator
from region_ops import subset_df_by_region




def find_overlapping_motif(motifs, snp):
    overlapping_motifs = subset_df_by_region(motifs,snp,by='1bp')
    if overlapping_motifs.shape[0]==0:
        return None, None, None
    overlapping_motifs_sorted = overlapping_motifs.sort_values(
        by=['chip_evidence', 'score'], ascending=[False, False])
    best_motif = overlapping_motifs_sorted.iloc[0,3]
    best_score= overlapping_motifs_sorted.iloc[0,4]
    chip_evidence = overlapping_motifs_sorted.iloc[0,6]
    return best_motif,best_score,chip_evidence



def annotate_one_region(region,gnomad_annotator,jaspar_annotator,remap_annotator):
    snps=gnomad_annotator.annotate(region)
    motifs=jaspar_annotator.annotate(region)
    motifs=remap_annotator.annotate(motifs,region)
    snps[['motif', 'score',"chip_evidence"]] = snps.apply(lambda snp: find_overlapping_motif(motifs, snp), axis=1).apply(pd.Series)
    return snps

    

def annotate_one_bed_file(file_name,gnomad_annotator):
    df_regions=pr.read_bed(f"/isdata/alab/people/pcr980/DeepCompare/Pd4_promoters_enhancers_and_featimp/{file_name}.bed").df
    df_regions.iloc[:,2]=df_regions.iloc[:,2]-1
    
    # load tools
    jaspar_annotator=JasparAnnotator()
    if "k562" in file_name:
        remap_annotator = ReMapAnnotator("/isdata/alab/people/pcr980/Resource/ReMap2022/ReMap2022_K-562.bed")
    if "hepg2" in file_name:
        remap_annotator = ReMapAnnotator("/isdata/alab/people/pcr980/Resource/ReMap2022/ReMap2022_Hep-G2.bed")


    for i in range(df_regions.shape[0]):
        if i%1000==0:
            print(f"Processing region {i}")
        region=df_regions.iloc[i,0:3].tolist()
        snps=annotate_one_region(region,gnomad_annotator,jaspar_annotator,remap_annotator)
        if os.path.exists(f"tfbs_maf_{file_name}.csv"):
            snps.to_csv(f"tfbs_maf_{file_name}.csv",mode="a",index=False,header=True)
        else:
            snps.to_csv(f"tfbs_maf_{file_name}.csv",mode="w",index=False)
        





if __name__=="__main__":
    gnomad_annotator=gnomADSNPAnnotator()
    
    logger.info("Annotating promoters hepg2")
    annotate_one_bed_file("promoters_hepg2",gnomad_annotator)
    logger.info("Annotating promoters k562")
    annotate_one_bed_file("promoters_k562",gnomad_annotator)
    logger.info("Annotating enhancers hepg2")
    annotate_one_bed_file("enhancers_hepg2",gnomad_annotator)
    logger.info("Annotating enhancers k562")
    annotate_one_bed_file("enhancers_k562",gnomad_annotator)
    
    
# nohup python3 get_TFBS_of_SNP.py > get_TFBS_of_SNP.out &