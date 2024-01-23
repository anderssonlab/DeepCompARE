import pyBigWig
import pyranges as pr
import pandas as pd
from loguru import logger

def write_conservation_score(bed_path, bw_path, out_path):
    bw=pyBigWig.open(bw_path)
    gr=pr.read_bed(bed_path)
    for i in range(gr.df.shape[0]): 
        interval=gr.df.iloc[i].to_dict()
        chrom, start, end = interval['Chromosome'], interval['Start'], interval['End']
        values = bw.values(chrom, start, end)
        pd.DataFrame(values).T.to_csv(out_path, mode='a', header=False, index=False)


for conservation_score in ["phyloP","phastCons"]:
    for modality in ["CAGE","DHS","STARR","SuRE"]:
        for cell in ["HepG2","K562"]:
            logger.info(f"Processing {conservation_score} {modality} {cell}")
            bw_path=f"/isdata/alab/people/pcr980/Resource/hg38.{conservation_score}20way.bw"
            bed_path=f"/isdata/alab/people/pcr980/DeepCompare/Pd1_bed_processed/resize_600bp_{modality}_{cell}.bed"
            out_path=f"/isdata/alab/people/pcr980/DeepCompare/BioApp_evolutionary_conservation/{conservation_score}_{modality}_{cell}.csv"
            write_conservation_score(bed_path, bw_path, out_path)
            logger.info(f"Done {conservation_score} {modality} {cell}")