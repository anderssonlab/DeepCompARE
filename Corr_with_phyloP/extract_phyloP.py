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


for conservation_score in ["100way","241way","447way","447wayLRT"]:
    for file in ["enhancers_hepg2"]:
    #for file in ["promoters_hepg2","promoters_k562","enhancers_k562"]:
        logger.info(f"Processing {conservation_score} {file}")
        bw_path=f"/isdata/alab/people/pcr980/Resource/Conservation/hg38.phyloP{conservation_score}.bw"
        bed_path=f"/isdata/alab/people/pcr980/DeepCompare/Pd4_promoters_enhancers_and_featimp/{file}.bed"
        out_path=f"/isdata/alab/people/pcr980/DeepCompare/Corr_with_phyloP/{conservation_score}_{file}.csv"
        write_conservation_score(bed_path, bw_path, out_path)
        logger.info(f"Done {conservation_score} {file}")
        
        
        
# nohup python3 extract_phyloP.py > extract_phyloP.out &