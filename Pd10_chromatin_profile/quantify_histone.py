import pandas as pd
import numpy as np
import pyBigWig
import os

from loguru import logger


# TODO: add true CAGE/DHS/STARR/CAGE as well

def calculate_average_bw(bw, region):
    chrom, start, end = region
    values = bw.values(chrom, start, end, numpy=True)
    # calculate averge signal on non-NaN values
    values = values[~np.isnan(values)]
    return np.log1p(np.mean(values) if len(values) > 0 else 0)


def annotate(df,cell_line,out_path):
    bw_path="/isdata/alab/people/pcr980/DeepCompare/Pd10_chromatin_profile/Bigwig_data"
    bw_files=os.listdir(bw_path)
    # subset for files ending with f"_{cell_line}.bigwig"
    bw_files=[bw_file for bw_file in bw_files if bw_file.endswith(f"_{cell_line}.bigwig")]
    # remove f"_{cell_line}.bigwig"
    chromatin_signatures = [bw_file.split("_")[0] for bw_file in bw_files]
    for chromatin_signature in chromatin_signatures:
        bw = pyBigWig.open(f"Bigwig_data/{chromatin_signature}_{cell_line}.bigwig")
        df[f'log1p_{chromatin_signature}'] = df.apply(lambda row: calculate_average_bw(bw, (row[0], row[1], row[2])), axis=1)
        logger.info(f'Done with {chromatin_signature}')
    df.to_csv(out_path,index=False)


for re in ["proximal", "distal"]:
    for cell_line in ["hepg2","k562"]:
        df=pd.read_csv(f'/isdata/alab/people/pcr980/DeepCompare/Pd4_promoters_enhancers_and_featimp/dhs_{re}_{cell_line}.tsv', sep=' ')
        # retain first 3 columns
        df = df.iloc[:, :3]
        df.columns = ['chrom', 'start', 'end']
        annotate(df,cell_line,f'{re}_{cell_line}.csv')
        
        
        
for re in ["promoters", "enhancers"]:
    for cell_line in ["hepg2","k562"]:
        df=pd.read_csv(f'/isdata/alab/people/pcr980/DeepCompare/Pd4_promoters_enhancers_and_featimp/{re}_{cell_line}.bed', sep='\t', header=None)
        # retain first 3 columns
        df = df.iloc[:, :3]
        df.columns = ['chrom', 'start', 'end']
        annotate(df,cell_line,f'{re}_{cell_line}.csv')
        
        
        


# nohup python3 quantify_histone.py > quantify_histone.log &