import sys
sys.path.insert(1,"/isdata/alab/people/pcr980/DeepCompare/Scripts_python")
from utils import SeqExtractor
import pandas as pd
from loguru import logger

def bed2sequences(bed_file,out_path):
    """
    Args:
        bed_file: path to bed file
    Returns:
        a list of sequences
    """
    extractor = SeqExtractor()
    bed_df=pd.read_csv(bed_file,header=None,sep="\t")
    seqs = [extractor.get_seq(*interval) for interval in zip(bed_df.iloc[:,0],
                                                            bed_df.iloc[:,1],
                                                            bed_df.iloc[:,2]-1)]
    assert len(seqs)==len(bed_df)
    assert len(seqs[0])==600
    pd.DataFrame({"sequence":seqs}).to_csv(out_path,index=False)
    
if __name__=="__main__":
    for modality in ["CAGE","DHS","STARR","SuRE"]:
        for cell_type in ["HepG2","K562"]:
            bed_file=f"/isdata/alab/people/pcr980/DeepCompare/Pd1_bed_processed/resize_600bp_{modality}_{cell_type}.bed"
            out_path=f"/isdata/alab/people/pcr980/DeepCompare/Pd2_sequences/seqs_{modality}_{cell_type}.csv"
            bed2sequences(bed_file,out_path)
            logger.info(f"Done {modality} {cell_type}")