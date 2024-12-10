import pandas as pd
import argparse
from loguru import logger
import pandas as pd
from pyfaidx import Fasta


from tf_cooperativity import write_pair_mutation
from seq_annotators import JasparAnnotator






class SeqExtractor:
    """
    Convert bed regions to fasta sequences
    """
    def __init__(self,fasta_path):
        self.fasta = Fasta(fasta_path)
    #
    def get_seq(self,region):
        if len(region) == 3:
            # deal with format like (chr,start,end) or [chr,start,end]
            chr, start, end = region
        elif isinstance(region, str):
            # deal with format like chr:start-end
            chr, start_end = region.split(":")
            start, end = start_end.split("-")
            start, end = int(start), int(end)
        else:
            raise ValueError("Invalid input.")
        return str(self.fasta.get_seq(chr,start,end).seq).upper()





#-----------------------
# Functions
#-----------------------

def analysis(file_name,device,sep="\t"):
    # load data and tools
    seq_extractor = SeqExtractor("/isdata/alab/people/pcr980/Resource/hg38.fa")
    jaspar_annotator=JasparAnnotator("/isdata/alab/people/pcr980/Resource/JASPAR2022_tracks/JASPAR2022_hg38.bb",
                                    score_thresh=500,
                                    rna_file="/isdata/alab/people/pcr980/DeepCompare/RNA_expression/expressed_tf_list_k562.tsv",
                                    subset="rna")
    df_regions = pd.read_csv(f"{file_name}.bed",sep=sep,header=None)
    df_regions.iloc[:,2]=df_regions.iloc[:,2]-1
    out_path = f"mutate_pairs_{file_name}_subset_by_jaspar.csv"
    write_pair_mutation(df_regions,seq_extractor,jaspar_annotator,device,out_path)


if __name__ == "__main__":
    parser=argparse.ArgumentParser()
    parser.add_argument("--file_name",type=str)
    parser.add_argument("--device",type=str)
    args=parser.parse_args()
    analysis(args.file_name,args.device)
    logger.info(f"Done with {args.file_name}")