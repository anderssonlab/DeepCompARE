from write_ism_A_seq import write_ism_score,create_dir
import torch
import argparse
import pandas as pd
import os
from seq_ops import SeqExtractor
    
    
    
def write_ism_from_bed(bed_file,out_path,device):
    # get sequences
    extractor = SeqExtractor()
    bed_df=pd.read_csv(bed_file,header=None,sep="\t")
    intervals = [(chrom,start,end) for chrom,start,end in zip(bed_df.iloc[:,0],bed_df.iloc[:,1],bed_df.iloc[:,2]-1)] # -1 to not include the end
    seqs = [extractor.get_seq(*interval) for interval in intervals]
    
    # write sequences to file
    dir_temp=create_dir()
    pd.DataFrame(seqs,columns=["sequence"]).to_csv(os.path.join(dir_temp,"seqs.csv"),index=False)
    
    # call function to write ism score
    write_ism_score(os.path.join(dir_temp,"seqs.csv"),"sequence",out_path,dir_temp,device)




if __name__=="__main__":
    argparser = argparse.ArgumentParser()
    argparser.add_argument("--bed_file",type=str,help="path to sequence file")
    argparser.add_argument("--out_path",type=str,help="path to output file")
    argparser.add_argument("--device",type=str,help="gpu index")

    args=argparser.parse_args()
    write_ism_from_bed(args.bed_file,args.out_path,torch.device(f"cuda:{args.device}"))
    