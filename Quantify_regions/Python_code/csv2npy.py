#!/usr/bin/python3
import numpy as np
import pandas as pd
import torch
import os

import sys
import getopt
sys.path.insert(1,"/isdata/alab/people/pcr980/DeepCompare/Scripts_python")
from utils import encode




def csv2npy(data_dir,fname,model_type,seq_colname,y_collist,z_collist):
    df=pd.read_csv(os.path.join(data_dir,fname),index_col=0)
    seqs=list(df.loc[:,seq_colname])

    X=np.zeros((len(seqs),len(seqs[0]),4))
    X=np.array(list(map(encode,seqs)))
    X=torch.from_numpy(X)
    X=torch.permute(X, (0, 2, 1))
    f_X = np.memmap(os.path.join(data_dir,fname[:-4]+"_X.npy"), dtype='float32', mode='w+', shape=X.shape)
    f_X[:]=X[:]
    f_X.flush()
    
    if model_type in ["classification","CR","Bin"]:
        print("classification column selected:")
        print(",".join(y_collist))
        Y=np.array(df.loc[:,y_collist])
        f_Y = np.memmap(os.path.join(data_dir,fname[:-4]+"_Y.npy"), dtype='float32', mode='w+', shape=Y.shape)
        f_Y[:]=Y[:]
        f_Y.flush()

    if model_type in ["regression","CR","Bin"]:
        print("regression column selected:")
        print(",".join(z_collist))
        Z=np.array(df.loc[:,z_collist])
        f_Z = np.memmap(os.path.join(data_dir,fname[:-4]+"_Z.npy"), dtype='float32', mode='w+', shape=Z.shape)
        f_Z[:]=Z[:]
        f_Z.flush()


def add_suffix(alist,suffix):
    return [element+suffix for element in alist]

if __name__=="__main__":
    argv=sys.argv[1:]
    opts, args=getopt.getopt(argv,"d:m:f:j:")
    for opt,arg in opts:
        if opt=="-d":
            data_dir=arg
        if opt=="-m":
            model_type=arg
        if opt=="-f":
            if arg=="all":  # process all files?
                files=os.listdir(data_dir)
                files=[f for f in files if f.endswith(".csv")]
            else:
                files=[arg]
        if opt=="-j":      # process all columns? 
            if arg=="yes":
                column="all" 
            if arg=="no":
                column="infer"
    
    key_list=["cage_hepg2","cage_k562","dhs_hepg2","dhs_k562","starr_hepg2","starr_k562","sure_hepg2","sure_k562"]

    for f in files:
        print("Process file "+f)
        if column=="infer":
            name_elements=f.split("_")
            y_collist=["_".join([name_elements[0],name_elements[1],"class"])]
            z_collist=["_".join([name_elements[0],name_elements[1],"intensity"])]
        
        if column=="all":
            y_collist=add_suffix(key_list,"_class")
            if model_type in ["CR","regression"]:
                z_collist=add_suffix(key_list,"_intensity")
        
        csv2npy(data_dir,f,model_type,"seq",y_collist,z_collist)

print("End!")


