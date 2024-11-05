#!/usr/bin/python3

import sys
sys.path.insert(1,"/isdata/alab/people/pcr980/DeepCompare/Train/Utils/")
sys.path.insert(1,"/isdata/alab/people/pcr980/Scripts_python/")
from multitasking_models import AstigCRConv5D, AstigCRConv6D
from trainer import *

import json
import datetime
import os
import argparse
import shutil
import numpy as np
import pandas as pd
import torch
import torch.nn.functional as F
import torch.nn as nn
from torch.utils.data import Dataset



class npyDataset3(Dataset):
    """
    A numpy dataset created from npy arrays stores on hard disk, used for training
    """
    def __init__(self,x_path,y_path,z_path,x_shape,y_shape,z_shape):
        # x: seq
        # y: 0/1
        # z: log1p(signal)
        self.X= np.memmap(x_path, dtype='float32', mode='r',shape=x_shape)
        self.Y= np.memmap(y_path, dtype='float32', mode='r',shape=y_shape)
        self.Z= np.memmap(z_path, dtype='float32', mode='r',shape=z_shape)
        self.len=len(self.Y)       
    
    def __getitem__(self, idx):
        return self.X[idx],self.Y[idx],self.Z[idx]

    def __len__(self):
        return self.len



class npyDataset2(Dataset):
    """
    A numpy dataset created from npy arrays stores on hard disk, used for training
    """
    def __init__(self,x_path,z_path,x_shape,z_shape):
        # x: seq
        # y: 0/1
        # z: log1p(signal)
        self.X= np.memmap(x_path, dtype='float32', mode='r',shape=x_shape)
        self.Z= np.memmap(z_path, dtype='float32', mode='r',shape=z_shape)
        self.len=len(self.Z)       
    
    def __getitem__(self, idx):
        return self.X[idx],self.Z[idx]

    def __len__(self):
        return self.len




def get_data_size_dict(this_dir):
    d={}
    with open(os.path.join(this_dir,"file_info.txt"),"r") as f:
        content = f.read().splitlines()
        for line in content:
            line=line.strip()
            words=line.split(" ")
            d[words[1][:-4]]=int(words[0])-1
    return d



def report_training_summary(record_train,record_val,model_dir,param_dict):
    """
    Write learning curves into .csv files
    Write parameters in to .json
    """
    # report regression loss
    pd.DataFrame(record_train.losses).to_csv(os.path.join(model_dir,"lc_train_loss.csv"))
    pd.DataFrame(record_val.losses).to_csv(os.path.join(model_dir,"lc_val_loss.csv"))
    # report classification accuracy
    pd.DataFrame(record_train.accs).to_csv(os.path.join(model_dir,"lc_train_acc.csv"))
    pd.DataFrame(record_val.accs).to_csv(os.path.join(model_dir,"lc_val_acc.csv"))

    param_dict["final_train_loss"]=float(record_train.losses[-1].mean())
    param_dict["final_val_loss"]=float(record_val.losses[-1].mean())
    param_dict["final_train_acc"]=float(record_train.accs[-1].mean())
    param_dict["final_val_acc"]=float(record_val.accs[-1].mean())
    
    with open(os.path.join(model_dir,'params.json'), 'w') as f:
        json.dump(param_dict, f)



def main(model_dir,param_dict):

    data_dir=param_dict["data_dir"]
    data_size_dict=get_data_size_dict(data_dir)
    data=param_dict["file_prefix"]  # for multitasking models, data="dat"
    train_size=data_size_dict[data+"_train"]
    val_size=data_size_dict[data+"_val"]
    device=torch.device(param_dict["gpu"])
    # create datasets
    try:
        train_dat=npyDataset3(os.path.join(data_dir,data+"_train_X.npy"),
                                os.path.join(data_dir,data+"_train_Y.npy"),
                                os.path.join(data_dir,data+"_train_Z.npy"),
                                torch.Size([train_size, 4, 600]),
                                torch.Size([train_size, param_dict["num_tasks"]]),
                                torch.Size([train_size, param_dict["num_tasks"]]))

        val_dat=npyDataset3(os.path.join(data_dir,data+"_val_X.npy"),
                            os.path.join(data_dir,data+"_val_Y.npy"),
                            os.path.join(data_dir,data+"_val_Z.npy"),
                            torch.Size([val_size, 4, 600]),
                            torch.Size([val_size, param_dict["num_tasks"]]),
                            torch.Size([val_size, param_dict["num_tasks"]]))
    except FileNotFoundError:
        print("Cannot find all X,Y,Z files. Try load X and Z only")
        train_dat=npyDataset2(os.path.join(data_dir,data+"_train_X.npy"),
                                os.path.join(data_dir,data+"_train_Z.npy"),
                                torch.Size([train_size, 4, 600]),
                                torch.Size([train_size, param_dict["num_tasks"]]))

        val_dat=npyDataset2(os.path.join(data_dir,data+"_val_X.npy"),
                            os.path.join(data_dir,data+"_val_Z.npy"),
                            torch.Size([val_size, 4, 600]),
                            torch.Size([val_size, param_dict["num_tasks"]]))


    # train model
    model=eval(param_dict["model"])()
    model.to(device)
    optimizer=torch.optim.Adam(model.parameters(),lr=param_dict["lr"],weight_decay=1e-5)
    trainer=eval(param_dict["trainer"])(device,model,train_dat,val_dat,param_dict,optimizer)
    record_train,record_val=trainer.train()

    # save model
    model.cpu()
    torch.save(model,os.path.join(model_dir,"model.h5"))
    
    #  report learning result
    report_training_summary(record_train,record_val,model_dir,param_dict)



if __name__=="__main__":
    # parse parameters
    param_dict={} 
    parser = argparse.ArgumentParser(description="Input parameters for training")
    parser.add_argument("-m", "--model", type=str, required=True, help="Model name")
    parser.add_argument("-d", "--data_directory", type=str, required=True, help="Data directory")
    parser.add_argument("-f", "--file_prefix", type=str, required=True, help="Prefix of file name")
    parser.add_argument("-g", "--gpu_index", type=str, required=True, help="GPU index")
    parser.add_argument("-t", "--trainer", type=str, required=True, help="Trainer")
    parser.add_argument("-r", "--num_tasks", type=int,required=True,help="Number of tasks")
    parser.add_argument("-l", "--learning_rate", type=float,required=True,help="Learning rate")
    parser.add_argument("-n", "--n_epochs", type=int,required=True,help="Number of epochs")
    parser.add_argument("-b", "--batch_size",type=int,required=True,help="Batch_size")
    parser.add_argument("-w", "--weight", type=str,help="Weights for 8 modality regression tasks")
    parser.add_argument("-o", "--output_suffix", type=str, help="Suffix to output directory")

    args=parser.parse_args()
    param_dict["model"]=args.model
    param_dict["data_dir"]=os.path.join("/isdata/alab/people/pcr980/DeepCompare/Datasets",args.data_directory)
    param_dict["file_prefix"]=args.file_prefix
    param_dict["trainer"]=args.trainer
    param_dict["num_tasks"]=args.num_tasks
    param_dict["gpu"]="cuda:"+args.gpu_index
    param_dict["lr"]=args.learning_rate
    param_dict["n_epochs"]=args.n_epochs
    param_dict["batch_size"]=args.batch_size
    if args.weight:
        param_dict["weight"]=eval(args.weight)
    
    # create model directory
    model_dir_suffix=""
    if args.output_suffix:
        model_dir_suffix=args.output_suffix
    model_dir="_".join([param_dict["model"],args.data_directory,args.file_prefix,args.trainer])+model_dir_suffix
    model_dir=os.path.join("/isdata/alab/people/pcr980/DeepCompare/Models",model_dir)
    os.makedirs(model_dir)

    # run
    now = datetime.datetime.now()
    print("{} ".format(now)+"Create directory {}\n".format(model_dir))
    try:
        main(model_dir,param_dict)
        now = datetime.datetime.now()
        print("{} ".format(now)+"Succeed with "+" ".join(sys.argv))
    except:
        print("Error!\nRemove directory "+model_dir)
        shutil.rmtree(model_dir)
        raise
