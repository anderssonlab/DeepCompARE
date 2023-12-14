# Bench mark models on a common test


import os
import sys
import pandas as pd
import numpy as np

sys.path.insert(1,"/isdata/alab/people/pcr980/DeepCompare/Scripts_python")
from write_prediction import write_predictions

dir_models="/isdata/alab/people/pcr980/DeepCompare/Models"
dir_predictions="/isdata/alab/people/pcr980/DeepCompare/Test/"
#------------------------------------------------------
# Task1: write predictions on unified standard test set
#------------------------------------------------------
# write predictions if the model is well trained
for model_name in os.listdir(dir_models):
    model_path = os.path.join(dir_models, model_name)
    if os.path.isdir(model_path) and os.listdir(model_path): # if model is trained
        out_dir=os.path.join(dir_predictions,model_name)     # directory to output predictions
        if not os.path.exists(out_dir):
            os.makedirs(out_dir)
            print(f"Process {model_name}.")
            write_predictions(gpu_idx="7",
                            model_dir=model_path,
                            data_path="/isdata/alab/people/pcr980/DeepCompare/Datasets/Dataset_final_rep/dat_test.csv",
                            seq_colname="seq",
                            out_path=os.path.join(out_dir,"pred_dat_test.csv"))
        else:
            print(f"{model_name} already processed. Skip." )

            
# remove empty directories in Test/ if prediction fails in some way
for model_name in os.listdir(dir_predictions):
    pred_path=os.path.join(dir_predictions,model_name)
    if os.path.isdir(pred_path):
        if len(os.listdir(pred_path)) == 0:
            print(f"Remove {model_name}")
            os.rmdir(pred_path)
        
        
    
#------------------------------------------------------
# Task2: calculate either accuracy and classification
#------------------------------------------------------
# functions
def split(name):
    info_dict={}
    parts = name.split('_')
    info_dict["model_name"]=parts[0]
    info_dict["task_type"]=parts[-1]
    info_dict["model_type"]=parts[-2]
    if parts[-3]=="dat":
        info_dict["file"]=parts[-3]
        info_dict["data_dir"]="_".join(parts[1:-3])
    else:
        info_dict["file"]="_".join(parts[-4:-2])
        info_dict["data_dir"]="_".join(parts[1:-4])
    return info_dict

def calc_acc(pred,truth):
    pred=np.array(pred).squeeze()
    pred=(pred>0)
    truth=np.array(truth).squeeze()
    mask=(truth!=-1)
    pred=pred[mask]
    truth=truth[mask]
    return float((pred==truth).mean())
    
def calc_pcc(pred,truth):
    pred=np.array(pred).squeeze()
    truth=np.array(truth).squeeze()
    return float(np.corrcoef(pred,truth)[0][1])

def append_list_in_dict(dict_database,dict_new):
    for key,value in dict_new.items():
        dict_database[key].append(value)
        








# analysis
df_truth=pd.read_csv("/isdata/alab/people/pcr980/DeepCompare/Datasets/Dataset_final_rep/dat_test.csv")

        
res_dict={
    "model_name":[],
    "task_type":[],
    "model_type":[],
    "file":[],
    "data_dir":[],
    "pcc":[],
    "acc":[]
}
    

for model_name in os.listdir(dir_models):
    print(model_name)
    info_dict=split(model_name)
    if info_dict["task_type"]=="ST":
        pred=pd.read_csv(os.path.join(dir_predictions,model_name,"pred_dat_test.csv"),header=None)
        if info_dict["model_type"]=="CR":
            # the order is regression, classification
            pcc=calc_pcc(pred.iloc[:,0],df_truth.loc[:,info_dict["file"]+"_intensity"])
            acc=calc_acc(pred.iloc[:,1],df_truth.loc[:,info_dict["file"]+"_class"])
        elif info_dict["model_type"]=="Reg":
            pcc=calc_pcc(pred.iloc[:,0],df_truth.loc[:,info_dict["file"]+"_intensity"])
            acc=0
        elif info_dict["model_type"]=="Class":
            pcc=0
            acc=calc_acc(pred.iloc[:,0],df_truth.loc[:,info_dict["file"]+"_class"])
        else: 
            print("Model type not recognized")
        
        res_dict["pcc"].append(pcc)
        res_dict["acc"].append(acc)
        append_list_in_dict(res_dict,info_dict)            

for key,value in res_dict.items():
    print(key)
    print(len(value))
        
res_df=pd.DataFrame(res_dict)

res_df.to_csv(os.path.join(dir_predictions,"ST_metrics.csv"))


