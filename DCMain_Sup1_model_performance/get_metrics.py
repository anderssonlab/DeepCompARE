# Bench mark models on a common test set


import os
import sys
import pandas as pd
from utils_metrics import calc_acc,calc_pcc,calc_columnwise_metrics,create_colnames,add_suffix

sys.path.insert(1,"/isdata/alab/people/pcr980/DeepCompare/Scripts_python/")
from prediction import write_predictions

dir_models="/isdata/alab/people/pcr980/DeepCompare/Models"
dir_predictions="Pd1_model_predictions/"


#------------------------------------------------------
# Task1: generate Pd1_model_predictions/, write predictions on unified standard test set
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
            print(f"Remove empty directory {model_name}")
            os.rmdir(pred_path)
        
        
    
#------------------------------------------------------
# Task2: Generate Pd2_metrics/ST_metrics.csv
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

def append_list_in_dict(dict_database,dict_new):
    for key,value in dict_new.items():
        dict_database[key].append(value)
        

# analysis
df_truth=pd.read_csv("/isdata/alab/people/pcr980/DeepCompare/Datasets/Dataset_final_rep/dat_test.csv")

res_dict_ST={
    "model_name":[],
    "task_type":[],
    "model_type":[],
    "file":[],
    "data_dir":[],
    "pcc":[],
    "acc":[]
}
for model_name in os.listdir(dir_predictions):
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
        # append metrics
        
        res_dict_ST["pcc"].append(pcc)
        res_dict_ST["acc"].append(acc)
        append_list_in_dict(res_dict_ST,info_dict)            

res_df=pd.DataFrame(res_dict_ST)
res_df.to_csv("Pd2_benchmark_metrics/ST_metrics.csv")

            


#------------------------------------------------------
# Task3: Generate Pd2_metrics/MT_metrics.csv
#------------------------------------------------------
df_truth=pd.read_csv("/isdata/alab/people/pcr980/DeepCompare/Datasets/Dataset_final_rep/dat_test.csv")

df_truth_reg=df_truth.loc[:,add_suffix(create_colnames(),"_intensity")]
df_truth_class=df_truth.loc[:,add_suffix(create_colnames(),"_class")]

# write header
df_metric=pd.DataFrame(columns=['file', 'pcc', 'acc', 'model_name', 'task_type', 'model_type', 'data_dir'])
df_metric.to_csv("Pd2_benchmark_metrics/MT_metrics.csv",index=False,header=True)
for model_name in os.listdir(dir_predictions):
    info_dict=split(model_name)
    if info_dict["task_type"]=="MT":
        print(model_name)
        df_pred=pd.read_csv(os.path.join(dir_predictions,model_name,"pred_dat_test.csv"),header=None)
        if info_dict["model_type"]=="CR":
            df_pred_reg=df_pred.iloc[:,0:8]
            df_pred_class=df_pred.iloc[:,8:16]
            df_pred_reg.columns=create_colnames()
            df_pred_class.columns=create_colnames()
            df_metric=calc_columnwise_metrics(df_pred_reg,df_truth_reg,'pcc',calc_pcc)
            df_acc=calc_columnwise_metrics(df_pred_class,df_truth_class,'acc',calc_acc)
            df_metric.loc[:,"acc"]=df_acc.acc
        elif info_dict["model_type"]=="Reg":
            df_pred.columns=create_colnames()
            df_metric=calc_columnwise_metrics(df_pred,df_truth_reg,'pcc',calc_pcc)
            df_metric.loc[:,"acc"]=0
        elif info_dict["model_type"]=="Class":
            df_pred.columns=create_colnames()
            df_metric=calc_columnwise_metrics(df_pred,df_truth_class,'acc',calc_acc)
            df_metric.loc[:,"pcc"]=0
        else: 
            print("Model type not recognized")
            
        df_metric.loc[:,"model_name"]=info_dict["model_name"]
        df_metric.loc[:,"task_type"]=info_dict["task_type"]
        df_metric.loc[:,"model_type"]=info_dict["model_type"]
        df_metric.loc[:,"data_dir"]=info_dict["data_dir"]   
        df_metric=df_metric.loc[:,['file', 'pcc', 'acc', 'model_name', 'task_type', 'model_type', 'data_dir']]
        df_metric.to_csv("Pd2_benchmark_metrics/MT_metrics.csv",mode='a',index=False,header=False)

