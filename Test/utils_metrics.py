import numpy as np
import pandas as pd




def create_colnames():
    colnames=[]
    for modality in ["cage","dhs","starr","sure"]:
        for cell in ["hepg2","k562"]:
            colnames.append("_".join([modality,cell]))
    return colnames

def add_suffix(my_list,suffix):
    return [item+suffix for item in my_list]


def mask_neg1(pred,truth):
    # pred and truth are numpy arrays
    # mask out the -1 values in truth, and corresponding values in pred
    # return the masked pred and truth
    mask=(truth!=-1)
    pred=pred[mask]
    truth=truth[mask]
    return pred,truth


def calc_acc(pred,truth):
    pred=np.array(pred).squeeze()
    pred=(pred>0)
    truth=np.array(truth).squeeze()
    pred,truth=mask_neg1(pred,truth)
    return float((pred==truth).mean())
    
def calc_pcc(pred,truth):
    pred=np.array(pred).squeeze()
    truth=np.array(truth).squeeze()
    return float(np.corrcoef(pred,truth)[0][1])

def calc_columnwise_metrics(df_pred,df_truth,metric,metric_func):
    res=[]
    assert df_pred.shape==df_truth.shape
    for i in range(df_pred.shape[1]):
        res.append(metric_func(df_pred.iloc[:,i],df_truth.iloc[:,i]))
    df_res=pd.DataFrame({"file":df_pred.columns,metric:res})
    return df_res

