import pandas as pd
import numpy as np
import pynvml




def find_available_gpu():
    """
    Find first available GPU
    Returns:
        str: GPU id
    """
    pynvml.nvmlInit()
    deviceCount = pynvml.nvmlDeviceGetCount()
    for i in range(deviceCount):
        handle = pynvml.nvmlDeviceGetHandleByIndex(i)
        mem = pynvml.nvmlDeviceGetMemoryInfo(handle)
        if mem.free/1024**3>2:
            return str(i)
    raise ValueError("No available GPU found!")


def remove_nan_inf(x,y):
    assert len(x)==len(y)
    mask_x_nan = np.isnan(x)
    mask_x_inf = np.isinf(x)
    mask_y_nan = np.isnan(y)
    mask_y_inf = np.isinf(y)
    mask_either = mask_x_nan | mask_x_inf | mask_y_nan | mask_y_inf
    return x[~mask_either],y[~mask_either]



    
def read_featimp(featimp_file,track_num):
    """
    Read featimp from featimp_file, subset by track_num
    Featimp is either gradxinp or ism, no header
    """
    featimp_df=pd.read_csv(featimp_file,header=None,index_col=0)
    # Given that indices are composed of "SeqX_TrackY", we can subset to contain only "_Track{track_num}"
    featimp_df=featimp_df[featimp_df.index.str.contains(f"_Track{track_num}$")]
    return featimp_df



def get_track_num(dataset,classification=False):
    if classification:
        if "hepg2" in dataset:
            return [0,2,4,6,8,10,12,14]
        elif "k562" in dataset:
            return [1,3,5,7,9,11,13,15]
        elif "common" in dataset:
            return [0,1,2,3,4,5,6,7,8,9,10,11,12,13,14,15]
        else:
            return None
    else:
        if "hepg2" in dataset:
            return [0,2,4,6]
        elif "k562" in dataset:
            return [1,3,5,7]
        elif "common" in dataset:
            return [0,1,2,3,4,5,6,7]
        else:
            return None



def change_track_name(df):
    track_dict={"track0":"cage","track1":"cage",
                "track2":"dhs","track3":"dhs",
                "track4":"starr","track5":"starr",
                "track6":"sure","track7":"sure",
                "_0":"cage","_1":"cage",
                "_2":"dhs","_3":"dhs",
                "_4":"starr","_5":"starr",
                "_6":"sure","_7":"sure"}
    for col in df.columns:
        col_suffix=col.split("_")[-1]
        if col_suffix in track_dict:
            df=df.rename(columns={col:col.replace(col_suffix,track_dict[col_suffix])})
    return df





def split_dimer(tf_list):
    res_list=[]
    for tf in tf_list:
        if "::" not in tf:
            res_list.append(tf)
        else:
            res_list.extend(tf.split("::"))
    return list(set(res_list))


def abs_max(vals):
    vals[max(enumerate(vals), key = lambda x: abs(x[1]))[0]]



def remove_intersection(list1,list2):
    intersection=set(list1).intersection(set(list2))
    list1_new=[item for item in list1 if item not in intersection]
    list2_new=[item for item in list2 if item not in intersection]
    return list1_new,list2_new

