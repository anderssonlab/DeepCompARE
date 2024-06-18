import pandas as pd
import numpy as np
import pynvml
import re




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
        if mem.free/1024**3>20:
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



def extract_numbers(s):
    """
    Solely for sorting row names
    """
    return list(map(int, re.findall(r'\d+', s)))
    
    
def read_featimp(featimp_file,track_num):
    """
    Read featimp from featimp_file, subset by track_num
    Featimp is either gradxinp or ism, no header
    """
    featimp_df=pd.read_csv(featimp_file,header=None,index_col=0)
    # Given that indices are composed of "SeqX_TrackY", we can subset to contain only "_Track{track_num}"
    featimp_df=featimp_df[featimp_df.index.str.contains(f"_Track{track_num}$")]
    return featimp_df



def get_track_num_list_from_file_name(file):
    if "k562" in file:
        return list(range(1,8,2))
    elif "hepg2" in file:
        return list(range(0,8,2))
    return None