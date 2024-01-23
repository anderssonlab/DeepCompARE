import pandas as pd
import numpy as np
from scipy.stats import pearsonr

track_info={0: "CAGE_HepG2",
            1: "CAGE_K562",
            2: "DHS_HepG2",
            3: "DHS_K562",
            4: "STARR_HepG2",
            5: "STARR_K562",
            6: "SuRE_HepG2",
            7: "SuRE_K562"}

corrs=[]
pvals=[]
tracks=[]
for i in range(8):
    print(f"Track {i}: {track_info[i]}")
    df_ism=pd.read_csv(f"ism_{track_info[i]}.csv",index_col=0)
    df_gradxinp=pd.read_csv(f"gradxinp_{track_info[i]}.csv",header=None,index_col=0)
    print(df_ism.index)
    print(df_gradxinp.index)
    assert np.all(df_ism.index==df_gradxinp.index)
    ism=df_ism.values.reshape(-1,1).squeeze()
    gradxinp=df_gradxinp.values.reshape(-1,1).squeeze()
    assert len(ism)==len(gradxinp)
    mask_ism_nan = np.isnan(ism)
    mask_ism_inf = np.isinf(ism)
    mask_gradxinp_nan = np.isnan(gradxinp)
    mask_gradxinp_inf = np.isinf(gradxinp)
    mask_either = mask_ism_nan | mask_ism_inf | mask_gradxinp_nan | mask_gradxinp_inf
    corr,pval=pearsonr(ism[~mask_either],gradxinp[~mask_either])
    corrs.append(corr)
    pvals.append(pval)
    tracks.append(track_info[i])
    
df=pd.DataFrame({"track":tracks,"corr":corrs,"pval":pvals})
df.to_csv("correlation_ism_gradxinp.csv",index=False)
    
    