import pandas as pd
from scipy.stats import pearsonr

def read_file(file_suffix):
    df=pd.read_csv(f'/isdata/alab/people/pcr980/DeepCompare/Pd5_motif_info/motif_info_thresh_500_{file_suffix}.csv',index_col=0)
    df["max_af"] = df["gnomad_af"].apply(lambda x: max([float(i) for i in str(x).split(":")]))
    df["max_241way"] = df["phylop_241way"].apply(lambda x: max([float(i) for i in str(x).split(":")]))
    for track_num in range(8):
        df[f"max_gradxinp_{track_num}"] = df[f"gradxinp_{track_num}"].apply(lambda x: max([float(i) for i in str(x).split(":")]))
        df[f"min_gradxinp_{track_num}"] = df[f"gradxinp_{track_num}"].apply(lambda x: min([float(i) for i in str(x).split(":")]))
        df[f"max_ism_{track_num}"] = df[f"ism_{track_num}"].apply(lambda x: max([float(i) for i in str(x).split(":")]))
        df[f"max_isa_{track_num}"] = df[f"isa_{track_num}"].apply(lambda x: max([float(i) for i in str(x).split(":")]))
    # remove columns starting with "ism" and "gradxinp" and "isa"
    cols_to_remove = [col for col in df.columns if col.startswith("ism") or col.startswith("gradxinp")]
    df.drop(cols_to_remove, axis=1, inplace=True)
    df["dataset"]=file_suffix
    return df


df=read_file("promoters_hepg2")

df.columns

# pearson correlation
pearsonr(df["max_gradxinp_0"],df["min_gradxinp_0"])
pearsonr(df["max_gradxinp_0"],df["isa_track0"])