import pandas as pd
import numpy as np
import os
from loguru import logger


from prediction import compute_predictions
from utils import get_track_num
from seq_ops import ablate_motifs


#-------------------------------
# For all pair permutations
#-------------------------------


def get_permutations(n):
    temp=[(i,j) for i in range(n) for j in range(n)]
    temp=[(i,j) for i,j in temp if i!=j]
    return temp



def extract_motif(jaspar_annotator,region):
    """
    Given a genomic region, extract the motif information.
    Args:
        region: A tuple of (chromosome, start, end).
    Returns:
        A data frame of location, TF, Jaspar score and strand of ChIP-supported motif
    """
    df_motif=jaspar_annotator.annotate(region)
    if df_motif.shape[0]==0:
        return df_motif
    # df_motif.drop(columns=['rna_evidence'],inplace=True)
    df_motif["start_rel"]=df_motif["start"]-region[1]
    df_motif["end_rel"]=df_motif["end"]-region[1]
    return df_motif




def get_predictions(seq_col,device,prefix):
    pred=compute_predictions(seq_col,device=device)
    df_pred=pd.DataFrame(pred,columns=[f"pred_{prefix}_track{i}" for i in range(16)])
    return df_pred




def get_isas(df):
    """
    Given a dataframe of mutation results, calculate isa score.
    Args:
        df: A dataframe of mutation results.
    Returns:
        A dataframe of mutation results with isa score.
    """
    for i in range(16):
        df[f"isa1_track{i}"]=df[f"pred_orig_track{i}"]-df[f"pred_mut1_track{i}"]
        df[f"isa2_track{i}"]=df[f"pred_orig_track{i}"]-df[f"pred_mut2_track{i}"]
        df[f"isa_both_track{i}"]=df[f"pred_orig_track{i}"]-df[f"pred_mut_both_track{i}"]
        df[f"isa1_track{i}"]=df[f"isa1_track{i}"].round(3)
        df[f"isa2_track{i}"]=df[f"isa2_track{i}"].round(3)
        df[f"isa_both_track{i}"]=df[f"isa_both_track{i}"].round(3)
    return df




def write_pair_mutation_res_for_one_region(seq_extractor,jaspar_annotator,region,out_path,region_idx,device):

    df_motif=extract_motif(jaspar_annotator,region)
    if df_motif.shape[0]==0:
        return
    
    # get indices of all possible motif pairs
    combination_indices=get_permutations(df_motif.shape[0])
    df_mutation = pd.DataFrame(combination_indices, columns=['idx1', 'idx2'])
    df_mutation["region_idx"]=f"Region{region_idx}"
    df_mutation = df_mutation.merge(df_motif, left_on='idx1', right_index=True)
    df_mutation = df_mutation.merge(df_motif, left_on='idx2', right_index=True, suffixes=('', '2'))
    df_mutation.rename(columns={
        'chromosome':'chromosome1',
        'start':'start1',
        'end':'end1',
        'protein':'protein1',
        'score':'score1',
        'strand':'strand1',
        'chip_evidence':'chip_evidence1',
        'start_rel':'start_rel1',
        'end_rel':'end_rel1'},inplace=True) 

    df_mutation.drop(columns=['region2'],inplace=True)
    # move "region_idx" and "region" to the first 2 column
    df_mutation=df_mutation[["region_idx","region"]+[col for col in df_mutation.columns if col not in ["region_idx","region"]]]
    # remove overlapping motifs
    idx_olap_row=df_mutation.apply(lambda row: max(row['start_rel1'],row['start_rel2'])<=min(row['end_rel1'],row['end_rel2']),axis=1)
    df_mutation=df_mutation[~idx_olap_row].reset_index(drop=True)
    if df_mutation.shape[0]==0:
        return
    
    # get sequences
    seq = seq_extractor.get_seq(region)
    df_mutation["seq_orig"]=seq
    df_mutation['seq_mut1'] = df_mutation.apply(lambda row: ablate_motifs(seq, row['start_rel1'], row['end_rel1']), axis=1)
    df_mutation['seq_mut2'] = df_mutation.apply(lambda row: ablate_motifs(seq, row['start_rel2'], row['end_rel2']), axis=1)
    df_mutation['seq_mut_both'] = df_mutation.apply(lambda row: ablate_motifs(seq, [row['start_rel1'], row['start_rel2']],[row['end_rel1'], row['end_rel2']]), axis=1)
 
    # get predictions
    df_pred_orig=get_predictions(df_mutation["seq_orig"],device=device,prefix="orig").round(3)
    df_pred_mut1=get_predictions(df_mutation["seq_mut1"],device=device,prefix="mut1").round(3)
    df_pred_mut2=get_predictions(df_mutation["seq_mut2"],device=device,prefix="mut2").round(3)
    df_pred_mut_both=get_predictions(df_mutation["seq_mut_both"],device=device,prefix="mut_both").round(3)

    # concatenate predictions with df_mutation by columns
    df_mutation=pd.concat([df_mutation,df_pred_orig,df_pred_mut1,df_pred_mut2,df_pred_mut_both],axis=1)

    # positive isa score indicate the original motif is contributing positively to RE activity.
    df_mutation=get_isas(df_mutation)
    
    columns_to_drop=[col for col in df_mutation.columns if col.startswith("seq")]+\
                    [col for col in df_mutation.columns if col.startswith("pred_mut")]+\
                    ["idx1","idx2"]
    df_mutation.drop(columns=columns_to_drop,inplace=True)
    # if out_path does not exist, create it using mode="w", include header
    if not os.path.exists(out_path):
        df_mutation.to_csv(out_path,mode="w",header=True,index=False)
        return
    else:
        df_mutation.to_csv(out_path,mode="a",header=False,index=False)
        return 




def write_pair_mutation(df_regions,
                        seq_extractor,
                        jaspar_annotator,
                        device,
                        out_path):
    """
    First three columns of df_regions should be chromosome, start, end
    """
    for idx in range(df_regions.shape[0]):
        if idx%100==0:
            logger.info(f"{idx} regions processed.")
        region=df_regions.iloc[idx,[0,1,2]].tolist()
        write_pair_mutation_res_for_one_region(seq_extractor,
                                               jaspar_annotator,
                                               region,
                                               out_path,
                                               idx,
                                               device)













#--------------------------------------------
# Functions to aggregate df_mutate_pair
#--------------------------------------------


def _read_file(file_name,independent_threshold, cooperative_threshold,track_nums=None):
    """
    Read in file, remove confusing tracks, calculate cooperativity (c), and return the dataframe with cooperativity information
    Possible values of cooperativity: 
        * independent
        * redundancy
        * synergy
        * gray_zone: independent_threshold <= c <= cooperative_threshold. Won't be used for calling redundant or synergistic nor calculating cooperativity_index
        * None: invalid track gives negative isa score or isa2_wo_protein1
    """
    # read file
    df=pd.read_csv(file_name)
    # check if there are completely duplicated rows
    if df.duplicated().sum()>0:
        logger.info(f"# duplicated rows: {df.duplicated().sum()}")
        df.drop_duplicates(inplace=True)    
    # get track number
    if track_nums is None:
        track_nums=get_track_num(file_name)
    tracks_to_remove=[i for i in range(16) if i not in track_nums]
    # remove all columns with track{i} and i is in tracks_to_remove
    col_suffix_to_remove=[f"_track{i}" for i in tracks_to_remove]
    col_names_to_remove=[col for col in df.columns if any([col.endswith(suffix) for suffix in col_suffix_to_remove])]
    df.drop(columns=col_names_to_remove,inplace=True)
    # calculate cooperativity index
    for i in track_nums:
        df[f"isa2_wo_protein1_track{i}"]=df[f'isa_both_track{i}']-df[f'isa1_track{i}']
        # remove potentially mistaking tracks, retain other tracks
        df.loc[(df[f"isa1_track{i}"]<0) | (df[f"isa2_track{i}"]<0),[col for col in df.columns if col.endswith(f"track{i}")]]=None
        df.loc[(df[f"isa2_wo_protein1_track{i}"]<0),[col for col in df.columns if col.endswith(f"track{i}")]]=None
        df[f"i_track{i}"]=df[f"isa2_track{i}"]-df[f"isa2_wo_protein1_track{i}"]
        # if c_track{i} is None, set cooperativity_track{i} to None
        df[f"cooperativity_track{i}"]=None
        df.loc[(df[f"i_track{i}"].abs()<=independent_threshold),f"cooperativity_track{i}"]="independent"
        df.loc[df[f"i_track{i}"]< -cooperative_threshold,f"cooperativity_track{i}"]="redundancy"
        df.loc[df[f"i_track{i}"]> cooperative_threshold,f"cooperativity_track{i}"]="synergy"
        df.loc[(df[f"i_track{i}"].abs()>independent_threshold) & (df[f"i_track{i}"].abs()<=cooperative_threshold),f"cooperativity_track{i}"]="gray_zone"
    # drop columns starting with "pred" and isa, and motif_length
    col_names_to_remove=[col for col in df.columns if col.startswith("pred") or col.startswith("isa")]
    df.drop(columns=col_names_to_remove,inplace=True)
    return df




def _sum_cooperativity_by_row(df):
    cooperativity_cols=[col for col in df.columns if "cooperativity" in col]
    df["redundancy_count"]=0
    df["synergy_count"]=0
    df["independent_count"]=0
    df["num_gray_zone"]=0
    for col in cooperativity_cols:
        df["redundancy_count"]+= (df[col]=="redundancy")
        df["synergy_count"]+= (df[col]=="synergy")
        df["independent_count"]+= (df[col]=="independent")
        df["num_gray_zone"]+= (df[col]=="gray_zone")
    df["num_valid_profiles"]=df["redundancy_count"]+df["synergy_count"]+df["independent_count"]+df["num_gray_zone"]
    df["num_cooperative_profiles"]=df["redundancy_count"]+df["synergy_count"]
    return df




def _remove_confusing_pairs(df):
    row_idx=df[(df["redundancy_count"]>0) & (df["synergy_count"]>0)].index
    logger.info(f"# Confusing tf pairs: {len(row_idx)}")    
    df.drop(row_idx,inplace=True)
    return df




def _preprocess_nonlinear_pairs(df,track_nums):
    # if cooperativity_track{i}==gray_zone or independent or None, set c_track{i}=0, 
    # only c with cooperativity_track{i}==redundancy or synergy will be considered
    for i in track_nums:
        df.loc[(df[f"cooperativity_track{i}"]=="gray_zone") | 
               (df[f"cooperativity_track{i}"]=="independent") |
                (df[f"cooperativity_track{i}"].isnull()),
                f"i_track{i}"]=0
    # aggregate cooperativity index from different tracks
    i_cols=[col for col in df.columns if "i_track" in col]
    df["i"]=df[i_cols].sum(axis=1)/df["num_cooperative_profiles"]
    df["synergy"]=(df["i"]>0).astype(int)
    cooperativity_cols=[col for col in df.columns if "cooperativity" in col]
    df.drop(columns=i_cols+cooperativity_cols+["redundancy_count","synergy_count","num_gray_zone","num_valid_profiles"],inplace=True)
    return df





def read_cooperativity(file_name,independent_threshold=0.01, cooperative_threshold=0.1,track_nums=None):
    df=_read_file(file_name,independent_threshold, cooperative_threshold,track_nums)
    df=_sum_cooperativity_by_row(df)
    df=_remove_confusing_pairs(df)
    # select only rows with at least 1 valid profile
    df=df[df["num_valid_profiles"]>=1].reset_index(drop=True)
    return df 



def calculate_tf_pair_synergy_score(df,track_nums):
    """
    df should be the c of various tf pairs
    For each tf pair,
    summarize the c into redundancy and synergy,
    calculate cooperativity_index
    """
    # 1. calculate extra features
    df["distance"]=np.abs(df["start1"]-df["start2"])
    df["distance_iqr"]=df["distance"]
    df["count"]=1
    # remove distance > 255 (receptive field)
    df=df[df["distance"]<=255].reset_index(drop=True)
    # 2. split independent and non-independent
    df_independent=df[(df["independent_count"]==df["num_valid_profiles"])].reset_index(drop=True)
    df_coop=df[df["num_cooperative_profiles"]>0].reset_index(drop=True)
    df_coop=_preprocess_nonlinear_pairs(df_coop,track_nums)
    
    # 3. for independent pairs: group by tf pair, count, calculate distance and distance IQR
    df_independent=df_independent.groupby(["protein1","protein2"]).agg({"count":"sum",
                                                                "distance":"median",
                                                                "distance_iqr": lambda x: np.percentile(x, 75) - np.percentile(x, 25)
                                                                }).reset_index()    
    df_independent.rename(columns={"count":"independent_count",
                                "distance":"independent_distance",
                                "distance_iqr":"independent_distance_iqr"
                                },inplace=True)


    # 4. for cooperative pairs, group by tf pair, get i_redundancy, i_synergy, count_redundancy, count_synergy
    df_coop1=df_coop.groupby(["protein1","protein2","synergy"]).agg({"i":"sum","count":"sum"}).reset_index()
    df_coop1["tf_pair"]=df_coop1["protein1"]+"_"+df_coop1["protein2"]
    df_coop1=df_coop1.pivot(index="tf_pair",columns="synergy",values=["i","count"]).reset_index()
    df_coop1.columns = ['_'.join(map(str, col)).strip('_') for col in df_coop1.columns]
    df_coop1.columns = [col.replace('_0', '_redundancy').replace('_1', '_synergy') for col in df_coop1.columns]
    df_coop1["protein1"]=df_coop1["tf_pair"].apply(lambda x: x.split("_")[0])
    df_coop1["protein2"]=df_coop1["tf_pair"].apply(lambda x: x.split("_")[1])
    df_coop1.drop("tf_pair",axis=1,inplace=True)
    df_coop1.fillna(0,inplace=True)
    # for cooperative pairs, group by tf pair, get cooperative distance and distance IQR
    df_coop2=df_coop.groupby(["protein1","protein2"]).agg({"distance":"median",
                                                                        "count":"sum",
                                                                        "distance_iqr": lambda x: np.percentile(x, 75) - np.percentile(x, 25)}).reset_index()
    # rename columns
    df_coop2.rename(columns={"distance":"cooperative_distance",
                                    "distance_iqr":"cooperative_distance_iqr",
                                    "count":"cooperative_count"
                                },inplace=True)
    # merge df_coop1 and df_coop2
    df_coop=pd.merge(df_coop1,df_coop2,on=["protein1","protein2"],how="inner")
    # calculate cooperativity
    df_coop["i_sum"]=df_coop["i_redundancy"].abs()+df_coop["i_synergy"]
    df_coop["synergy_score"]=df_coop["i_synergy"].abs()/df_coop["i_sum"]
    df_coop["synergy_fraction"]=df_coop["count_synergy"]/(df_coop["count_redundancy"]+df_coop["count_synergy"])
    
    # 5. merge independent and cooperative by protein1 and protein2
    df_res=pd.merge(df_independent,df_coop,on=["protein1","protein2"],how="outer")
    # fill na in "independent_count" and "cooperative_count" with 0
    df_res["independent_count"].fillna(0,inplace=True)
    df_res["cooperative_count"].fillna(0,inplace=True)
    # calculate linearity_index
    df_res["independence_score"]=df_res["independent_count"]/(df_res["independent_count"]+df_res["cooperative_count"])
    df_res=df_res.round(2)
    return df_res



def assign_cooperativity(df,c_sum_thresh,linearity_thresh,redun_thresh,codep_thresh):
    """
    Also remove tf pairs with:
        * total_count <= 10
        * low i_sum & low linearity index, indicating low total counts.
        * Each cooperative TF will have a cooperativity index and linearity index
        # Some independent TFs originally have cooperativity index, set them to NaN
    TF pairs params: c_sum_thresh=1, linearity_thresh=0.9
    TF params: c_sum_thresh=5, linearity_thresh=0.95
    """
    df["total_count"]=df["independent_count"]+df["cooperative_count"]
    df=df[df["total_count"]>10].reset_index(drop=True)
    # assign cooperativity
    df["cooperativity"]=None
    # assign cooperative tf pairs:  set i_sum>c_sum_thresh, cooperativity_index > codep_thresh, or cooperativity_index < redun_thresh
    df.loc[(df["i_sum"]>=c_sum_thresh)&(df["synergy_score"]>codep_thresh),"cooperativity"]="Synergistic"
    df.loc[(df["i_sum"]>=c_sum_thresh)&(df["synergy_score"]<redun_thresh),"cooperativity"]="Redundant"
    df.loc[(df["i_sum"]>=c_sum_thresh)&(df["synergy_score"]>=redun_thresh)&(df["synergy_score"]<=codep_thresh),"cooperativity"]="Intermediate"
    # assign independent TF pairs:  cooperativity_index is nan, or linearity index > 0.9 and i_sum<=1
    df.loc[df["synergy_score"].isna(),"cooperativity"]="Independent"
    df.loc[(df["i_sum"]<c_sum_thresh) & (df["independence_score"]>linearity_thresh),"cooperativity"]="Independent"
    # set the cooperativity_index of independent TFs to NaN
    df.loc[df["cooperativity"]=="Independent","synergy_score"]=np.nan
    # remove cooperativity is nan, because they can't be determined due to too few cases
    df=df[~df["cooperativity"].isna()].reset_index(drop=True)
    df["cooperativity"]=pd.Categorical(df["cooperativity"],categories=["Independent","Redundant","Intermediate","Synergistic"],ordered=True)
    return df





def calculate_tf_synergy_score(df):
    """
    df should be tf pair cooperativity index (pre-filter)
    """
    # get linearity index
    df_independence_info=df.groupby("protein2").agg({"independent_count":"sum","cooperative_count":"sum"}).reset_index()
    df_independence_info["total_count"]=df_independence_info["cooperative_count"]+df_independence_info["independent_count"]
    df_independence_info=df_independence_info[df_independence_info["total_count"]>10].reset_index(drop=True)
    df_independence_info["independence_score"]=df_independence_info["independent_count"]/(df_independence_info["independent_count"]+df_independence_info["cooperative_count"])
    logger.info("Min total count: "+str(df_independence_info["total_count"].min()))
    # process cooperative pairs
    df_coop=df[~df["synergy_score"].isnull()].reset_index(drop=True)
    df_coop=df_coop.groupby("protein2").agg({"i_redundancy":"sum","i_synergy":"sum","count_redundancy":"sum","count_synergy":"sum"}).reset_index()
    df_coop["i_sum"]=df_coop["i_redundancy"].abs()+df_coop["i_synergy"]
    df_coop["count_sum"]=df_coop["count_redundancy"]+df_coop["count_synergy"]
    df_coop["synergy_score"]=df_coop["i_synergy"]/(df_coop["i_redundancy"].abs()+df_coop["i_synergy"])
    df_coop["synergy_fraction"]=df_coop["count_synergy"]/df_coop["count_sum"]
    df_tf=pd.merge(df_coop,df_independence_info,on="protein2",how="outer")
    # fill NaN in column i_sum with 0
    df_tf["i_sum"]=df_tf["i_sum"].fillna(0)
    return df_tf


