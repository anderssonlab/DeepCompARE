import pandas as pd
import seaborn as sns
import matplotlib.pyplot as plt
from loguru import logger
import json

import os
import sys
sys.path.insert(1,"/isdata/alab/people/pcr980/DeepCompare/Scripts_python")
from stat_tests import match_by_decile
from scipy.stats import mannwhitneyu

def calculate_gc_bias(df,gc_col):
    gc_true = df[df['chip_evidence'] == True][gc_col]
    gc_false = df[df['chip_evidence'] == False][gc_col]
    _, pval = mannwhitneyu(gc_true, gc_false)
    if pval > 0.05:
        return "No_sig_diff"
    else:
        if gc_true.median() > gc_false.median():
            return "GC_rich"
        else:
            return "AT_rich"
    


def analyze_remap_gc_context(tf):
    df=pd.read_csv(f"Pd1_motif_loc_and_context_gc/{tf}.csv",header=None,skiprows=1)
    df.columns=['start', 'end', 'protein', 'score', 'strand', 'chromosome','chip_evidence', 'motif_seq', 'motif_gc', 'context_gc_2bp', 'context_gc_10bp', 'context_gc_50bp', 'context_gc_100bp', 'context_gc_300bp']
    # make sure df contain both true and false for chip_evidence column
    df['chip_evidence']=df['chip_evidence'].astype(bool)
    if len(df['chip_evidence'].unique())!=2:
        logger.info(f"Dataframe for {tf} does not contain both true and false for chip_evidence column")
        return
    df=match_by_decile(df,'score','chip_evidence')
    df_melted = df.melt(id_vars=['chip_evidence'], 
                        value_vars=['motif_gc', 'context_gc_2bp', 'context_gc_10bp', 'context_gc_50bp', 'context_gc_100bp', 'context_gc_300bp'],
                        var_name='gc_type', value_name='gc_content')

    # plot
    plt.figure(figsize=(12, 6))
    custom_palette = {True: "#76b041", False: "#408941"}
    sns.violinplot(x='gc_type', 
                   y='gc_content', 
                   hue='chip_evidence', 
                   data=df_melted, 
                   split=True, 
                   inner='quart',
                   palette=custom_palette)
    plt.title(f'{tf}: Distribution of GC Content')
    plt.ylabel('GC percentage')
    plt.xlabel('Context Length')
    plt.xticks(rotation=45)
    plt.tight_layout()
    plt.savefig(f'Plots_remap_gc/{tf}.png')
    plt.close()
    
    # Mann-Whitney U test
    # compare the distribution of GC content between groups with and without ChIP evidence for each context length
    # if 4/5 of the tests are significantly and consistently bias towards GC (test via median), then the TF is GC-rich
    # if 4/5 of the tests are significantly and consistently bias towards AT (test via median), then the TF is AT-rich
    # else, no significant difference
    test_res=[]
    for colname in ['context_gc_2bp', 'context_gc_10bp', 'context_gc_50bp', 'context_gc_100bp', 'context_gc_300bp']:
        test_res.append(calculate_gc_bias(df,colname))
    if test_res.count("GC_rich")>=4:
        return "GC_rich"
    elif test_res.count("AT_rich")>=4:
        return "AT_rich"
    else:
        return "No_sig_diff"
        

    

    
if __name__ == "__main__":
    
    # list all the files in /isdata/alab/people/pcr980/DeepCompare/Motif_gc_context/Pd1_motif_loc_and_context_gc
    files=os.listdir("/isdata/alab/people/pcr980/DeepCompare/Motif_gc_context/Pd1_motif_loc_and_context_gc")
    tfs=[f.split('.')[0] for f in files]
    
    res={"GC_rich":[],"AT_rich":[],"No_sig_diff":[]}
    for tf in tfs:
        tf_type=analyze_remap_gc_context(tf)
        if tf_type is not None:
            res[tf_type].append(tf)
        
    with open('gc_context_res.json', 'w') as f:
        json.dump(res, f)
        
        
        