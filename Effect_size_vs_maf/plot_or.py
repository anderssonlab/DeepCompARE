import pandas as pd
import numpy as np
from loguru import logger

import sys
sys.path.insert(1,"/isdata/alab/people/pcr980/DeepCompare/Scripts_python")
from stat_tests import bin_two_columns, calc_or, plot_or, plot_or_jitter
from utils import get_track_num_list_from_file_name


#----------------------
# x axis: effect size
#----------------------
# for debug purpose
# file="enhancers_k562"
# track_num=1
for file in ["enhancers_k562","enhancers_hepg2","promoters_k562","promoters_hepg2"]:
    for track_num in get_track_num_list_from_file_name(file):
        logger.info(f"Processing {file}, track {track_num}")
        df=pd.read_csv(f"maf_with_effect_size_{file}.csv",header=None,index_col=0)
        df.reset_index(drop=True,inplace=True)
        df.columns=["chromosome","start","end","ID","REF","ALT","AF",'Name','Score','Strand']+ ["track_"+str(i) for i in range(16)]
        df["log10_AF"]=np.log10(df["AF"])
        df=bin_two_columns(df,
                f"track_{track_num}",
                [-np.inf,-0.2,-0.1,0,0.1,0.2,np.inf],
                "log10_AF",
                {"-inf - -3": "rare", "-3 - -2": "low", "-2 - 0":"common"},
                [-np.inf,-3,-2,0],
                "variant_type")
        df=df[df.sum(axis=1)>10]
        df_plot=calc_or(df,"Predicted_effect_size","variant_type",out_group="low")
        plot_or(df_plot, 'Predicted_effect_size', 'odds_ratio', "variant_type",
                f"Odds ratio ({file}, track {track_num})",
                {'rare': "#1f77b4", 'common': '#ff7f0e'},
                f"Plots/or_{file}_track{track_num}.pdf")
        logger.info(f"Done {file}, track {track_num}")






#---------------------------
# x axis: allele frequency
#---------------------------


# for file in ["enhancers_k562","enhancers_hepg2","promoters_k562","promoters_hepg2"]:
#     for track_num in get_track_num_list_from_file_name(file):
#         logger.info(f"Processing {file}, track {track_num}")
#         df=pd.read_csv(f"maf_with_effect_size_{file}.csv",header=None,index_col=0)
#         df.reset_index(drop=True,inplace=True)
#         df.columns=["chromosome","start","end","ID","REF","ALT","AF",'Name','Score','Strand']+ ["track_"+str(i) for i in range(16)]
#         df["log10_AF"]=np.log10(df["AF"])
#         df=bin_two_columns(df,
#             "log10_AF",
#             [-np.inf,-4,-2,0],
#             f"track_{track_num}",
#             {"-inf - -0.2":"large_negative","-0.2 - 0.2":"small","0.2 - inf":"large_positive"},
#             [-np.inf,-0.2,0.2,np.inf],
#             "variant_type")
#         df=df[df.sum(axis=1)>10]
#         df=df.loc[:, (df != 0).all(axis=0)]
#         df=calc_or(df)
#         df_plot = calc_or(df,'log10(Allele frequency)', 'variant_type')
#         plot_or_jitter(df_plot,'log10(Allele frequency)', 'odds_ratio',
#                        f"Odds ratio ({file}, track {track_num})",
#                        f"Plot_or_reverse/or_{file}_track{track_num}.pdf",
#                        {"large_negative": "#1f77b4","small": "#2ca02c","large_positive": "#9467bd"}
#                        )
#         logger.info(f"Done {file}, track {track_num}")









# nohup python3 plot_or.py > plot_or.out &
