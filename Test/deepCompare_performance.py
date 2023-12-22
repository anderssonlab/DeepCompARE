#!/usr/bin/python3

import pandas as pd
from sklearn.metrics import roc_auc_score,roc_curve,auc,precision_recall_curve
import matplotlib.pyplot as plt
from matplotlib.pyplot import figure


import sys
sys.path.insert(1,"/isdata/alab/people/pcr980/DeepCompare/Scripts_python")
from prediction import write_predictions 
from utils_metrics import *

#-------------------------------------------------------
# Task1: Write predictions
#------------------------------------------------------
# for suffix in ["train", "test", "val"]:
#     write_predictions(data_path=f"/isdata/alab/people/pcr980/DeepCompare/Datasets/Dataset_final_rep/dat_{suffix}.csv",
#                     seq_colname="seq",
#                     out_path=f"/isdata/alab/people/pcr980/DeepCompare/Test/Pd3_DeepCompare_performance/pred_dat_{suffix}.csv")


#-------------------------------------------------------
# Task2: Write metadata
#------------------------------------------------------

# # create empty df
# colnames=["seqnames","start","end","width","strand","seq",
#           "cage_hepg2_class","cage_k562_class","dhs_hepg2_class","dhs_k562_class",
#           "starr_hepg2_class","starr_k562_class","sure_hepg2_class","sure_k562_class",
#           "cage_hepg2_intensity","cage_k562_intensity","dhs_hepg2_intensity","dhs_k562_intensity",
#           "starr_hepg2_intensity","starr_k562_intensity","sure_hepg2_intensity","sure_k562_intensity",
#           "data_type",
#           'signal_cage_hepg2','signal_cage_k562','signal_dhs_hepg2','signal_dhs_k562',
#           'signal_starr_hepg2','signal_starr_k562', 'signal_sure_hepg2','signal_sure_k562',
#           'z_cage_hepg2','z_cage_k562','z_dhs_hepg2','z_dhs_k562',
#           'z_starr_hepg2','z_starr_k562', 'z_sure_hepg2','z_sure_k562']

# df_res=pd.DataFrame(columns=colnames)
# df_res.to_csv("/isdata/alab/people/pcr980/DeepCompare/Metadata_and_prediction/metadata_all_seqs.csv", mode='a')

# id_start=0
# for data_type in ["train", "val", "test"]:
#     # read in files
#     df_orig=pd.read_csv(f"/isdata/alab/people/pcr980/DeepCompare/Datasets/Dataset_final_rep/dat_{data_type}.csv")
#     df_pred=pd.read_csv(f"/isdata/alab/people/pcr980/DeepCompare/Test/Pd3_DeepCompare_performance/pred_dat_{data_type}.csv",header=None)
#     df_gr=pd.read_csv(f"/isdata/alab/people/pcr980/DeepCompare/Datasets/Dataset_final_rep/gr_{data_type}.bed",sep='\t',header=None)
#     assert df_pred.shape[0]==df_orig.shape[0]
#     assert df_pred.shape[0]==df_gr.shape[0]
    
#     df_gr.columns=["seqnames","start","end","unknown","width","strand"]
#     df_pred.columns=['signal_cage_hepg2','signal_cage_k562','signal_dhs_hepg2','signal_dhs_k562',
#                     'signal_starr_hepg2','signal_starr_k562', 'signal_sure_hepg2','signal_sure_k562',
#                     'z_cage_hepg2','z_cage_k562','z_dhs_hepg2','z_dhs_k562',
#                     'z_starr_hepg2','z_starr_k562', 'z_sure_hepg2','z_sure_k562']
#     # create unique identifier and set as index
#     ids=list(range(id_start,id_start+df_orig.shape[0]))
#     identifiers=["Seq" + str(i) for i in ids]
#     id_start+=df_orig.shape[0]
#     df_orig["identifier"]=identifiers
#     df_pred["identifier"]=identifiers
#     df_gr["identifier"]=identifiers
#     df_orig.set_index("identifier",inplace=True)
#     df_pred.set_index("identifier",inplace=True)
#     df_gr.set_index("identifier",inplace=True)
    
#     # adding information from df_orig
#     columns=["seqnames","start","end","width","strand","seq",
#             "cage_hepg2_class","cage_k562_class","dhs_hepg2_class","dhs_k562_class",
#             "starr_hepg2_class","starr_k562_class","sure_hepg2_class","sure_k562_class",
#             "cage_hepg2_intensity","cage_k562_intensity","dhs_hepg2_intensity","dhs_k562_intensity",
#             "starr_hepg2_intensity","starr_k562_intensity","sure_hepg2_intensity","sure_k562_intensity",
#             "data_type"]
    
    
#     df_orig["data_type"]=data_type
#     df_orig["start"]=df_gr["start"]
#     df_orig["end"]=df_gr["end"]
#     df_orig["width"]=600
#     df_orig=df_orig.loc[:,columns]
    
#     # merge and write as csv
#     df_res=pd.merge(df_orig, df_pred,left_index=True, right_index=True)
#     df_res["identifier"]=identifiers
#     df_res.set_index("identifier",inplace=True)
#     assert colnames==df_res.columns.to_list()
#     df_res.to_csv("/isdata/alab/people/pcr980/DeepCompare/Metadata_and_prediction/metadata_all_seqs.csv", mode='a',header=False)


#----------------------------------
# Write metrics
#----------------------------------


# df_truth=pd.read_csv("/isdata/alab/people/pcr980/DeepCompare/Datasets/Dataset_final_rep/dat_test.csv")

# df_truth_reg=df_truth.loc[:,add_suffix(create_colnames(),"_intensity")]
# df_truth_class=df_truth.loc[:,add_suffix(create_colnames(),"_class")]

# df_pred=pd.read_csv("/isdata/alab/people/pcr980/DeepCompare/Test/Pd3_DeepCompare_performance/pred_dat_test.csv",header=None)
# df_pred_reg=df_pred.iloc[:,0:8]
# df_pred_class=df_pred.iloc[:,8:16]
# df_pred_reg.columns=create_colnames()
# df_pred_class.columns=create_colnames()

# # write metrics
# df_metric=calc_columnwise_metrics(df_pred_reg,df_truth_reg,'pcc',calc_pcc)
# df_acc=calc_columnwise_metrics(df_pred_class,df_truth_class,'acc',calc_acc)
# df_metric.loc[:,"acc"]=df_acc.acc
# df_pcc_class=calc_columnwise_metrics(df_pred_class,df_truth_reg,'pcc',calc_pcc)
# df_metric.loc[:,"pcc_class"]=df_pcc_class.pcc
# df_metric.to_csv("/isdata/alab/people/pcr980/DeepCompare/Test/Pd3_DeepCompare_performance/metrics.csv")



#----------------------------------
# Plot PRC and ROC
#----------------------------------

colors_list=["tab:blue","tab:cyan","tab:orange","tab:red",
             "tab:green","tab:olive","tab:purple","tab:brown"]


# def plot_test_roc(df_truth_class,df_pred_class):
#     figure(figsize=(6, 6), dpi=80)
#     for i in range(8):
#         pred,truth=mask_neg1(df_pred_class.iloc[:,i],df_truth_class.iloc[:,i])
#         fpr,tpr,_=roc_curve(truth,pred)
#         roc_auc = auc(fpr, tpr)
#         dataset=df_pred_class.columns[i]
#         plt.plot(fpr, tpr, label = dataset+'_AUC = %0.2f' % roc_auc,color=colors_list[i])

#     plt.plot([0,1],[0,1],linestyle="dotted")
#     plt.legend(loc = 'lower right')
#     plt.title("Receiver operating characteristic")
#     plt.xlabel("False positive rate")
#     plt.ylabel("True positive rate")
#     plt.savefig("/isdata/alab/people/pcr980/DeepCompare/Test/Pd3_DeepCompare_performance/ROC.pdf",format="pdf")



# def plot_test_prc(df_truth_class,df_pred_class):
#     figure(figsize=(6, 6), dpi=80)
#     for i in range(8):
#         pred,truth=mask_neg1(df_pred_class.iloc[:,i],df_truth_class.iloc[:,i])
#         precision,recall,_=precision_recall_curve(truth,pred)
#         auprc=auc(recall,precision)
#         dataset=df_pred_class.columns[i]
#         plt.plot(recall, precision, label = dataset+'_AUC = %0.2f' % auprc,color=colors_list[i])

#     plt.plot([0,1],[1,0],linestyle="dotted")
#     plt.legend(loc = 'lower left')
#     plt.title("Precision recall curve")
#     plt.xlabel("Recall")
#     plt.ylabel("Precision")
#     plt.savefig("/isdata/alab/people/pcr980/DeepCompare/Test/Pd3_DeepCompare_performance/PRC.pdf",format="pdf")


# plot_test_roc(df_truth_class,df_pred_class)
# plot_test_prc(df_truth_class,df_pred_class)


#-------------------------------------------
# Plot PCC of classification and regression
#-------------------------------------------

df=pd.read_csv("/isdata/alab/people/pcr980/DeepCompare/Test/Pd3_DeepCompare_performance/metrics.csv")

plt.figure(figsize=(7, 4))

# Plot each point with its corresponding color, shape, and label
for i, row in df.iterrows():
    plt.scatter(row['file'], row['pcc'], color=colors_list[i], marker='D', label='Regression' if i == 0 else "")
    plt.scatter(row['file'], row['pcc_class'], color=colors_list[i], marker='o', label='Classification (z value)' if i == 0 else "")

plt.ylim(bottom=0)
plt.xticks(rotation=45)
plt.xlabel('Dataset')
plt.ylabel('Pearson correlation')
plt.title('')
plt.legend(bbox_to_anchor=(1.05, 1), loc='upper left')
plt.tight_layout()
plt.savefig("/isdata/alab/people/pcr980/DeepCompare/Test/Pd3_DeepCompare_performance/Reg_vs_class.pdf",format="pdf")


