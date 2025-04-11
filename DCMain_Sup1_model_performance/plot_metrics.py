#!/usr/bin/python3

import pandas as pd
from sklearn.metrics import roc_curve,auc,precision_recall_curve
import matplotlib.pyplot as plt
from matplotlib.pyplot import figure


import sys
sys.path.insert(1,"/isdata/alab/people/pcr980/Scripts_python")
from prediction import write_predictions 
from utils_metrics import *



import matplotlib
matplotlib.rcParams['pdf.fonttype']=42



#-------------------------------------------------------
# Task1: Generate Pd3_DeepCompare_performance/: Write DeepCompare predictions
#------------------------------------------------------
for suffix in ["train", "test", "val"]:
    write_predictions(data_path=f"/isdata/alab/people/pcr980/DeepCompare/Datasets/Dataset_final_rep/dat_{suffix}.csv",
                      seq_colname="seq",
                      out_path=f"/isdata/alab/people/pcr980/DeepCompare/Test/Pd3_DeepCompare_performance/pred_dat_{suffix}.csv")



#----------------------------------
# Task2: Write DeepCompARE metrics
#----------------------------------


df_truth=pd.read_csv("/isdata/alab/people/pcr980/DeepCompare/Datasets/Dataset_final_rep/dat_test.csv")

df_truth_reg=df_truth.loc[:,add_suffix(create_colnames(),"_intensity")]
df_truth_class=df_truth.loc[:,add_suffix(create_colnames(),"_class")]

df_pred=pd.read_csv("Pd3_DeepCompare_performance/pred_dat_test.csv",header=None)
df_pred_reg=df_pred.iloc[:,0:8]
df_pred_class=df_pred.iloc[:,8:16]
df_pred_reg.columns=create_colnames()
df_pred_class.columns=create_colnames()

# write metrics
df_metric=calc_columnwise_metrics(df_pred_reg,df_truth_reg,'pcc',calc_pcc)
df_acc=calc_columnwise_metrics(df_pred_class,df_truth_class,'acc',calc_acc)
df_metric.loc[:,"acc"]=df_acc.acc
df_pcc_class=calc_columnwise_metrics(df_pred_class,df_truth_reg,'pcc',calc_pcc)
df_metric.loc[:,"pcc_class"]=df_pcc_class.pcc
df_metric.to_csv("Pd3_DeepCompare_performance/metrics.csv")



#----------------------------------
# Plot PRC and ROC
#----------------------------------
color_list = [
    "#8B0000", "#F08080",  # Dark Red and Light Red
    "#4169E1", "#ADD8E6",  # Lighter Dark Blue (Royal Blue) and Light Blue
    "#006400", "#90EE90",  # Dark Green and Light Green
    "#8B008B", "#D8BFD8"   # Redder Dark Purple (Dark Magenta) and Light Purple
]

def plot_roc(df_truth_class,df_pred_class):
    figure(figsize=(3, 3))
    # thin frame
    ax = plt.gca()
    ax.spines['top'].set_linewidth(0.5)
    ax.spines['bottom'].set_linewidth(0.5)
    ax.spines['left'].set_linewidth(0.5)
    ax.spines['right'].set_linewidth(0.5)
    # Set tick parameters to make ticks thinner
    ax.tick_params(axis='both', which='both', width=0.5)
    for i in range(8):
        pred,truth=mask_neg1(df_pred_class.iloc[:,i],df_truth_class.iloc[:,i])
        fpr,tpr,_=roc_curve(truth,pred)
        roc_auc = auc(fpr, tpr)
        dataset=df_pred_class.columns[i]
        plt.plot(fpr, tpr, label = dataset+': %0.2f' % roc_auc,color=color_list[i],linewidth=0.5)
    #
    plt.plot([0,1],[0,1],linestyle="dotted",color="black",linewidth=0.5)
    plt.legend(loc = 'lower right',fontsize=5, title_fontsize=5, title="Area under curve")
    plt.title("Classification performance:\nReceiver operating characteristic",fontsize=7)
    plt.xlabel("False positive rate",fontsize=7)
    plt.ylabel("True positive rate",fontsize=7)
    plt.xticks(fontsize=5)
    plt.yticks(fontsize=5)
    plt.tight_layout()
    plt.savefig("roc.pdf",dpi=300)
    plt.close()


plot_roc(df_truth_class,df_pred_class)





def plot_prc(df_truth_class,df_pred_class):
    figure(figsize=(3, 3))
    # thin frame
    ax = plt.gca()
    ax.spines['top'].set_linewidth(0.5)
    ax.spines['bottom'].set_linewidth(0.5)
    ax.spines['left'].set_linewidth(0.5)
    ax.spines['right'].set_linewidth(0.5)
    # Set tick parameters to make ticks thinner
    ax.tick_params(axis='both', which='both', width=0.5)
    for i in range(8):
        pred,truth=mask_neg1(df_pred_class.iloc[:,i],df_truth_class.iloc[:,i])
        precision,recall,_=precision_recall_curve(truth,pred)
        auprc=auc(recall,precision)
        dataset=df_pred_class.columns[i]
        plt.plot(recall, precision, label = dataset+': %0.2f' % auprc,color=color_list[i],linewidth=0.5)
    #
    plt.plot([0,1],[1,0],linestyle="dotted",color="black",linewidth=0.5)
    plt.legend(loc = 'lower left',fontsize=5, title_fontsize=5, title="Area under curve")
    plt.title("Classification performance:\nPrecision recall curve",fontsize=7)
    plt.xlabel("Recall",fontsize=7)
    plt.ylabel("Precision",fontsize=7)
    plt.xticks(fontsize=5)
    plt.yticks(fontsize=5)
    plt.tight_layout()
    plt.savefig("prc.pdf",dpi=300)
    plt.close()


plot_prc(df_truth_class,df_pred_class)


#-------------------------------------------
# Plot PCC of classification and regression
#-------------------------------------------

df=pd.read_csv("Pd3_DeepCompare_performance/metrics.csv")



figure(figsize=(3, 2.5))
# thin frame
ax = plt.gca()
ax.spines['top'].set_linewidth(0.5)
ax.spines['bottom'].set_linewidth(0.5)
ax.spines['left'].set_linewidth(0.5)
ax.spines['right'].set_linewidth(0.5)
# Set tick parameters to make ticks thinner
ax.tick_params(axis='both', which='both', width=0.5)
# Plot each point with its corresponding color, shape, and label
for i, row in df.iterrows():
    plt.scatter(row['file'], row['pcc'], color=color_list[i], s=5)
    # plt.scatter(row['file'], row['pcc_class'], facecolors='none', edgecolors=color_list[i], marker='o')

plt.ylim(bottom=0,top=1)
plt.xticks(rotation=45,fontsize=5)
plt.yticks(fontsize=5)
plt.xlabel('Dataset',fontsize=7)
plt.ylabel('Pearson correlation',fontsize=7)
plt.title("Regression Performance",fontsize=7)
# add invisible dots for legend
# plt.scatter([], [], facecolors='none', edgecolors='black', label='Regression', s=5)
# plt.scatter([], [], facecolors='none', edgecolors='black', marker='o', label='Classification (z value)')
# plt.legend(fontsize=5)
plt.tight_layout()
plt.savefig("pcc.pdf",dpi=300)
plt.close()


