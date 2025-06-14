import matplotlib
import matplotlib.pyplot as plt
import matplotlib.patches as mpatches
import pandas as pd
import matplotlib.ticker as ticker
from loguru import logger
from scipy.stats import pearsonr

import sys
sys.path.insert(1,"/isdata/alab/people/pcr980/Scripts_python")
from seq_ops import SeqExtractor
from region_ops import resize_region
from seq_annotators import JasparAnnotator
from in_silico_mutagenesis import compute_mutagenesis_score, get_motif_isa


matplotlib.rcParams['pdf.fonttype']=42

#---------------
# Load annotators
#---------------
seq_extractor=SeqExtractor("/isdata/alab/people/pcr980/Resource/hg38.fa")
jaspar_hepg2_annotator=JasparAnnotator("/isdata/alab/people/pcr980/Resource/JASPAR2022_tracks/JASPAR2022_hg38.bb",
                                       chip_file="/isdata/alab/people/pcr980/Resource/ReMap2022/ReMap2022_hg38_hepg2.bed",
                                       rna_file="/isdata/alab/people/pcr980/DeepCompare/RNA_expression/expressed_tf_list_hepg2.tsv"
                                       )
jaspar_k562_annotator=JasparAnnotator("/isdata/alab/people/pcr980/Resource/JASPAR2022_tracks/JASPAR2022_hg38.bb",
                                       chip_file="/isdata/alab/people/pcr980/Resource/ReMap2022/ReMap2022_hg38_k562.bed",
                                       rna_file="/isdata/alab/people/pcr980/DeepCompare/RNA_expression/expressed_tf_list_k562.tsv"
                                       )

#---------------
# Helper functions
#---------------



def get_motifs(seq_extractor,jaspar_annotator,region,track_num,score_threshold):
    df_motif=jaspar_annotator.annotate(region)
    # remove uncertain proteins
    df_motif=df_motif[~df_motif["protein"].str.contains("::")].copy().reset_index(drop=True)
    df_chip=df_motif.loc[df_motif["chip_evidence"]==True,:].reset_index(drop=True)
    df_rest=df_motif.loc[df_motif["chip_evidence"]==False,:].reset_index(drop=True)
    df_rna=df_rest.loc[df_rest["rna_evidence"]==True,:].reset_index(drop=True)
    df_rna=df_rna.loc[df_rna["score"]>=score_threshold,:].reset_index(drop=True)
    df_motif=pd.concat([df_chip,df_rna],axis=0).reset_index(drop=True)
    # sort df_motif by "start"
    df_motif=df_motif.sort_values(by="start").reset_index(drop=True)
    # get motif isa
    df_motif=get_motif_isa(seq_extractor,df_motif,track_num)
    df_motif.rename(columns={f"isa_track{track_num}":"isa"},inplace=True)
    return df_motif



def reduce_protein_names(protein_list):
    protein_list=list(set(protein_list))
    # order alphabetically
    protein_list.sort()
    # if there are more than 2 proteins sharing same prefix of length > 4
    # only keep the prefix, followed by "s"
    # eg: hoxa9, hoxa9b, hoxa9c -> hoxa9s
    protein_dict={}
    for protein in protein_list:
        prefix=protein[:4]
        if prefix in protein_dict:
            protein_dict[prefix].append(protein)
        else:
            protein_dict[prefix]=[protein]
    prefix_list=[]
    for prefix in protein_dict:
        if len(protein_dict[prefix])>1:
            prefix_list.append(prefix+"s")
        else:
            prefix_list.append(protein_dict[prefix][0])
    # concatenate by "\n"
    return "\n".join(prefix_list)




def reduce_motifs(df_motif,window=4):
    # if start is within 3bp of another start
    # choose the top 3 based on "score"
    # concatenate protein with "\n", use largest "end" as end
    df_res=pd.DataFrame()
    while df_motif.shape[0]>0:
        current_start=df_motif.loc[0,"start_rel"]
        df_temp=df_motif[(df_motif["start_rel"]>=current_start) & (df_motif["start_rel"]<=(current_start+window))].copy().reset_index(drop=True)
        df_temp=df_temp.sort_values(by="score",ascending=False).reset_index(drop=True)
        df_temp=df_temp.iloc[:2,:]
        df_temp["protein"]=reduce_protein_names(df_temp["protein"])
        df_temp["isa"]=df_temp["isa"].mean()
        df_temp["end_rel"]=df_temp["end_rel"].max()
        df_temp["start_rel"]=df_temp["start_rel"].min()
        df_temp["start"]=df_temp["start"].min()
        df_temp["end"]=df_temp["end"].max()
        df_res=df_res.append(df_temp.iloc[0,:],ignore_index=True)
        # remove the rows in df_temp from df_motif
        df_motif=df_motif[df_motif["start_rel"]>current_start+window].copy().reset_index(drop=True)
    return df_res



def plot_motif_imp(df_motif,ax,ylim):
    # only relative position matters
    # Hide spines
    for spine in ax.spines.values():
        spine.set_visible(False)
    # Plot ISA score
    ax.bar(df_motif["start_rel"], df_motif["isa"], width=df_motif["end_rel"] - df_motif["start_rel"], color="#1f77b4", alpha=0.5, align='edge')
    ax.axhline(0, color='black', lw=0.5)
    # Add text labels for proteins
    prev_text_pos=0
    for idx, row in df_motif.iterrows():
        current_text_pos=row["start_rel"]
        if current_text_pos-prev_text_pos<5:
            current_text_pos=prev_text_pos+5
            if current_text_pos>row["end_rel"]:
                pass
                # raise ValueError("Annotation overlap cannot be resolved")
        ax.text(current_text_pos,row["isa"], row["protein"], rotation=90, fontsize=5)
        prev_text_pos=current_text_pos
    #  set labels
    ax.set_title("Motif ISA",fontsize=7)
    ax.tick_params(axis='y', which='major', labelsize=5)
    ax.set_xticks([])
    ax.set_ylim(ylim)



def plot_base_imp(df,ax,markersize=1):
    """
    df have columns "position", "base", "imp"
    """
    # Hide spines
    for spine in ax.spines.values():
        spine.set_visible(False)
    ax.bar(df["position"], df["imp"], color="#1f77b4", alpha=0.5)
    ax.axhline(0, color='black', lw=0.5)
    # Create color map for bases
    color_dict = {"A": "#1f77b4", "C": "#ff7f0e", "G": "#2ca02c", "T": "#d62728"}
    # Plot markers for each base
    for row in df.itertuples():
        ax.plot(row.position, row.imp, marker="o", color=color_dict[row.base], markersize=markersize)




def plot_region_4panels(seq_extractor,jaspar_annotator, df_truth, element_name, region, track_num, score_threshold, markersize):
    region_resized=resize_region(region,599,fix="center")
    left_shift=region[1]-region_resized[1]
    # get seq
    seq=seq_extractor.get_seq(region)
    seq_resized=seq_extractor.get_seq(region_resized)
    # get base importance
    isa=compute_mutagenesis_score(seq_resized,"isa","mean").loc[f"Seq0_Track{track_num}",:]
    isa=isa[left_shift:(left_shift+region[2]-region[1]+1)].reset_index(drop=True)
    ism=compute_mutagenesis_score(seq_resized,"ism","mean").loc[f"Seq0_Track{track_num}",:]
    ism=ism[left_shift:(left_shift+region[2]-region[1]+1)].reset_index(drop=True)
    # get df_motif
    df_motif=get_motifs(seq_extractor,jaspar_annotator,region_resized,track_num,score_threshold)
    df_motif=reduce_motifs(df_motif)
    # subset to region
    df_motif=df_motif[(df_motif.loc[:,"start"]>=region[1]) & (df_motif.loc[:,"end"]<=region[2])].copy().reset_index(drop=True)
    df_motif["start_rel"]-=left_shift
    df_motif["end_rel"]-=left_shift
    ymax=max(max(ism),max(isa))
    ymin=min(min(ism),min(isa))
    fig, (ax1, ax2, ax3, ax4) = plt.subplots(4, 1, sharex=True, figsize = (180/25.4, 110/25.4))
    # ax1: truth
    df_truth=df_truth.loc[df_truth["Element"]==element_name,:].reset_index(drop=True)
    df_truth["position"]=df_truth["position"]-region[1]
    plot_base_imp(df_truth,ax1)
    ax1.xaxis.set_major_locator(ticker.MaxNLocator(nbins=6))
    ax1.tick_params(axis='both', which='major', labelsize=5)
    ax1.set_title("MPRA experiment",fontsize=7)
    # ax2: base level ism
    df_a=pd.DataFrame({"position":list(range(len(ism))),"base":list(seq),"imp":ism})
    plot_base_imp(df_a,ax2)
    ax2.set_title("ISM (Average of 3 alternative bases)",fontsize=7)
    ax2.tick_params(axis='y', which='major', labelsize=5)
    ax2.set_ylim((ymin,ymax))
    # ax3: base level isa
    df_n= pd.DataFrame({"position":list(range(len(isa))),"base":list(seq),"imp":isa})
    plot_base_imp(df_n,ax3)
    ax3.set_title("ISA (Base replaced by N)",fontsize=7)
    ax3.tick_params(axis='y', which='major', labelsize=5)
    ax3.set_ylim((ymin,ymax))
    # Plot motif ISA on the last axis
    plot_motif_imp(df_motif, ax4,(ymin,ymax))
    # add legend
    color_dict = {"A": "#1f77b4", "C": "#ff7f0e", "G": "#2ca02c", "T": "#d62728"}
    handles = [mpatches.Patch(color=color, label=base) for base, color in color_dict.items()]
    ax4.legend(handles=handles, title="Bases", title_fontsize=5, fontsize=5, loc="upper right")
    fig.text(0.5, 0.05, "Relative position", ha='center', fontsize=7)
    fig.text(0.05, 0.5, "Feature importance (base or motif)", va='center', rotation='vertical', fontsize=7)
    plt.savefig(f"{element_name}_track{track_num}.pdf",dpi=300)
    plt.close()
    # calculate correlation between truth and ism
    df_merged=pd.merge(df_truth, df_a, on='position', suffixes=('_truth', '_ism'))
    df_merged=df_merged.dropna()
    r,p=pearsonr(df_merged.imp_truth, df_merged.imp_ism)
    logger.info(f"correlation between truth and ISM: {r}")
    # calculate correlation between truth and isa
    df_merged=pd.merge(df_truth, df_n, on='position', suffixes=('_truth', '_isa'))
    df_merged=df_merged.dropna()
    r,p=pearsonr(df_merged.imp_truth, df_merged.imp_isa)
    logger.info(f"correlation between truth and ISA: {r}")









def plot_region_2panels(seq_extractor,df_truth, element_name, region, track_num,markersize):
    region_resized=resize_region(region,599,fix="center")
    left_shift=region[1]-region_resized[1]
    # get seq
    seq=seq_extractor.get_seq(region)
    seq_resized=seq_extractor.get_seq(region_resized)
    # get base importance
    isa=compute_mutagenesis_score(seq_resized,"isa","mean").loc[f"Seq0_Track{track_num}",:]
    isa=isa[left_shift:(left_shift+region[2]-region[1]+1)].reset_index(drop=True)
    fig, (ax1, ax2) = plt.subplots(2, 1, sharex=True, figsize = (130/25.4, 60/25.4))
    # plot isa on the first axis
    df_n= pd.DataFrame({"position":list(range(len(isa))),"base":list(seq),"imp":isa})
    plot_base_imp(df_n,ax1,markersize)
    # customize the first axis
    ax1.set_title(element_name,fontsize=7)
    ax1.set_xlabel("ISA (Base replaced by N)",fontsize=7,labelpad=0)
    ax2.set_yticks([0, 0.25])
    ax1.tick_params(axis='both', which='major', labelsize=7)
    # plot truth on the second axis
    df_truth=df_truth.loc[df_truth["Element"]==element_name,:].reset_index(drop=True)
    df_truth["position"]=df_truth["position"]-region[1]
    plot_base_imp(df_truth,ax2,markersize)
    # add supra ticks and labels
    ax2.set_xlabel("MPRA experiment",fontsize=7,labelpad=0)
    ax2.set_yticks([0, 2.5])
    ax2.tick_params(axis='both', which='major', labelsize=7)
    # add legend
    color_dict = {"A": "#1f77b4", "C": "#ff7f0e", "G": "#2ca02c", "T": "#d62728"}
    handles = [mpatches.Patch(color=color, label=base) for base, color in color_dict.items()]
    ax2.legend(handles=handles, title="Bases", title_fontsize=5, fontsize=5, loc="upper right")
    fig.text(0.5, 0.01, "Relative position", ha='center', fontsize=7)
    plt.savefig(f"{element_name}_track{track_num}.pdf",dpi=300)
    plt.close()










# truth
df_truth=pd.read_csv("GRCh38_TERT-GAa_TERT-GSc_TERT-GBM_TERT-HEK_LDLR_F9_PKLR-24h_SORT1.tsv",sep="\t")
# select P-Value<0.05
df_truth=df_truth.loc[df_truth["P-Value"]<0.05,:].reset_index(drop=True)
# group by Position, Ref, Element, and select the largest Value
# method 1: aggregate by mean
df_truth=df_truth.groupby(["Position","Ref","Element"],as_index=False).agg({"Value":"mean"})
# method 2: aggregate by max(key=abs)
# df_truth = df_truth.groupby(["Position", "Ref", "Element"], as_index=False).agg({"Value": lambda x: x.loc[x.abs().idxmax()]})
df_truth=df_truth.rename(columns={"Position":"position","Value":"imp","Ref":"base"})
df_truth["imp"]=-df_truth["imp"]
df_truth.Element.unique()




# F9 promoter: chrX:139530463-139530765
plot_region_4panels(seq_extractor,jaspar_hepg2_annotator,df_truth,"F9",("chrX",139530463,139530765),6,360,1) # window 4, text space 5


# LDLR promoter: chr19:11,089,231-11,089,548
# plot_region_2panels(seq_extractor,jaspar_hepg2_annotator,df_truth,"LDLR",("chr19",11089231,11089548),6,500,1) 
# PKLR promoter: chr1:155,301,395-155,301,864
# plot_region_2panels(seq_extractor,df_truth,"PKLR-24h",("chr1",155301395,155301864),7,0.8)






#----------------------------------
# Understand SNV on SORT1 enhancer
#----------------------------------

# SORT1 enhancer: chr1:109,274,652-109,275,251


def plot_ism_ref_vs_mut(imp_ref,imp_mut,title_prefix):
    fig, (ax0, ax1) = plt.subplots(2, 1, sharex=True, figsize=(160/25.4, 70/25.4))
    # plot ism_ref on the first axis
    df_ref= pd.DataFrame({"position":list(range(len(imp_ref))),"base":list(seq_ref),"imp":imp_ref})
    plot_base_imp(df_ref,ax0)
    ax0.add_patch(mpatches.Rectangle((315, -0.1), 9, 0.2, edgecolor='red', facecolor='none', lw=1))
    ax0.set_title(f"{title_prefix} (reference)", fontsize=7)
    ax0.tick_params(axis='y', which='major', labelsize=7)
    ax0.set_yticks([-0.05, 0, 0.05])
    # plot ism_mut on the second axis
    df_mut=pd.DataFrame({"position":list(range(len(imp_mut))),"base":list(seq_ref),"imp":imp_mut})
    plot_base_imp(df_mut,ax1)
    # add a red box to position 315-323
    ax1.add_patch(mpatches.Rectangle((315, -0.1), 9, 0.2, edgecolor='red', facecolor='none', lw=1))
    ax1.xaxis.set_major_locator(ticker.MaxNLocator(nbins=6))
    ax1.tick_params(axis='both', which='major', labelsize=7)
    ax1.set_title(f"{title_prefix} (109274968:G>T)", fontsize=7)
    ax1.set_yticks([-0.05, 0, 0.05])
    # add legend
    color_dict = {"A": "#1f77b4", "C": "#ff7f0e", "G": "#2ca02c", "T": "#d62728"}
    handles = [mpatches.Patch(color=color, label=base) for base, color in color_dict.items()]
    ax1.legend(handles=handles, title="Bases", title_fontsize=5, fontsize=5, loc="upper right")
    # add supra labels and title
    fig.text(0.5, 0.05, "Relative position", ha='center', fontsize=7)
    fig.text(0.05, 0.5, "Base importance", va='center', rotation='vertical', fontsize=7)
    plt.suptitle("SORT1 enhancer", fontsize=7)
    plt.tight_layout(rect=[0, 0.01, 1, 0.95])
    plt.savefig(f"SORT1_enhancer_mutation_{title_prefix}.pdf",dpi=300)
    plt.close()



region=("chr1",109274652,109275251)
seq_ref=seq_extractor.get_seq(region)

# get relative position of the mutation 1-109274968-G-T
mut_pos=109274968-109274652
seq_ref[mut_pos]
seq_mut=seq_ref[:mut_pos]+"T"+seq_ref[(mut_pos+1):]

# calculate ism and isa
track_num=4
# select first 400bp
isa_ref=compute_mutagenesis_score(seq_ref,"isa","mean").loc[f"Seq0_Track{track_num}",:399]
ism_ref=compute_mutagenesis_score(seq_ref,"ism","mean").loc[f"Seq0_Track{track_num}",:399]
isa_mut=compute_mutagenesis_score(seq_mut,"isa","mean").loc[f"Seq0_Track{track_num}",:399]
ism_mut=compute_mutagenesis_score(seq_mut,"ism","mean").loc[f"Seq0_Track{track_num}",:399]
seq_ref=seq_ref[:400]
seq_mut=seq_mut[:400]
# plot


# plot_ism_ref_vs_mut(ism_ref,ism_mut,"ISM")
plot_ism_ref_vs_mut(isa_ref,isa_mut,"ISA")




