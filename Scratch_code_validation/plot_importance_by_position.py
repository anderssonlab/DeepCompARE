import matplotlib
import matplotlib.pyplot as plt
import matplotlib.patches as mpatches
import pandas as pd

import sys
sys.path.insert(1,"/isdata/alab/people/pcr980/Scripts_python")
from seq_ops import SeqExtractor
from in_silico_mutagenesis import compute_mutagenesis_score, get_motif_isa
from seq_annotators import JasparAnnotator, gnomadAnnotator, phylopAnnotator, baseImpAnnotator

matplotlib.rcParams['pdf.fonttype']=42

#---------------
# Load annotators
#---------------
seq_extractor=SeqExtractor("/isdata/alab/people/pcr980/Resource/hg38.fa")
jaspar_annotator=JasparAnnotator("/isdata/alab/people/pcr980/Resource/JASPAR2022_tracks/JASPAR2022_hg38.bb",
                                       chip_file="/isdata/alab/people/pcr980/Resource/ReMap2022/ReMap2022_hg38_hepg2.bed",
                                       rna_file="/isdata/alab/people/pcr980/DeepCompare/RNA_expression/expressed_tf_list_hepg2.tsv"
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
    # if there are more than 3 proteins sharing same prefix of length > 4
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



def plot_motif_imp(df_motif, ax,ylim=None):
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
    # Set title and labels
    ax.set_title("Motif ISA score",fontsize=7)
    ax.tick_params(axis='y', which='major', labelsize=7)
    if ylim:
        ax.set_ylim(ylim)



def plot_base_imp(df,ax,title,ylim=None,xlabel=False):
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
        ax.plot(row.position, row.imp, marker="o", color=color_dict[row.base], markersize=1) # 1 for F9, 2 for SORT1
    # Set title and labels
    ax.set_title(title,fontsize=7)
    if xlabel:
        ax.set_xlabel("Position",fontsize=7)
    ax.tick_params(axis='both', which='major', labelsize=7)
    #
    # Create legend
    handles = [mpatches.Patch(color=color, label=base) for base, color in color_dict.items()]
    ax.legend(handles=handles, title="Bases", title_fontsize=5, fontsize=5, loc="upper right")
    if ylim:
        ax.set_ylim(ylim)


# everything calculated fresh
region=("chr18",49492297,49492896)
track_num=0
score_threshold=500
# get seq
seq=seq_extractor.get_seq(region)
# get base importance
isa=compute_mutagenesis_score(seq,"isa","mean").loc[f"Seq0_Track{track_num}",:]
isa=isa[0:(region[2]-region[1]+1)].reset_index(drop=True)
ism=compute_mutagenesis_score(seq,"ism","mean").loc[f"Seq0_Track{track_num}",:]
ism=ism[0:(region[2]-region[1]+1)].reset_index(drop=True)
# get df_motif
df_motif=get_motifs(seq_extractor,jaspar_annotator,region,track_num,score_threshold)
df_motif=reduce_motifs(df_motif)
# subset to region
df_motif=df_motif[(df_motif.loc[:,"start"]>=region[1]) & (df_motif.loc[:,"end"]<=region[2])].copy().reset_index(drop=True)

ymax=max(max(ism),max(isa))
ymin=min(min(ism),min(isa))
fig, (ax1, ax2, ax3) = plt.subplots(3, 1, sharex=True, figsize = (180/25.4, 120/25.4))
# Plot ISA score on the first axis
plot_motif_imp(df_motif, ax1,ylim=(ymin,ymax))
# plot isa on the second axis
df_n= pd.DataFrame({"position":list(range(len(isa))),"base":list(seq),"imp":isa})
plot_base_imp(df_n,ax2,"ISA (Base replaced by N)",ylim=(ymin,ymax))
df_a=pd.DataFrame({"position":list(range(len(ism))),"base":list(seq),"imp":ism})
plot_base_imp(df_a,ax3,"ISM (Average of 3 alternative bases)",ylim=(ymin,ymax))
plt.tight_layout()
plt.savefig(f"example.pdf",dpi=300)
plt.close()





# what's written in the file
df_file=pd.read_csv("/isdata/alab/people/pcr980/DeepCompare/Pd5_motif_info/motif_info_thresh_500_enhancers_hepg2.csv")
# select region
df_file=df_file[df_file["region"]=="chr18:49492297-49492896"].reset_index(drop=True)


df_file.protein
gradxinp_annotator=baseImpAnnotator(f"/isdata/alab/people/pcr980/DeepCompare/Pd4_promoters_enhancers_and_featimp/gradxinp_enhancers_hepg2.csv")
ism_annotator=baseImpAnnotator(f"/isdata/alab/people/pcr980/DeepCompare/Pd4_promoters_enhancers_and_featimp/ism_enhancers_hepg2.csv")
isa_annotator=baseImpAnnotator(f"/isdata/alab/people/pcr980/DeepCompare/Pd4_promoters_enhancers_and_featimp/isa_enhancers_hepg2.csv")



df_enhancers_hepg2=pd.read_csv("/isdata/alab/people/pcr980/DeepCompare/Pd4_promoters_enhancers_and_featimp/enhancers_hepg2.bed",sep="\t",header=None)
# check which regio is chr18:49492297-49492896
df_enhancers_hepg2[(df_enhancers_hepg2[0]=="chr18") & (df_enhancers_hepg2[1]==49492297) & (df_enhancers_hepg2[2]==49492897)]