import pandas as pd
import os
import pyBigWig
from region_ops import subset_df_by_region
from loguru import logger



def find_binding_evidence(motif_df, protein_col, tf_list, new_col_name):
    """
    motif_df: A DataFrame with containing motifs spanning a region
    tf_list: A list of transcription factors
    new_col_name: A string of the new column name
    """
    if len(tf_list)==0:
        motif_df[new_col_name] = False
        return motif_df
    elif motif_df.shape[0]==0:
        motif_df[new_col_name] = pd.Series(dtype='object')
        return motif_df
    else:
        protein1 = motif_df[protein_col].str.split("::", expand=True).iloc[:, 0]
        protein2 = motif_df[protein_col].str.split("::", expand=True).iloc[:, -1]
        protein2 = protein2.fillna(protein1)
        # require both proteins to be in the list
        motif_df[new_col_name] = True
        motif_df.loc[~protein1.isin(tf_list), new_col_name] = False
        motif_df.loc[~protein2.isin(tf_list), new_col_name] = False
    return motif_df



class rnaAnnotator:
    def __init__(self, reference):
        """
        Annotate whether each row in motif_df have binding evidence
        reference: file path
            "/isdata/alab/people/pcr980/DeepCompare/RNA_expression/expressed_tf_list_hepg2.tsv"
            "/isdata/alab/people/pcr980/DeepCompare/RNA_expression/expressed_tf_list_k562.tsv
        """
        self.df_rna = pd.read_csv(reference, sep='\t', header=None)
        self.df_rna.columns = ['TF']
    #
    #
    def annotate(self, motif_df):
        """
        motif_df: A DataFrame with containing motifs spanning a region
        Return: motif_df with a new column "rna_evidence"
        """
        motif_df=find_binding_evidence(motif_df,"protein",self.df_rna['TF'].values, 'rna_evidence')
        return motif_df
        






class chipAnnotator:
    def __init__(self, reference):
        """
        Annotate whether each row in motif_df have binding evidence
        reference: file path, see below
            "/isdata/alab/people/pcr980/Resource/ReMap2022/ReMap2022_hg38_hepg2.bed",
            "/isdata/alab/people/pcr980/Resource/ReMap2022/ReMap2022_hg38_k562.bed",
        Return: True or False for each row
        """
        self.df_chip = pd.read_csv(reference, sep='\t', header=None)
        self.df_chip.columns = ['chromosome', 'start', 'end', 'detail', 'x', 'strand', 'thick_start', 'thick_end', 'score']
        self.df_chip['TF'] = self.df_chip['detail'].str.split(':', expand=True).iloc[:,0].str.upper()
        self.df_chip.drop(columns=['detail', 'x', 'thick_start', 'thick_end'], inplace=True)
    #
    #
    def annotate(self, motif_df):
        """
        motif_df: A DataFrame output by JasparAnnotator.annotate()
        """
        # if motif_df is empty, return the dataframe with a new column "chip_evidence"
        if motif_df.shape[0]==0:
            motif_df['chip_evidence'] = pd.Series(dtype='object')
            return motif_df
        # get unique regions, retain order
        df_res=pd.DataFrame()
        for region in motif_df["region"].unique():
            # subset df_chip by region
            chr=region.split(":")[0]
            start=int(region.split(":")[1].split("-")[0])
            end=int(region.split(":")[1].split("-")[1])
            df_chip_subset=subset_df_by_region(self.df_chip,(chr,start,end),by="1bp")
            # subset motif_df by region
            motif_df_subset=motif_df[motif_df["region"]==region].reset_index(drop=True)
            motif_df_subset=find_binding_evidence(motif_df_subset, "protein",set(df_chip_subset["TF"].values), 'chip_evidence')
            df_res=pd.concat([df_res,motif_df_subset],ignore_index=True)
            return df_res













class JasparAnnotator:
    def __init__(self, 
                 jaspar_path, 
                 by="contained", 
                 score_thresh=0,
                 chip_file=None, 
                 rna_file=None,
                 subset="rna"):
        """
        Get Jaspar motif information from a genomic region
        by=
            "1bp": return motifs that overlap at least 1bp with the region 
            "contained": return motifs that are completely contained in the region
            "reverse": return motifs that do not overlap with the region
        chip_file=
            "/isdata/alab/people/pcr980/Resource/ReMap2022/ReMap2022_hg38_hepg2.bed",
            "/isdata/alab/people/pcr980/Resource/ReMap2022/ReMap2022_hg38_k562.bed",
        rna_file=
            "/isdata/alab/people/pcr980/DeepCompare/RNA_expression/expressed_tf_list_hepg2.tsv"
            "/isdata/alab/people/pcr980/DeepCompare/RNA_expression/expressed_tf_list_k562.tsv"
        subset=
            None: do not subset
            "rna": subset by rna_evidence==True
            "chip": subset by chip_evidence==True
        
        """
        self.jaspar_path = jaspar_path
        self.jaspar = pyBigWig.open(jaspar_path)
        self.by = by
        self.score_thresh = score_thresh
        self.chip_file=chip_file
        self.rna_file=rna_file
        self.subset=subset
        if chip_file is not None:
            self.be_chip = chipAnnotator(chip_file)
        if rna_file is not None:
            self.be_rna = rnaAnnotator(rna_file)

        #
        #
    def _annotate_region(self, region):
        """
        Args:
            region: A tuple of (chromosome, start, end).
        Returns:
            A DataFrame of chrom, start, end, score, strand and protein
        """
        start = region[1]
        end = region[2]
        # Retrieve entries from the JASPAR file
        try:
            tf_info_list = self.jaspar.entries(region[0], start, end)
        except:
            logger.warning(f"Region {region} is not in the JASPAR file")
            return pd.DataFrame()
        df = pd.DataFrame(tf_info_list, columns=['start', 'end', 'details'])
        df['chromosome'] = region[0]
        df = df.loc[:, ['chromosome', 'start', 'end', 'details']]
        df = subset_df_by_region(df, (region[0], start, end), by=self.by)
        if df.shape[0] == 0:
            return df
        # Parse details based on the genome version
        if "hg38" in self.jaspar_path:
            df[['protein', 'score', 'strand']] = df['details'].str.split('\t', expand=True)
            df.drop(columns='details', inplace=True)
        elif "hg19" in self.jaspar_path:
            df[['motif_id', 'score', 'strand', 'protein']] = df['details'].str.split('\t', expand=True)
            df.drop(columns=['details', 'motif_id'], inplace=True)
        # Process the protein and score columns
        df['protein'] = df['protein'].str.upper()
        df['region'] = f"{region[0]}:{region[1]}-{region[2]}"
        df['score'] = df['score'].astype(int)
        df = df.drop_duplicates(subset=['chromosome', 'start', 'end', 'strand', 'protein']).reset_index(drop=True)
        # Filter by score threshold
        df = df[df['score'] > self.score_thresh].reset_index(drop=True)
        return df
    #
    #
    #
    def _add_or_subset_binding_evidence(self,df):
        if self.chip_file is not None:
            df = self.be_chip.annotate(df)
            if self.subset=="chip":
                df=df[df["chip_evidence"]].reset_index(drop=True)
        if self.rna_file is not None:
            df = self.be_rna.annotate(df)
            if self.subset=="rna":
                df=df[df["rna_evidence"]].reset_index(drop=True)
        return df
    #
    #
    #
    def annotate(self, regions, outpath=None):
        """
        Args:
            region(s): A tuple or list of (chromosome, start, end) or a DataFrame of regions
            outpath: Optional. A string path to write the results to
        """
        # if region is a tuple or list, call annotate_region
        if isinstance(regions, tuple) or isinstance(regions, list):
            df=self._annotate_region(regions)
            df=self._add_or_subset_binding_evidence(df)
            return df
        # if region is a data frame, call annotate_region for each tuple and 
        # if outpath is None, concatenate the results and return 
        # if outpath is not None, write the results to outpath in mode='a'
        if isinstance(regions, pd.DataFrame) and outpath is None:
            df_res = pd.DataFrame()
            for idx, region in regions.iterrows():
                df = self._annotate_region((region[0], region[1], region[2]))
                df["seq_idx"]=f"Seq{idx}"
                df_res = pd.concat([df_res, df], ignore_index=True)
            df_res=self._add_or_subset_binding_evidence(df_res)
            return df_res
        
        if isinstance(regions, pd.DataFrame) and outpath is not None:
            for idx, region in regions.iterrows():
                df = self._annotate_region((region[0], region[1], region[2]))
                df["seq_idx"]=f"Seq{idx}"
                df=self._add_or_subset_binding_evidence(df)
                if not os.path.isfile(outpath):
                    df.to_csv(outpath,index=False,mode="w",header=True)
                else:
                    df.to_csv(outpath,index=False,mode="a",header=False)











def read_maf_file(chrom_num):
    df_maf=pd.read_csv(f"/isdata/alab/people/pcr980/Resource/Gnomad_vcf/Pd1_MAF/chrom{chrom_num}.csv")
    df_maf=df_maf[df_maf["REF"].str.len()==1]
    df_maf=df_maf[df_maf["ALT"].str.len()==1]
    df_maf=df_maf.reset_index(drop=True)
    df_maf["End"]=df_maf["POS"]+1
    df_maf.columns=["Chromosome","Start","ID","REF","ALT","AF","End"]
    df_maf=df_maf[["Chromosome","Start","End","ID","REF","ALT","AF"]]
    return df_maf
    
    




class gnomadAnnotator:
    def __init__(self,chrom_num):
        # Read MAF of all chromosomes into one dictionary
        self.maf_file=read_maf_file(chrom_num)
    #
    def _annotate_one_motif(self, motif,maf_subset_by_region):
        maf_subset_by_motif=subset_df_by_region(maf_subset_by_region,motif,by="contained")
        res=":".join(maf_subset_by_motif["AF"].astype(str).tolist())
        return res
    #
    def annotate(self, motif_df):
        df_res=pd.DataFrame()
        for region in motif_df.region.unique():
            chr=region.split(":")[0]
            start=int(region.split(":")[1].split("-")[0])
            end=int(region.split(":")[1].split("-")[1])
            maf_subset_by_region=subset_df_by_region(self.maf_file,(chr,start,end),by="contained")
            motif_df_subset=motif_df[motif_df["region"]==region].reset_index(drop=True)
            motif_df_subset["gnomad_af"]=motif_df_subset.apply(lambda row: self._annotate_one_motif((row["chromosome"],row["start"],row["end"]),maf_subset_by_region),axis=1)
            df_res=pd.concat([df_res,motif_df_subset],ignore_index=True)
        return df_res






class phylopAnnotator:
    def __init__(self,phylop_path):
        """
        specify the path to the phylop file
        """
        self.phylop = pyBigWig.open(phylop_path)
        
    def annotate(self, region):
        """
        region format: (chromosome, start, end)
        Output: a string of phylop scores separated by ":"
        """
        res=self.phylop.values(region[0], region[1], region[2],numpy=True).round(3)
        if res is None:
            return ""
        res=":".join(res.astype(str).tolist())
        return res





class baseImpAnnotator:
    def __init__(self,base_imp_path):
        self.base_imp = pd.read_csv(base_imp_path,index_col=0,header=None)
        
    def annotate(self,region,track_num,start_rel,end_rel):
        res=self.base_imp.loc[f"{region}_Track{track_num}"].round(3)
        res=res[start_rel:end_rel]
        res=":".join(res.astype(str).tolist())
        return res
    