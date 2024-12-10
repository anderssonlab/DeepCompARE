import pandas as pd
import os
import pyBigWig
import sys

sys.path.insert(1,"/isdata/alab/people/pcr980/Scripts_python")
from region_ops import subset_df_by_region




def find_binding_evidence(motif_df, tf_list, new_col_name):
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
        protein1 = motif_df['protein'].str.split("::", expand=True).iloc[:, 0]
        print(protein1)
        protein2 = motif_df['protein'].str.split("::", expand=True).iloc[:, -1]
        protein2 = protein2.fillna(protein1)
        print(protein2)
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
        motif_df=find_binding_evidence(motif_df, self.df_rna['TF'].values, 'rna_evidence')
        return motif_df
        









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
        # Retrieve entries from the JASPAR file
        start, end = region[1], region[2]
        tf_info_list = self.jaspar.entries(region[0], start, end)
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
        # if region is a tuple or list, it's only one region, call annotate_region
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
        #
        if isinstance(regions, pd.DataFrame) and outpath is not None:
            for idx, region in regions.iterrows():
                df = self._annotate_region((region[0], region[1], region[2]))
                df["seq_idx"]=f"Seq{idx}"
                df=self._add_or_subset_binding_evidence(df)
                if not os.path.isfile(outpath):
                    df.to_csv(outpath,index=False,mode="w",header=True)
                else:
                    df.to_csv(outpath,index=False,mode="a",header=False)




