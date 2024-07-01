import pandas as pd
import pyranges as pr
import pyBigWig
from region_ops import subset_df_by_region




class JasparAnnotator:
    def __init__(self,jaspar_path="/binf-isilon/alab/people/pcr980/Resource/JASPAR2022_tracks/JASPAR2022_hg38.bb"):
        self.jaspar = pyBigWig.open(jaspar_path)
        assert self.jaspar.isBigBed()

    def annotate(self, region,by="contained"):
        """Get putative motif information from a genomic region.
        Args:
            region: A tuple of (chromosome, start, end).
        Returns:
            A DataFrame of start, end, TF, score, and strand.
        """
        tf_info_list = self.jaspar.entries(region[0], region[1], region[2])
        df = pd.DataFrame(tf_info_list, columns=['start', 'end', 'details'])
        df['chromosome'] = region[0]
        df=df.loc[:, ['chromosome', 'start', 'end', 'details']]
        df = subset_df_by_region(df, region, by=by)
        if df.shape[0]==0:
            return df
        df[['protein', 'score', 'strand']] = df['details'].str.split('\t', expand=True)
        df.drop(columns='details', inplace=True)
        df['protein'] = df['protein'].str.upper()
        df['score'] = df['score'].astype(int)
        df = df.drop_duplicates(subset=['chromosome', 'start', 'end', 'strand', 'protein']).reset_index(drop=True)
        return df





class ReMapAnnotator:
    def __init__(self, file_path):
        # Reading the remap object once during initialization
        self.remap = pd.read_csv(file_path, sep='\t', header=None)
        self.remap.columns = ['chromosome', 'start', 'end', 'detail', 'x', 'strand', 
                              'thick_start', 'thick_end', 'score']
        self.remap_gr=pr.PyRanges(self.remap.rename(columns={"chromosome":"Chromosome","start":"Start","end":"End"}))
    
    def clean_df(self, df):
        # Extract TF from "detail" column
        df['TF'] = df['detail'].str.split(':', expand=True).iloc[:,0].str.upper()
        df.drop(columns=['detail', 'x', 'thick_start', 'thick_end', 'score'], inplace=True)
        return df

    def find_chip_evidence(self,chip_df,region):
        """
        Args: 
            chip_df: A data frame of ChIP-Seq data
            region: A row of motif_df
        Returns:
            A boolean value indicating whether the region has ChIP-Seq evidence
        """
        if region["protein"].find("::")==-1:
            if region["protein"] in chip_df['TF'].values:
                return True
            else:
                return False
        else: # dimer: require only one of the TFs to have ChIP evidence
            tf1=region["protein"].split("::")[0]
            tf2=region["protein"].split("::")[1]
            if tf1 in chip_df['TF'].values or tf2 in chip_df['TF'].values:
                return True
            else:
                return False

    def annotate(self, motif_df, region):
        # Step 1: Subset remap to the region of interest
        chip_seq=subset_df_by_region(self.remap,region,by="1bp")
        # if no Chip-Seq data in the region, add a column of False
        if chip_seq.shape[0]==0:
            motif_df['chip_evidence'] = False
            return motif_df
        chip_seq = self.clean_df(chip_seq)
        motif_df['chip_evidence'] = motif_df.apply(lambda row: self.find_chip_evidence(chip_seq, row), axis=1)
        return motif_df
    
    
    
def read_maf_file(chrom_num):
    df_maf=pd.read_csv(f"/isdata/alab/people/pcr980/Resource/Gnomad_vcf/Pd1_MAF/chrom{chrom_num}.csv")
    df_maf=df_maf[df_maf["REF"].str.len()==1]
    df_maf=df_maf[df_maf["ALT"].str.len()==1]
    df_maf=df_maf.reset_index(drop=True)
    df_maf["End"]=df_maf["POS"]+1
    df_maf.columns=["Chromosome","Start","ID","REF","ALT","AF","End"]
    df_maf=df_maf[["Chromosome","Start","End","ID","REF","ALT","AF"]]
    return df_maf
    


class gnomADSNPAnnotator:
    def __init__(self):
        # Read MAF of all chromosomes into one dictionary
        self.maf_dict = {}
        for chrom_num in list(range(1, 23))+["X","Y"]:
            self.maf_dict["chr"+str(chrom_num)] = read_maf_file(chrom_num)

    def annotate(self, region):
        return subset_df_by_region(self.maf_dict[region[0]],
                                   region,
                                   by="contained")
        