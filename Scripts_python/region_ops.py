from pybedtools import BedTool



def subset_df_by_region(df,region,by):
    """
    Args:
        df: A data frame of bed file, first 3 columns should be chrom,start,end
        region: A tuple of (chromosome, start, end)
    Returns:
        if by=="1bp": return data frame of bed file with only the rows that overlap (at least 1bp) with the region 
        if by=="contained": return the rows that are completely contained in the region
        if by=="reverse": return the rows that do not overlap with the region
    """
    if by=="1bp":
        return df[(df.iloc[:,0]==region[0]) & (df.iloc[:,2]>=region[1]) & (df.iloc[:,1]<=region[2])].copy().reset_index(drop=True)
    if by=="contained":
        return df[(df.iloc[:,0]==region[0]) & (df.iloc[:,1]>=region[1]) & (df.iloc[:,2]<=region[2])].copy().reset_index(drop=True)
    if by=="reverse":
        return df[(df.iloc[:,0]!=region[0]) | (df.iloc[:,2]<region[1]) | (df.iloc[:,1]>region[2])].copy().reset_index(drop=True)



def resize_region(region,width,fix="center"):
    """
    Args:
        region: A tuple of (chromosome, start, end)
        width: An integer of the new width
        fix: A string of "center" or "left"
    Returns:
        A tuple of (chromosome, start, end)
    """
    if fix=="center":
        center=(region[1]+region[2])//2
        return (region[0],center-width//2, center+width//2)
    if fix=="left":
        return (region[0],region[1],region[1]+width)
    if fix=="expand":
        return (region[0],region[1]-width,region[2]+width)  
    
    
def resize_df(df_region,width,fix="center"):
    """
    Args:
        df_region: A data frame containing columns "start", "end" (have to be this name)
        width: An integer of the new width
        fix: A string of "center" or "left"
    Returns:
        A data frame, and the column "start" and "end" are modified
    """
    df_region=df_region.copy()
    if fix=="center":
        centers=(df_region["start"]+df_region["end"])//2
        df_region["start"]=centers-width//2
        df_region["end"]=centers+width//2
        return df_region
    if fix=="left":
        df_region["end"]=df_region["start"]+width
        return df_region
    if fix=="expand":
        df_region["start"]=df_region["start"]-width
        df_region["end"]=df_region["end"]+width
        return df_region
    
    
    
def calc_gc_context(df_region,context_width,seq_extractor):
    """
    Given a region, return the GC content of the context
    """
    df_resized=resize_df(df_region,context_width,fix="expand")
    df_resized["seq"]=df_resized.apply(lambda row: seq_extractor.get_seq(row["chromosome"], row["start"], row["end"]), axis=1)
    df_resized["gc"]=df_resized["seq"].apply(lambda x: (x.count("G")+x.count("C"))/len(x))
    return df_resized["gc"]


# This is probably problematic
def merge_intervals(df, 
                    other_cols=['protein'],    #['protein', 'max_gradxinp','mean_gradxinp','mean_abs_gradxinp'],
                    operations=["mode"],  #["mode","mean","mean","mean"]
                    ):
    bed = BedTool.from_dataframe(df)
    col_idxs = [df.columns.get_loc(col) + 1 for col in other_cols]  # +1 because BedTool columns are 1-indexed
    col_str = ','.join(map(str, col_idxs))
    op_str = ','.join(operations)
    merged = bed.merge(c=col_str, o=op_str).to_dataframe(names=['chromosome', 'start', 'end'] + other_cols)
    return merged



