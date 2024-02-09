import pandas as pd
from region_ops import *

region=("chr1",300,400)
assert resize_region(region,600)==("chr1",50,650)

region=("chr1",300,309)
assert resize_region(region,600)==("chr1",4,604)

df=pd.DataFrame([["chr1",300,400],["chr10",11300,14000]],columns=["chromosome","start","end"])
df_resized=resize_df(df,600)
df_resized_answer=pd.DataFrame([["chr1",50,650],["chr10",12350,12950]],columns=["chromosome","start","end"])
assert df_resized.equals(df_resized_answer)


def test_merge_intervals():
    df=pd.DataFrame([["chr1",1,10,"a"],["chr1",5,15,"a"],["chr1",20,30,"a"]],columns=['chrom','start','end','protein'])
    merged=merge_intervals(df)
    assert merged.equals(pd.DataFrame([["chr1",1,15,"a"],["chr1",20,30,"a"]],columns=['chrom','start','end','protein']))
