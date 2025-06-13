import pandas as pd

import sys
sys.path.insert(1,"/isdata/alab/people/pcr980/Scripts_python")
from region_ops import resize_df

def skip_comment_lines(file_path):
    # skip the rows starting with '#'
    with open(file_path, 'r') as file:
        # Read lines and determine which to skip
        lines = file.readlines()
        return [i for i, line in enumerate(lines) if line.startswith('#')]



#------------------------------------------------------------------------------------------------
# Promoter
#------------------------------------------------------------------------------------------------
# Specify the path to your file
file_path = '41586_2022_4877_MOESM2_ESM.txt'
# Identify rows to skip
rows_to_skip = skip_comment_lines(file_path)
# Load the file, skipping the identified rows
df = pd.read_csv(file_path, sep='\t', skiprows=rows_to_skip)
df=df[df["chr"].str.contains("chr")].reset_index(drop=True)
for i in [0,1,2]:
    df[df['promoter_exp_class']==i].reset_index(drop=True).to_csv(f'P{i}_k562.tsv', sep='\t', index=False)





#------------------------------------------------------------------------------------------------
# Enhancer
#------------------------------------------------------------------------------------------------
file_path = '41586_2022_4877_MOESM3_ESM.txt'
rows_to_skip = skip_comment_lines(file_path)
df = pd.read_csv(file_path, sep='\t', skiprows=rows_to_skip)
df=df[df["chr"].str.contains("chr")].reset_index(drop=True)
for i in [0,1,2]:
    df[df['enhancer_exp_class']==i].reset_index(drop=True).to_csv(f'E{i}_k562.tsv', sep='\t', index=False)

