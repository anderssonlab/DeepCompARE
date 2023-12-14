from Bio import SeqIO
import sys
sys.path.insert(1,"/isdata/alab/people/pcr980/DeepCompare/Scripts_python/")
from ism import *

# Define genomic regions
regions = [
    ("chr10", 132396935, 132397534),
    ("chr10", 8054655, 8055254),
    ("chr15", 37098712, 37099311),
    ("chr12", 132828699, 132829298),
    ("chr13", 114234732, 114235331),
    ("chr17", 36544867, 36545466)
]


def get_sequence(chrom, start, end):
    sequence = hg38_genome[chrom].seq[start-1:end]
    return str(sequence).upper()

with open("/isdata/alab/people/pcr980/Hg38/hg38.fa", "r") as file:
    hg38_genome = SeqIO.to_dict(SeqIO.parse(file,"fasta"))
    seqs=[get_sequence(*region) for region in regions]
    

    


