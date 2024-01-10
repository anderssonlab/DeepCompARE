import unittest
import pandas as pd
import random
import sys
sys.path.insert(1,"/isdata/alab/people/pcr980/DeepCompare/Scripts_python")
from motif_annotation import *
from utils_for_test import generate_random_seqs
from utils import SeqExtractor

class TestSeqExtraction(unittest.TestCase):

    def setUp(self):
        self.seq_extractor = SeqExtractor()

    def test_sequence_extraction(self):
        df = pd.DataFrame.from_dict(
            {
                "chromosome": ["chr1", "chr1", "chr1"],
                "start_motif": [11100, 111200, 300],
                "end_motif": [11110, 111215, 309],
                "start_seq": [11050, 111150, 250],
                "end_seq": [11150, 111250, 350],
            }
        )
        df["sequence"] = df.apply(lambda row: self.seq_extractor.get_seq(row["chromosome"], row["start_seq"], row["end_seq"]), axis=1)
        for i in range(3):
            start_rel, end_rel = compute_relative_location((df['chromosome'][i], df['start_motif'][i], df['end_motif'][i]),
                                                           (df['chromosome'][i], df['start_seq'][i], df['end_seq'][i]))
            motif_seq1 = self.seq_extractor.get_seq(df['chromosome'][i], df['start_motif'][i], df['end_motif'][i])
            motif_seq2 = df["sequence"][i][start_rel:end_rel+1]
            self.assertEqual(motif_seq1, motif_seq2, "Extracted sequences do not match")
    
    def test_add_feat_imp(self):
        df = pd.DataFrame.from_dict(
            {
                "chromosome": ["chr1", "chr1", "chr1"],
                "start": [11100, 11200, 11300],
                "end": [11110, 11215, 11309]
            }
        )
        region=("chr1", 11000, 11500)
        df=add_feat_imp(df,region,np.random.rand(500))
        self.assertEqual(df.shape[0],3,"The number of rows should not change")
        self.assertEqual(df.shape[1],7,"The number of columns should increase by 4")
            

if __name__ == '__main__':
    unittest.main()
