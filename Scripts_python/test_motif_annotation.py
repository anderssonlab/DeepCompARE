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
        self.df = pd.DataFrame.from_dict(
            {
                "chromosome": ["chr1", "chr1", "chr1"],
                "start": [11100, 11200, 11300],
                "end": [11110, 11215, 11309]
            }
        )
        self.region=("chr1", 11000, 11500)
        self.sequence=self.seq_extractor.get_seq(self.region[0],self.region[1],self.region[2])

    def test_sequence_extraction(self):
        df_res=add_feat_imp(self.df,self.region,np.random.rand(500))
        df_res["motif_sequence1"] = df_res.apply(lambda row: self.seq_extractor.get_seq(row["chromosome"], row["start"], row["end"]), axis=1)
        df_res["motif_sequence2"] = df_res.apply(lambda row: self.sequence[row["start_rel"]:row["end_rel"]+1], axis=1)
        self.assertEqual((df_res["motif_sequence1"]==df_res["motif_sequence2"]).all(),True,
                          "The sequences extracted from absolute coord and relative coord should be the same")


if __name__ == '__main__':
    unittest.main()
