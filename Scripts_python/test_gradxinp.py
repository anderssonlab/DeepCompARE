import os
import pandas as pd
from pandas.testing import assert_frame_equal
import numpy as np
import shutil
import unittest
import sys
sys.path.insert(1, "/isdata/alab/people/pcr980/DeepCompare/Scripts_python")
from gradxinp import *



class TestGradXinp(unittest.TestCase):

    def setUp(self):
        self.NUM_TRACKS = 16
        self.wd = os.getcwd()
        os.makedirs(os.path.join(self.wd, "Temp_test_graxinp"), exist_ok=True)
        self.bed_file_path = os.path.join(self.wd, "Temp_test_graxinp", "test.bed")
        self.bed_file = pd.DataFrame({"chr": ["chr1", "chr2", "chr3"], 
                                      "start": [1053, 20000, 300001]})
        self.bed_file["end"] = self.bed_file["start"] + 599
        self.bed_file.to_csv(self.bed_file_path, sep="\t", header=False, index=False)
        self.df_gradxinp=compute_gradxinp(self.bed_file_path)
        self.df_gradxinp=compute_gradxinp(self.bed_file)
                
    def tearDown(self):
        shutil.rmtree(os.path.join(self.wd, "Temp_test_graxinp"))

    def check_shape(self):
        self.assertEqual(self.df_gradxinp.shape, (self.NUM_TRACKS * len(self.bed_file), 600))

        
    def test_compute_gradxinp_indices(self):
        indices = [f"{self.bed_file.iloc[i,0]}:{self.bed_file.iloc[i,1]}-{self.bed_file.iloc[i,2]}_Track{target}" for i in range(len(self.bed_file)) for target in range(self.NUM_TRACKS)]
        self.assertListEqual(self.df_gradxinp.index.tolist(), indices)

    def test_compute_gradxinp_targets(self):
        # iter over rows
        for i in range(len(self.bed_file)):
            for target in [1, 3, 5, 7, 9]:
                bed_file_subset = self.bed_file.iloc[i:i+1,:]
                df_target = compute_gradxinp(bed_file_subset, targets=target)
                assert_frame_equal(df_target, self.df_gradxinp.loc[df_target.index,:])
if __name__ == '__main__':
    unittest.main()