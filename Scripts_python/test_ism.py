import unittest
import os
import pandas as pd
import shutil
from pandas.testing import assert_frame_equal
from in_silico_mutagenesis import compute_mutagenesis_score
from gradxinp import compute_gradxinp


class TestMutagenesisScore(unittest.TestCase):

    def setUp(self):
        self.NUM_TRACKS = 16
        self.wd = os.getcwd()
        os.makedirs(os.path.join(self.wd, "Temp_test_mutagenesis"), exist_ok=True)
        self.bed_file_path = os.path.join(self.wd, "Temp_test_mutagenesis", "test.bed")
        self.bed_file = pd.DataFrame({"chr": ["chr1", "chr2", "chr3"], 
                                      "start": [1053, 20000, 300001]})
        self.bed_file["end"] = self.bed_file["start"] + 599
        self.bed_file.to_csv(self.bed_file_path, sep="\t", header=False, index=False)
        
        self.df_isa = compute_mutagenesis_score(
            self.bed_file_path,
            mode="isa",
            aggregation_method="mean"
        )
        
        self.df_ism = compute_mutagenesis_score(
            self.bed_file_path,
            mode="ism",
            aggregation_method="mean"
        )
        
    def tearDown(self):
        shutil.rmtree(os.path.join(self.wd, "Temp_test_mutagenesis"))
        
    def test_compute_mutagenesis_shape(self):
        expected_shape = (self.NUM_TRACKS * len(self.bed_file), 600)
        self.assertEqual(self.df_isa.shape, expected_shape)

    def test_compute_mutagenesis_indices(self):
        indices = [
            f"{self.bed_file.iloc[i, 0]}:{self.bed_file.iloc[i, 1]}-{self.bed_file.iloc[i, 2]}_Track{target}"
            for i in range(len(self.bed_file)) for target in range(self.NUM_TRACKS)
        ]
        self.assertListEqual(self.df_isa.index.tolist(), indices)

    def test_compute_mutagenesis_targets(self):
        # Iterate over rows and validate targets
        for i in range(len(self.bed_file)):
            bed_file_subset = self.bed_file.iloc[i:i+1, :]
            df_target = compute_mutagenesis_score(
                bed_file_subset,
                mode="isa",
                aggregation_method="mean"
            )
            assert_frame_equal(
                df_target,
                self.df_isa.loc[df_target.index, :]
            )
    
    def test_corr_ism_isa(self):
        ism_flat = self.df_ism.values.flatten()
        isa_flat = self.df_isa.values.flatten()
        corr = pd.Series(ism_flat).corr(pd.Series(isa_flat))
        self.assertGreater(corr, 0.5)
    
    def test_corr_isa_gradxinp(self):
        isa_flat = self.df_isa.values.flatten()
        gradxinp = compute_gradxinp(self.bed_file_path)
        gradxinp_flat = gradxinp.values.flatten()
        corr = pd.Series(isa_flat).corr(pd.Series(gradxinp_flat))
        self.assertGreater(corr, 0.5)
    

###############################################################################

if __name__ == '__main__':
    unittest.main()
