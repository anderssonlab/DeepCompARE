import os
import unittest
import pandas as pd
import numpy as np
import shutil
from ism import calculate_ism_delta
from utils_for_test import generate_random_seq, generate_random_seqs



class TestISM(unittest.TestCase):

    @classmethod
    def setUpClass(cls):
        cls.wd = os.getcwd()
        os.makedirs(os.path.join(cls.wd, "Temp_test_ism"), exist_ok=True)

    def test_single_seq(self):
        # generate sequence
        seqs = generate_random_seq(600)
        pd.DataFrame([seqs], columns=["sequence"]).to_csv(os.path.join(self.wd, "Temp_test_ism", "seq_ref.csv"), index=False)
        
        # single mode
        df_ism = calculate_ism_delta(seqs)
        df_ism.to_csv(os.path.join(self.wd, "Temp_test_ism", "ism_single_seq_single_mode.csv"), header=False, index=False)
        
        # bulk mode
        calculate_ism_delta(os.path.join(self.wd, "Temp_test_ism", "seq_ref.csv"),
                            "sequence",
                            os.path.join(self.wd, "Temp_test_ism", "ism_single_seq_bulk_mode.csv"))

        df_single_mode = pd.read_csv(os.path.join(self.wd, "Temp_test_ism", "ism_single_seq_single_mode.csv"), header=None)
        df_bulk_mode = pd.read_csv(os.path.join(self.wd, "Temp_test_ism", "ism_single_seq_bulk_mode.csv"))
        self.assertTrue(np.allclose(df_single_mode.values, df_bulk_mode.iloc[:, 2:].values, atol=1e-6))
        
    def test_multi_seqs(self):
        # generate sequence
        seqs = generate_random_seqs(10, 600)
        pd.DataFrame(seqs, columns=["sequence"]).to_csv(os.path.join(self.wd, "Temp_test_ism", "seqs_ref.csv"), index=False)
        
        # single mode
        for seq in seqs:
            df_ism = calculate_ism_delta(seq)
            df_ism.to_csv(os.path.join(self.wd, "Temp_test_ism", "ism_multi_seqs_single_mode.csv"), mode="a",header=False, index=False)
        
        # bulk mode
        calculate_ism_delta(os.path.join(self.wd, "Temp_test_ism", "seqs_ref.csv"),
                            "sequence",
                            os.path.join(self.wd, "Temp_test_ism", "ism_multi_seqs_bulk_mode.csv"))

        df_single_mode = pd.read_csv(os.path.join(self.wd, "Temp_test_ism", "ism_multi_seqs_single_mode.csv"), header=None)
        df_bulk_mode = pd.read_csv(os.path.join(self.wd, "Temp_test_ism", "ism_multi_seqs_bulk_mode.csv"))
        self.assertTrue(np.allclose(df_single_mode.values, df_bulk_mode.iloc[:, 2:].values, atol=1e-5))
    
    # To do: test ism has OK correlation with gradxinput
    
    
    @classmethod
    def tearDownClass(cls):
        shutil.rmtree(os.path.join(cls.wd, "Temp_test_ism"))
        pass



if __name__ == '__main__':
    unittest.main()