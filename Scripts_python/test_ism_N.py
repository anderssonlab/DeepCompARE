import os
import unittest
import pandas as pd
import numpy as np
import shutil
from write_ism_N_seq import compute_ism_score, write_ism_score
from seq_ops import generate_random_seq, generate_random_seqs
from gradxinp import compute_gradxinp_from_seq


class TestISM_N(unittest.TestCase):

    @classmethod
    def setUpClass(cls):
        cls.wd = os.getcwd()
        os.makedirs(os.path.join(cls.wd, "Temp_test_ism"), exist_ok=True)

    def test_compute_ism_score_single_seq(self):
        seqs = generate_random_seq(600)
        df_ism = compute_ism_score(seqs)
        self.assertTrue(df_ism.shape == (16,600))
        self.assertTrue(np.all(df_ism.index==[f'Seq0_Track{i}' for i in range(16)]))
        
    
    
    def test_compute_ism_score_multi_seqs(self):
        seqs = generate_random_seqs(3, 600)
        # method 1: iterate over seqs
        df_ism1=pd.DataFrame()
        for seq in seqs:
            df_temp = compute_ism_score(seq)
            df_ism1 = pd.concat([df_ism1, df_temp], axis=0)
            
        # method 2: bulk mode
        df_ism2=compute_ism_score(seqs)
        self.assertTrue(np.allclose(df_ism1.values, df_ism2.values, atol=1e-5))
        self.assertTrue(np.all(df_ism2.index==[f'Seq{i}_Track{j}' for i in range(3) for j in range(16)]))
        
    def test_write_ism_score(self):
        seqs = generate_random_seqs(3, 600)
        pd.DataFrame(seqs, columns=["sequence"]).to_csv(os.path.join(self.wd, "Temp_test_ism", "seqs_ref.csv"), index=False)    
        
        # compute
        df_ism1 = compute_ism_score(seqs)        
        
        # write
        write_ism_score(os.path.join(self.wd, "Temp_test_ism", "seqs_ref.csv"),
                        "sequence",
                        os.path.join(self.wd, "Temp_test_ism", "ism_multi_seqs.csv"))
        df_ism2 = pd.read_csv(os.path.join(self.wd, "Temp_test_ism", "ism_multi_seqs.csv"), index_col=0, header=None) # no header!
        self.assertTrue(np.allclose(df_ism1.values, df_ism2.values, atol=1e-6))
        self.assertTrue(np.all(df_ism2.index==df_ism1.index))
        
    def test_corr_with_gradxinput(self):
        seqs=generate_random_seqs(3, 600)
        # compute gradxinput
        gradxinp=compute_gradxinp_from_seq(seqs)
        # compute ism
        df_ism=compute_ism_score(seqs)
        # compute correlation
        corr=np.corrcoef(gradxinp.values.flatten(),df_ism.values.flatten())[0,1]
        self.assertTrue(corr>0.1)
        
    def test_idx(self):
        seqs=generate_random_seqs(100, 600)
        df_ism=compute_ism_score(seqs)
        # assure index is ordered numerically
        self.assertTrue(np.all(df_ism.index==[f'Seq{i}_Track{j}' for i in range(100) for j in range(16)]))
    
    @classmethod
    def tearDownClass(cls):
        shutil.rmtree(os.path.join(cls.wd, "Temp_test_ism"))
        pass



if __name__ == '__main__':
    unittest.main()