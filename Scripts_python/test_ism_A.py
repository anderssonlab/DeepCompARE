import os
import shutil
import unittest

import numpy as np
import pandas as pd

from write_ism_A_seq import mutate_all_seqs, compute_ism_score, write_ism_score

###########
import sys
sys.path.insert(1, "/isdata/alab/people/pcr980/DeepCompare/Scripts_python/")

from seq_ops import generate_random_seq, generate_random_seqs
from gradxinp import compute_gradxinp_from_seq
###########


###############################################################################

class TestMutateAllSeqs(unittest.TestCase):

    def test_SingleSeq_colnames(self):
        df = mutate_all_seqs(['AAA'])
        correct_cols = [f'mutated_sequence{i}' for i in range(1, 5)]
        self.assertEqual(df.columns.tolist(), correct_cols)


    def test_SingleSeq_index(self):
        df = mutate_all_seqs(['AAA'])
        correct_idx = [f'Seq0_{i}' for i in range(3)]
        self.assertEqual(df.index.tolist(), correct_idx)


    def test_SingleSeq_shape(self):
        df = mutate_all_seqs(['AAA'])
        self.assertEqual(df.shape, (3,4))
    

    def test_SingleSeq_row(self):
        original_seq = "A"*10 + "C"*10 + "T"*10 + "G"*10 + "N"*10
        loc = np.random.choice(np.arange(len(original_seq)))
        row = mutate_all_seqs([original_seq]).iloc[loc, :].values

        loc_vals = {original_seq[loc]}
        seq_without_loc = original_seq[:loc] + original_seq[loc+1:]
        for elem in row:
            row_without_elem = elem[:loc] + elem[loc+1:]
            self.assertEqual(row_without_elem, seq_without_loc)
            
            loc_vals.add(elem[loc])
        
        self.assertEqual(loc_vals, {'A', 'T', 'C', 'G', 'N'})
    

    def test_MultipleSeq_index(self):
        df = mutate_all_seqs(['A'*20, 'C'*20, 'T'*20, 'G'*20, 'N'*20, 
                              'ACGTN'*4])
        correct_idx = [f'Seq{i}_{j}' for i in range(6) for j in range(20)]
        self.assertEqual(df.index.tolist(), correct_idx)


    def test_MultipleSeq_shape(self):
        df = mutate_all_seqs(['A'*20, 'C'*20, 'T'*20, 'G'*20, 'N'*20, 
                              'ACGTN'*4])
        self.assertEqual(df.shape, (120,4))
    

    def test_MultipleSeq_row(self):
        seqs = ['A'*20, 'C'*20, 'T'*20, 'G'*20, 'N'*20, 'ACGTN'*4]
        chosen_seq = np.random.choice(np.arange(6))
        loc = np.random.choice(np.arange(20))
        row = mutate_all_seqs(seqs).loc[f"Seq{chosen_seq}_{loc}", :].values

        original_seq = seqs[chosen_seq]
        loc_vals = {original_seq[loc]}
        seq_without_loc = original_seq[:loc] + original_seq[loc+1:]
        for elem in row:
            row_without_elem = elem[:loc] + elem[loc+1:]
            self.assertEqual(row_without_elem, seq_without_loc)
            
            loc_vals.add(elem[loc])
        
        self.assertEqual(loc_vals, {'A', 'T', 'C', 'G', 'N'})
        


class TestISM_A(unittest.TestCase):

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


###############################################################################

if __name__ == '__main__':
    unittest.main()
