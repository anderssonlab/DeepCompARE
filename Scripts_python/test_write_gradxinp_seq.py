import os
import pandas as pd
import numpy as np
import shutil
import unittest

from seq_ops import generate_random_seqs
from write_gradxinp_seq import write_gradxinp_from_seq_file


class TestGradXinp(unittest.TestCase):

    def setUp(self):
        self.NUM_TRACKS = 16
        self.NUM_SEQS = 1
        self.SEQ_LEN = 600
        self.wd = os.getcwd()
        os.makedirs(os.path.join(self.wd, "Temp_test_graxinp"), exist_ok=True)
        self.seqs = generate_random_seqs(self.NUM_SEQS, self.SEQ_LEN)
        df_seq = pd.DataFrame(self.seqs, columns=["seq"])
        df_seq.to_csv(os.path.join(self.wd, "Temp_test_graxinp", "seqs.csv"), index=False)

    def tearDown(self):
        shutil.rmtree(os.path.join(self.wd, "Temp_test_graxinp"))

    def check_idx(self, imp_df, seqs):
        idx_components = pd.Series(imp_df.index).str.split("_", expand=True)
        self.assertTrue(np.all(idx_components.iloc[:, 0] == [f"Seq{number}" for number in range(len(seqs)) for _ in range(self.NUM_TRACKS)]))
        self.assertTrue(np.all(idx_components.iloc[:, 1] == [f"Track{number}" for number in list(range(self.NUM_TRACKS)) * len(seqs)]))
        self.assertTrue(imp_df.index.is_unique)

    def test_write_gradxinp_from_seq_file(self):
        write_gradxinp_from_seq_file(os.path.join(self.wd, "Temp_test_graxinp", "seqs.csv"), 
                                "seq", 
                                os.path.join(self.wd, "Temp_test_graxinp", "gradxinp_bulk.csv"), 
                                batch_size=4)
        imp_df_bulk = pd.read_csv(os.path.join(self.wd, "Temp_test_graxinp", "gradxinp_bulk.csv"), index_col=0, header=None)
        self.check_idx(imp_df_bulk, self.seqs)


if __name__ == '__main__':
    unittest.main()