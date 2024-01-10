import numpy as np
import pandas as pd
import os
import unittest
from utils_for_test import generate_random_seqs
from prediction import compute_predictions, write_predictions


class TestPredictionMethods(unittest.TestCase):
    def setUp(self):
        self.wd = os.getcwd()
        self.test_dir = os.path.join(self.wd, "Temp_test_prediction")
        os.makedirs(self.test_dir, exist_ok=True)
        self.seqs = generate_random_seqs(20000, 600)
        self.res1 = compute_predictions(self.seqs)
        
    def test_dimensions(self):
        # Check if the output has correct dimensions
        self.assertEqual(self.res1.shape, (20000, 16))
        
    def test_prediction_results(self):
        # Save sequences to CSV
        pd.DataFrame(self.seqs, columns=["sequence"]).to_csv(os.path.join(self.test_dir, "seqs.csv"), index=False)

        # Write predictions
        write_predictions(os.path.join(self.test_dir, "seqs.csv"),
                          "sequence",
                          os.path.join(self.test_dir, "preds.csv"))

        # Read predictions from CSV
        res2 = pd.read_csv(os.path.join(self.test_dir, "preds.csv"), header=None).values

        # Check if shapes are equal
        self.assertEqual(self.res1.shape, res2.shape)

        # Check if values are almost equal
        self.assertTrue(np.allclose(self.res1, res2, atol=1e-6))

    def tearDown(self):
        # Clean up the directory after tests
        for file in os.listdir(self.test_dir):
            os.remove(os.path.join(self.test_dir, file))
        os.rmdir(self.test_dir)

if __name__ == '__main__':
    unittest.main()
