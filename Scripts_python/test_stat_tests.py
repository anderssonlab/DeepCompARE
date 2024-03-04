import pandas as pd
import unittest

from stat_tests import *

class TestBinAndLabel(unittest.TestCase):
    def test_bin_and_label(self):
        # Prepare test data and parameters
        data = {'Value': [0.01, 0.03, 0.07, 0.11, 0.15, 0.22, 0.25]}
        df = pd.DataFrame(data)
        column_name = 'Value'
        bin_edges = [0, 0.05, 0.1, 0.2, 0.3]
        expected_labels = ['0-0.05', '0-0.05', '0.05-0.1', '0.1-0.2', '0.1-0.2', '0.2-0.3', '0.2-0.3']
        expected_bin_column = pd.Series(pd.Categorical(expected_labels, categories=['0-0.05', '0.05-0.1', '0.1-0.2', '0.2-0.3'],
                                                       ordered=True),
                                        name='Bin')
        # Execute the function under test
        modified_df = bin_and_label(df, column_name, bin_edges)
    
        # Assert the outcomes
        pd.testing.assert_series_equal(modified_df['Bin'], expected_bin_column, check_dtype=True)

if __name__ == "__main__":
    unittest.main()