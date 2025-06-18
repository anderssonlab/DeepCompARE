import unittest
import pandas as pd
from seq_annotators import JasparAnnotator, BindingEvidenceAnnotator

class TestAnnotators(unittest.TestCase):

    def setUp(self):
        self.region = ("chr1", 109274652, 109275251)
        self.jaspar_file = "/isdata/alab/people/pcr980/Resource/JASPAR2022_tracks/JASPAR2022_hg38.bb"
        self.chip_file = "/isdata/alab/people/pcr980/Resource/ReMap2022/ReMap2022_hg38_hepg2.bed"
        self.rna_file = "/isdata/alab/people/pcr980/DeepCompare/RNA_expression/expressed_tf_list_hepg2.tsv"
        jaspar_annotator1 = JasparAnnotator(self.jaspar_file, "contained",
                                            chip_file=self.chip_file,
                                            rna_file=self.rna_file)
        self.df_motif1 = jaspar_annotator1.annotate(self.region)


    def test_jaspar_annotator_with_chip_and_rna(self):
        # Test JasparAnnotator with both chip and rna files
        unique_rna_evidence_values = self.df_motif1[self.df_motif1["chip_evidence"] == True]["rna_evidence"].unique()
        self.assertEqual(len(unique_rna_evidence_values), 1)


    def test_jaspar_annotator_with_chip_only(self):
        # Test JasparAnnotator with chip file only
        jaspar_annotator2 = JasparAnnotator(self.jaspar_file, "contained", chip_file=self.chip_file)
        df_motif2 = jaspar_annotator2.annotate(self.region)
        self.assertTrue("rna_evidence" not in df_motif2.columns)
        pd.testing.assert_frame_equal(df_motif2, self.df_motif1.drop(columns=["rna_evidence"]), check_dtype=False)


    def test_jaspar_annotator_with_rna_only(self):
        # Test JasparAnnotator with rna file only
        jaspar_annotator3 = JasparAnnotator(self.jaspar_file, "contained", rna_file=self.rna_file)
        df_motif3 = jaspar_annotator3.annotate(self.region)
        self.assertTrue("chip_evidence" not in df_motif3.columns)
        pd.testing.assert_frame_equal(df_motif3, self.df_motif1.drop(columns=["chip_evidence"]), check_dtype=False)


    def test_jaspar_annotator_without_chip_or_rna(self):
        # Test JasparAnnotator without chip or rna files
        jaspar_annotator4 = JasparAnnotator(self.jaspar_file, "contained")
        df_motif4 = jaspar_annotator4.annotate(self.region)
        self.assertTrue("chip_evidence" not in df_motif4.columns)
        self.assertTrue("rna_evidence" not in df_motif4.columns)
        pd.testing.assert_frame_equal(df_motif4, self.df_motif1.drop(columns=["chip_evidence", "rna_evidence"]), check_dtype=False)


    def test_binding_evidence_rna(self):
        # Test BindingEvidenceAnnotator with rna mode on df_chip
        df_motif=self.df_motif1.drop(columns=["rna_evidence", "chip_evidence"]).iloc[:10, :]
        be_annotator1 = BindingEvidenceAnnotator(self.chip_file, "chip")
        df_chip = be_annotator1.annotate(df_motif)
        be_annotator2 = BindingEvidenceAnnotator(self.rna_file, "rna")
        df_rna = be_annotator2.annotate(df_chip)
        pd.testing.assert_frame_equal(df_rna, self.df_motif1.iloc[:10, :], check_dtype=False)

if __name__ == '__main__':
    unittest.main()

