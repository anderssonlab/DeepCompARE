This directory subset Jaspar database by transcription factors expressed in HegP2 and K562 cells.
Then add ChIP-Seq info to distinguish putative binding from actual binding.

Raw_data/:
Total RNA-Seq of HepG2 is downloaded from https://www.encodeproject.org/experiments/ENCSR181ZGR/
Total RNA-Seq of K562 is downloaded from https://www.encodeproject.org/experiments/ENCSR792OIJ/

curate.R:
Input gene expression data, output lists of expressed and unexpressed TF, for HepG2 and K562 cell lines.