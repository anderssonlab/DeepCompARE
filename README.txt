All R code should run on KU server. All python code should run on styx server.

Raw_data/ is downloaded from ENCODE or produced in Andersson lab.

Filter_regions/ stores the code to clean up raw data and produce Pd1_bed_processed/

Quantify_regions/ stores the code to quantify regions and produce Datasets/. The python code should run on styx.

Train/ stores the code for training models.

Test/ stores the code for benchmarking models.

Naming convention:
    write_* : write to csv format, return NULL.
    compute_: returns numpy array
    feat_attr: direct output from captum (e.g. gradxinp)
    feat_imp: including feat_attr, ism_delta and importance derived from all other methods



