Using bigwig files in Raw_data, quantify the selected regions and output dataset for training and testing

preprocess_bw.R: 
    For DHS, SuRE, sum the signal on plus and minus strand, then average the signal per replicate
    For STARR, since the bigwig file downloaded from ENCODE already reflects replicate-merged bam, there is no need to do any preprocessing, simply copy paste the original data.
    For CAGE, go through CAGE_analysis process first, output.Rdata
    
    
Quant_signal_processed/: stores the processed .bw and .Rdata files

generate_csv_files.R: Generate 5-fold cross validation data, loose for regression data, as well as final-testing data, in csv form.

CV5_data_chunks/: intermediate data chunks to generate 5-fold cross validation. Binary labels are added. Loose regions are removed.

CV5_loose_chunks/: intermediate data chunks to generate 5-fold cross validation. No binary labels. Loose regions are maintained.


To convert csv files to npy for training, call csv2npy.py through csv2npy_*.sh using styx server

