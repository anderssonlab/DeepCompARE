Preprocess each data modality, including preprocessing and liftover Bigwig files of SuRE-Seq. Then output reproducible regions to Pd1_bed_processed/


----------------------
Preprocess: SuRE 
----------------------
code: start with SuRE_
Raw_data/SuRE
  * SuRE_SNP_collection.R --> Filter_regions/SuRE/SuRE_SNP/
  * SuRE_process_bw.R --> Filter_regions/SuRE/BW_hg19/ --> liftover_sure_hg19.sh --> Filter_regions/SuRE/BW_hg38/

Filter_regions/SuRE/BW_hg19/
  * MACS3 bdgpeakcall -c 3 -l 300 -g 150 (Performed on local computer)
  * -c has choices 2,3,4
  
Filter_regions/SuRE/NarrowPeak_hg19
  * SuRE_narrowPeak_EDA_filter_liftover.R --> Raw_data/SuRE_files/NarrowPeak_hg38_c4_RMSNP
 
 lift hg19 to hg38 using listover_sure_hg19.sh
 
  
-----------------------
Preprocess: CAGE
----------------------
code: CAGE_analysis.R 

Irreproducible: expr_thresh=6, #min_support=1
Filter_regions/CAGE/Peaks_irreproducible

Reproducible within modality: expr_thresh=6, #min_support=2
Filter_regions/CAGE/Peaks_reproducible_within_modality/

Loose: expr_thresh=1, #min_overlap=1
Pd1_bed_processed/





-----------------------
Preprocess: STARR
----------------------
Only for loose
code: STARR_loose_preprocessing.R




--------------------
Reproducibility DHS 
---------------------
code: generate_reproducible_regions.R

Irreproducible: directly downloaded from ENCODE (summit already centered at 75)
Filter_regions/DHS/Peaks_irreproducible

Reproducible within modality: #overlap=2, min_length=120
Filter_regions/DHS/Peaks_reproducible_within_modality/

Loose: None



----------------------
Reproducibility STARR
----------------------
code: generate_reproducible_regions.R


Irreproducible: directly downloaded from ENCODE (no recenter)
/maps/projects/ralab/people/pcr980/Raw_data/STARR/Bed, also at Filter_regions/STARR/Peaks_irreproducible

Reproducible within modality: #overlap=2, min_length=400
Filter_regions/STARR/Peaks_reproducible_within_modality/


Loose: on replicate-merged bam file, run starrpeaker, q value threshold=0.2
Filter_regions/STARR/STARRPeaker_res


----------------------
Reproducibility SuRE
----------------------
code: generate_reproducible_regions.R

Irreproducible: -c4 RMSNP


Reproducible within modality: -c4, #overlap=2, min_length=240
Filter_regions/SuRE/Peaks_reproducible_within_modality/

Loose: -c2



