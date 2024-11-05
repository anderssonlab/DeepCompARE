# Concatenate the two files and remove duplicates
cat expressed_tf_list_hepg2.tsv expressed_tf_list_k562.tsv | sort | uniq > "expressed_tf_list_hepg2|k562.tsv"
