# Concatenate the two files and remove duplicates
cat expressed_tf_list_Hep-G2.tsv expressed_tf_list_K-562.tsv | sort | uniq > "expressed_tf_list_Hep-G2|K-562.tsv"
