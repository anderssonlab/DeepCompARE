setwd("~/Documents/DeepCompare_assistent/Enformer")

df <- read.csv("targets_human.txt",sep="\t", header=TRUE)

# What assays are available:"DNASE","ATAC","CHIP","CAGE" 
assays <- sapply(df$description, function(txt) sub("^([^:]+):.*", "\\1", txt))
unique(assays)


# DNAse HepG2: 28,92,235
idx <- which(grepl("HepG2",df$description) & grepl("DNASE",df$description))
df$description[idx]

# CAGE HepG2: 5110
idx <- which(grepl("HepG2",df$description) & grepl("CAGE",df$description))
df$description[idx]

# DNAse K562: 34,35,36,122,123,124,626
idx <- which(grepl("K562",df$description) & grepl("DNASE",df$description))
df$description[idx]

# CAGE K562: 4829 5112
idx <- which(grepl("K562",df$description) & grepl("CAGE",df$description))
df$description[idx]


