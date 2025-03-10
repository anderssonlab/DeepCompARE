##########################################################################################
# R script exporting PPI information from InWeb and STRING databases 
# as present in Genoppi R package for TF-TF pairs analyzed in DeepCompARE model

# Input:
# tf_pair_cooperativity_index_0.5_0.1_k562.csv
# tf_pair_cooperativity_index_0.5_0.1_hepg2.csv

# Output:
# <date>_s10_PublishedPPIandProtComplexes_0.5_0.1_k562.txt
# <date>_s10_PublishedPPIandProtComplexes_0.5_0.1_hepg2.txt

# Author: Petra Palenikova, Xuening He
# Last updated: 09 Jan 2025 (Petra), 07 Mar 2025 (Xuening)
##########################################################################################

### load libraries ####
library(genoppi) ## version from Github (lagelab/Genoppi@85de482)
library(tidyverse) ## version 2.0.0


### functions ####
pullInts <- function(protein, stringTh = 0, inWebAllProts, stringAllProts){
  pInts <- NULL
  if(protein %in% inWebAllProts){
    get_inweb_list(protein) %>% 
      filter(significant) %>% 
      pull(gene) -> pInts
  }
  if(protein %in% stringAllProts){
    get_string_list(protein, score = stringTh) -> pString
    if(!is.null(pString)){
      pString %>% 
        filter(significant) %>% 
        pull(gene) -> pString
      pInts <- union(pInts,pString)
    }
  }
  return(pInts) 
}

### data processing ####
#### read cooperativity index data for both cell lines #####
dfCoop_k562 <- read_delim("tf_pair_cooperativity_index_k562_pe.csv",
                     show_col_types = FALSE)
dfCoop_hepg2 <- read_delim("tf_pair_cooperativity_index_hepg2_pe.csv",
                          show_col_types = FALSE)

#### filter and select required columns ####
### remove protein1 == protein2 rows 
## - since the aim is to compere TF-TF interaction in PPI databases these
##   would be difficult to interpret
### combine two datasets

dfCoop_k562 %>% 
  filter(protein1 != protein2) %>% 
  select(protein1,protein2) %>% 
  mutate(cellLine = "k562") -> dfCoop_k562_filtered

dfCoop_hepg2 %>% 
  filter(protein1 != protein2) %>% 
  select(protein1,protein2) %>% 
  mutate(cellLine = "hepg2") -> dfCoop_hepg2_filtered

dfCoop_filtered <- rbind(dfCoop_k562_filtered,dfCoop_hepg2_filtered)

#### create a new column with protein pairs ####
## some TF-TF pairs are present in both orientations
## for this new column sort protein1 protein2 alphabetically
dfCoop_filtered$protPair <- unlist(
  lapply(
    lapply(seq(1,nrow(dfCoop_filtered)),
           function(x) sort(c(dfCoop_filtered$protein1[x],dfCoop_filtered$protein2[x]))
    ),
    function(x) paste(x, collapse = "_")
  )
)

#### extract all unique protein pairs ####
## to add information about published PPIs
dfCoop_filtered %>% 
  select(protPair) %>% 
  unique() -> dfCoop_allPairs

#### add TF-TF interaction information ####

## get all proteins reported in InWeb
get_inweb_list("AKAP11") %>% pull(gene) -> allInWeb
allInWeb <- c(allInWeb,"AKAP11")
## get all proteins reported in STRING
get_string_list("AKAP11") %>% pull(gene) -> allString
allString <- c(allString,"AKAP11")
stringThreshold <- 0
dfCoop_withInts <- NULL
for(i in seq(1, nrow(dfCoop_allPairs))){
  #for(i in seq(495, 505)){
  if(i %% 500 == 0){
    print(paste("Done with", i, "pairs."))
  }
  doInteract <- "No"
  twoProts <- unlist(str_split(dfCoop_allPairs$protPair[i],"_")) #get protein pair
  if(sum(str_detect(twoProts,"::"))>1){ ## handle when both proteins are heterodimers
    pOne <- twoProts[1] ## take first heterodimer
    pOne_split <- unlist(str_split(pOne,"::")) ## split it into separate proteins
    ## pull interactors for both of these proteins
    pIntOne <- union(pullInts(pOne_split[1],stringTh = stringThreshold, allInWeb,allString),
                     pullInts(pOne_split[2],stringTh = stringThreshold, allInWeb,allString)) 
    ### look for either of the two prots from second heterodimer among interactors
    if(any(unlist(str_split(twoProts[2],"::")) %in% pIntOne)){
      doInteract <- "Yes"
    }
  } else if (sum(str_detect(twoProts,"::")) == 1) { ## when one is heterodimer
    pOne <- twoProts[!str_detect(twoProts,"::")] ## pick protein that is NOT a heterodimer
    ## pull interactors of this protein
    pIntOne <- pullInts(pOne,stringTh = stringThreshold, allInWeb,allString) 
    ### look for either of the two prots from heterodimer among interactors
    if(any(unlist(str_split(twoProts[str_detect(twoProts,"::")],"::")) %in% pIntOne)){
      doInteract <- "Yes"
    }
  } else { # when neither is heterodimer
    pOne <- twoProts[1] ## pull interactors for the first protein
    pIntOne <- pullInts(pOne,stringTh = stringThreshold, allInWeb,allString) 
    ### look for the second protein among interactors
    if(twoProts[2] %in% pIntOne){
      doInteract <- "Yes"
    }
  }
  dfCoop_withInts<- rbind(dfCoop_withInts,
                          cbind(dfCoop_allPairs[i,],
                                tibble("Reported_PPI" = doInteract)))
}

nrow(dfCoop_hepg2_filtered) + nrow(dfCoop_k562_filtered)
nrow(dfCoop_filtered)

#### merge PPI info with protein pairs ####
dfCoop_withInts_cellModels <- inner_join(dfCoop_filtered, dfCoop_withInts,
                                         by="protPair")

### Add information whether protein pair interacts with cofactor complexes ####

### A) mediator complex
### B) SWI/SNF complex 
### C) RNA pol II 
### D) TFII B,D,E,F,H
### pull PPIs for these protein complexes
dfComplexes <- read_delim("CofactorComplexes.txt")
allComplexes <- unique(dfComplexes$Complex)
complexInts <- NULL
for(cVal in allComplexes){
  dfComplexes %>% 
    filter(Complex == cVal) %>% 
    pull(Gene) -> complexGenes
  allInts <- NULL
  for(pVal in complexGenes){
    allInts <- union(allInts,pullInts(pVal, stringTh = 0,
                                      inWebAllProts = allInWeb,
                                      stringAllProts = allString))
  }
  dfTemp <- tibble("Complex" = cVal,
             "Interactors" = allInts)
  dfTemp$IsPPI <- dfTemp$Interactors %in% allInts
  complexInts <- rbind(complexInts,
                       dfTemp)
}


dfCoop_withInts_Complexes <- NULL
for(i in seq(1, nrow(dfCoop_allPairs))){
#for(i in seq(495, 505)){
  if(i %% 500 == 0){
    print(paste("Done with", i, "pairs."))
  }
  twoProts <- unlist(str_split(dfCoop_allPairs$protPair[i],"_"))
  dfTemp <- NULL
  for(cVal in allComplexes){
    complexInts %>% 
      filter(Complex == cVal, 
             IsPPI) %>% 
      pull(Interactors) -> complexIntProts
    if(sum(str_detect(twoProts,"::"))>1){ ## handle when both proteins are heterodimers
      tfPresent <- as.numeric(any(unlist(str_split(twoProts[1],"::")) %in% complexIntProts))

      tfPresent <- tfPresent + 
        as.numeric(any(unlist(str_split(twoProts[2],"::")) %in% complexIntProts))
      
    } else if (sum(str_detect(twoProts,"::")) == 1) { ## when one is heterodimer
      pOne <- twoProts[!str_detect(twoProts,"::")] ## take protein that is NOT a heterodimer
      tfPresent <- sum(pOne %in% complexIntProts)

      tfPresent <- tfPresent + 
        as.numeric(any(unlist(str_split(twoProts[str_detect(twoProts,"::")],"::")) %in% complexIntProts))
      
    } else { # when neither is heterodimer
      tfPresent <- sum(twoProts %in% complexIntProts)
    }
    if(is.null(dfTemp)){
    dfTemp <- data_frame(protPair = dfCoop_allPairs$protPair[i],
                         complexInt = factor(tfPresent, levels = c(0,1,2)))
    names(dfTemp) <- c("protPair",cVal)
    } else {
      dfTemp <- cbind(dfTemp,
                      data_frame(complexInt = factor(tfPresent, levels = c(0,1,2))))
      names(dfTemp)[ncol(dfTemp)] <- cVal
      
    }
  }
  dfCoop_withInts_Complexes <- rbind(dfCoop_withInts_Complexes, 
                                     dfTemp)
}

#### merge PPI complexes info with protein pairs ####
dfCoop_withInts_Complexes_cellModels <- inner_join(dfCoop_withInts_cellModels, dfCoop_withInts_Complexes,
                                         by="protPair")
 
dfCoop_withInts_Complexes_cellModels %>% 
  filter(cellLine == "k562") %>% 
  select(-cellLine) -> dfCoop_withInts_Complexes_k562

dfCoop_withInts_Complexes_cellModels %>% 
  filter(cellLine == "hepg2") %>% 
  select(-cellLine) -> dfCoop_withInts_Complexes_hepg2

### export files with PPI annotation ####
write.table(dfCoop_withInts_Complexes_k562,
            paste0(Sys.Date(),"_s10_PublishedPPIandProtComplexes_k562_pe.txt"),
            sep = "\t",row.names = FALSE, quote = FALSE)

write.table(dfCoop_withInts_Complexes_hepg2,
            paste0(Sys.Date(),"_s10_PublishedPPIandProtComplexes_hepg2_pe.txt"),
            sep = "\t",row.names = FALSE, quote = FALSE)