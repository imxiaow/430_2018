#==================================================
if (! require(biomaRt, quietly=TRUE)) {
  if (! exists("biocLite")) {
    source("https://bioconductor.org/biocLite.R")
  }
  biocLite("biomaRt")
  library(biomaRt)
}
if (!require(igraph, quietly=TRUE)) {
  install.packages("igraph")
  library(igraph)
}
if (!require(readr, quietly=TRUE)) {
  install.packages("readr")
  library(readr)
}
# Install Packages
#==================================================
# Load the data
#================================================================
myMart <- useMart("ensembl", dataset="hsapiens_gene_ensembl")
(filters <- listFilters(myMart))
(attributes <- listAttributes(myMart))          # View possible filters and attributes in Biomart
filters[grep("ENSG", filters$description), ]    # Specific filter for ensemble IDs (ensemble_gene_id)
attributes[grep("symbol", attributes$description, ignore.case=TRUE), ] # Attribute description for HUGO symbol for mapping

#=============================
test_result_temp <- list()
result <- c()
for (ID in shared_geneID_005_NM) {
  print(ID)
  # retrieves information from the Ensembl database. by  <getBM>.
  test_result_temp[[ID]] <- getBM(filters = "ensembl_gene_id",
                                  attributes = c("hgnc_symbol"),
                                  values = ID,
                                  mart = myMart)
  
  cat(sprintf("Ensemble gene ID - %s \n ", ID)) # Format output
  cat(sprintf("HGNC Symbol - %s \n ", test_result_temp[[ID]]["hgnc_symbol"][1,]))
  cat("\n")
  need_store <- test_result_temp[[ID]]["hgnc_symbol"][1,]
  #  print(test_result_temp[[ID]]["hgnc_symbol"][1,])
  result <- c(result, need_store)
}

length(result) #1874
length(which(is.na(result))) # 57 
# make it into dictionary
dict_TCGA_NM_result <- as.data.frame(cbind(shared_geneID_005_NM, result), stringsAsFactors=FALSE)
save(dict_TCGA_NM_result, file="dict_TCGA_NM_result_gene_ID_symbol.Rda")
# check number of NAs in the results, 
length(which(is.na(dict_TCGA_NM_result$result))) #57
# dealing with NAs. 
sel <-which(!is.na(dict_TCGA_NM_result$result))

newEnsemID <- dict_TCGA_NM_result$shared_geneID_005_NM[sel]
newHUGOid <- dict_TCGA_NM_result$result[sel]
newGeneSymToEnsm <- as.data.frame(cbind(newEnsemID, newHUGOid), stringsAsFactors = FALSE)
# compare how many gene symbols in the 259 df. 
length(which(duplicated(newGeneSymToEnsm$newHUGOid))) # 0 
length(unique(newGeneSymToEnsm$newHUGOid))# 1817
# compare result. 
length(which(newGeneSymToEnsm$newHUGOid %in% common_fd_g)) #224
length(which(common_fd_g%in% newGeneSymToEnsm$newHUGOid)) #224

length(which(!common_fd_g%in% newGeneSymToEnsm$newHUGOid))#35

# what are these common genes results between these 2 database? what are they?
# what are the different genes (e.g. in PWCAGE but not in TCGA)
common_result_2_database <- newGeneSymToEnsm$newHUGOid[which(newGeneSymToEnsm$newHUGOid %in% common_fd_g)]
result_in_PWCAG_not_TCGA <- common_fd_g[which(!common_fd_g %in% newGeneSymToEnsm$newHUGOid)]

save(common_result_2_database, file = "common_result_2_database.Rda")
save(result_in_PWCAG_not_TCGA, file= "result_in_PWCAG_not_TCGA.Rda")
save(newGeneSymToEnsm, file="result_dict_no_NA.Rda")
write.table(common_result_2_database, "common_result_2_database.txt", sep="\t")

#=========
library(gProfileR)
gprofiler(newGeneSymToEnsm$newHUGOid, organism = "hsapiens",ordered_query = TRUE)

newGeneSymToEnsm$newHUGOid
write.table(newGeneSymToEnsm$newHUGOid, file="/Users/xiaowang/Desktop/common_genelistforgprofiler_TCGA.txt", quote=FALSE, append=FALSE,sep='\n', row.names = FALSE, col.names = FALSE)


all_pathways <- readLines('gmt_ids.txt') #17335 GO ID
common_NM_pathways_TCGA <- read.table("/Users/xiaowang/Desktop/gprofiler_results_common_TCGA.txt", sep="\t", header = TRUE, quote = "")
#294 GO total

Final_common_NM_pathways_TCGA <- common_NM_pathways_TCGA[which(common_NM_pathways_TCGA$GO.ID %in% all_pathways),] #206 in total
write.table(Final_common_NM_pathways_TCGA, file="/Users/xiaowang/Desktop/final_common_genelistforgprofiler_TCGA.txt", sep='\t',quote=FALSE,row.names = FALSE, col.names = TRUE)

#==========
#END