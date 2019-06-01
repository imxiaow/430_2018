#============================================================================
# required package gProfileR

if (!require(gProfileR)) {
  install.packages("gProfileR")
  library(gProfileR) 
}

# installing/loading the package:
if(!require(installr)) {
  install.packages("installr")
  require(installr)
} #load / install+load installr
# using the package:

# -Jan 3rd, 
# There's error when installing the package, 
# This is becuase ORACLE-JAVA-7 is missing from the current system (Linux Ubuntu), 
# Java 8 is installed before because of requirements for cytoscape.
# Problem solved by manully download JDK-JAVA7 tar.gz file from JAVA web,
# By place it at /var/cache/oracle-jdk7-installer/ folder,
# untar it and follow the instructions in below link.
# http://www.webupd8.org/2012/01/install-oracle-java-jdk-7-in-ubuntu-via.html
# sudo add-apt-repository ppa:webupd8team/java
# sudo apt-get update
# sudo apt-get install oracle-java7-installer

# Check availble java verions by:
#       sudo update-alternatives --config java
# Change your perfered version by select its number.

# So that you can change back and forth between java 7 and java 8
# under different circumstances. 

#====================================

setwd("/home/xwang/Documents/BCB430/data_oct2nd/cor_results/NEAT1")
load("NEAT1_appear1.Rda")
load("NEAT1_appear2.Rda")
load("NEAT1_appear3.Rda")
load("NEAT1_appear4.Rda")
load("NEAT1_appear5.Rda")
load("NEAT1_appear6.Rda")
load("NEAT1_appear7.Rda")

load("NEAT1_RA_result.Rda")

setwd("/home/xwang/Documents/BCB430/data_oct2nd/cor_results/MALAT1")
load("MALAT1_appear1.Rda")
load("MALAT1_appear2.Rda")
load("MALAT1_appear3.Rda")
load("MALAT1_appear4.Rda")
load("MALAT1_appear5.Rda")
load("MALAT1_appear6.Rda")

load("MALAT1_RA_result.Rda")

#================================================
setwd("/home/xwang/Documents/BCB430/data_oct2nd/")
ra_MALAT1$Name <- as.character(ra_MALAT1$Name)
ra_NEAT1$Name <- as.character(ra_NEAT1$Name)



#=============================
# Which groups of genes are used into pathway enriment analysis?
# Do we have to set a cutoff of RA score or a threshold?
# Are we using the whole list?
# In this case, is gProfileR still a suitable method to use?
# Importantly, How does my analysis solve the ultimate goal of research?
#   "Find if the potential target genes have some ChIP-seq binding sites of transcription factors."
#   "gene-regulatory networks of lncRNAs in multiple tumor types"
#
#   "Identify TFs that are highlighted as target genes of different lncRNAs, and time permitting, analyse
#   binding sites of TFs in target gene promoters."
#
#
#===================================================================
# For now,just a try to use most appeared group of genes for pathway enrichment analysis
# E.g. NEAT1_appear7

#==============


#write(NEAT1_appear7, file="/home/xwang/Documents/BCB430/data_oct2nd/NEAT1_ap7_test.txt", sep="\t")
#NEAT1_a7 <- readLines("/home/xwang/Documents/BCB430/data_oct2nd/NEAT1_ap7_test.txt")
#NEAT1_a7_gp <- gprofiler(NEAT1_a7, ordered_query = FALSE)


NEAT1_RAed_glist <- ra_NEAT1$Name[which(ra_NEAT1$Score < 0.25)] # 0.25 threshold, 876 genes
NEAT1_gp <- gprofiler(NEAT1_RAed_glist, organism = "hsapiens", ordered_query = TRUE) # gprofiler success

MALAT1_RAed_glist <- ra_MALAT1$Name[which(ra_MALAT1$Score < 0.25)] # 0.25 threshold, 952 genes
MALAT1_gp <- gprofiler(MALAT1_RAed_glist, organism = "hsapiens", ordered_query = TRUE) # gprofiler 

# read this file into R, keep
pathways <- readLines("gmt_ids.txt")
NEAT1_gp_filter <- NEAT1_gp[NEAT1_gp$term.id %in% pathways, ]
MALAT1_gp_filter <- MALAT1_gp[MALAT1_gp$term.id %in% pathways, ]


#============
# 1. meeting every Tuesday with Juri
# 2. list 1 -7 in 7 tumor types, filter pcg with low-expressions? Check!
# 3. MALAT1 <-> NEAT1 are highly associated to each other, Did they share lots of common target genes? Check!

# 4. Compare differences in number of pathways and number of merged genes before and after filter out noise. 
        # First record number of pathways and number of merged genes right now. 
        # filter out noise, Is there a different? changes?
# 5. Noise: check if there's cases when,  positively correlated in 1 cantype at the same time negatively correlated with another cancer type.

# 6. Redo plotting, maybe you can Log the RA-score (treated as p-values),
                #    Add numbers on the plot for each group. 

# 7. RA-list that input into gprofiler, use 0.05 as cutoff, what's the number of genes? Answer: both of them about ~400 genes

#========
# 7 

NEAT1_RAed_glist_cutoff005 <- ra_NEAT1$Name[which(ra_NEAT1$Score < 0.05)] # 0.05 threshold, 463 genes
write.table(NEAT1_RAed_glist_cutoff005, file="NEAT1_genelistforgprofiler_PCAWG.txt", quote=FALSE,append=FALSE,sep='\n', row.names = FALSE, col.names = FALSE)


MALAT1_RAed_glist_cutoff005 <- ra_MALAT1$Name[which(ra_MALAT1$Score < 0.05)] # 0.05 threshold,  440 genes
write.table(MALAT1_RAed_glist_cutoff005, file="MALAT1_genelistforgprofiler_PCAWG.txt",quote=FALSE, append=FALSE,sep='\n', row.names = FALSE, col.names = FALSE)

library(gProfileR)
NEAT1_gprofiler_cutoff005_PCAWG <- gprofiler(NEAT1_RAed_glist_cutoff005, organism = "hsapiens", ordered_query = TRUE) # gprofiler success, 93 obs.
MALAT1_gprofiler_cutoff005_PCAWG <- gprofiler(MALAT1_RAed_glist_cutoff005, organism = "hsapiens", ordered_query = TRUE) # gprofiler, 88 obs.


common_fd_g <- intersect(NEAT1_RAed_glist_cutoff005, MALAT1_RAed_glist_cutoff005) # 259
write.table(common_fd_g, file="common_genelistforgprofiler_PCAWG.txt", quote=FALSE, append=FALSE,sep='\n', row.names = FALSE, col.names = FALSE)

common_gene_PCAWG_NM_PATH <- gprofiler(common_fd_g, organism = "hsapiens", ordered_query = TRUE) # 
save(common_gene_PCAWG_NM_PATH, file="common_gene_NM_gprofiler_result_PCAWG.txt")

#===========

all_pathways <- readLines('gmt_ids.txt') #17335 GO ID
common_NM_pathways <- read.table("/Users/xiaowang/Desktop/gprofiler_results_common_gene_PCAWG.txt", sep="\t", header = TRUE, quote = "")
#97 GOID in total 
NEAT1_pathway_result<-read.table("/Users/xiaowang/Desktop/gprofiler_results_NEAT1_PCAWG.txt", sep="\t", header = TRUE, quote = "")
#135

Final_common_NM_pathways <- common_NM_pathways[which(common_NM_pathways$GO.ID %in% all_pathways),]
write.table(Final_common_NM_pathways, file="/Users/xiaowang/Desktop/final_common_genelistforgprofiler_PCAWG.txt", sep='\t',quote=FALSE,row.names = TRUE, col.names = TRUE)



#========
# END