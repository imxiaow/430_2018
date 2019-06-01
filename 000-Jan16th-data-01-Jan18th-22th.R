#======================================================================================================================
# '''I downloaded the individual patient files from the latest TCGA release and merged them into one matrix. 
#    I tried downloading all cancers at once so one file has most cancer types with at least 90 patients in each one 
#     and then GBM and PAAD have their own files separately. ''' - By Karina 
#
# (1) Merged files downloaded from TCGA:  
#                                     - 9246rnaSEQfiles.rds
#
# (2) Protein coding gene expression processed by Karina where additional columns are added with patient information.
#                               - 19438_lncRNAs_tcga_all_cancers_Dec20_wclinical_data.rds
#
# (3)lncRNA gene expression data again processed the same way.
#                               - 5919_lncRNAs_tcga_all_cancers_Dec20_wclinical_data.rds
#
# (4) GBM files.
#             - 167rnaSEQfiles.rds
# 
# (5) PAAD files.
#             - 182rnaSEQfiles_PAAD.rds
#
# For (4) & (5), the UUIDs would need to be changed to TCGA IDs. Following files are helpful for doing that.
#
#
#   - Sample code, to see whether RNA-Seq sample is from tumor or normal tissue: 
#             TCGA_sample_codes.csv
#   - Confirm these patients are NOT also in PCAWG:           
#             TCGA_IDs_usedinPCAWG.txt
#   - TSS codes to figure out cancer code in the patient ID:  
#             TCGA_TissueSourceSite_Codes2017.csv
#   - Clinical data:                        
#             all_clin_XML_tcgaSept2017.csv
#
#======================================================================================================
# Objective: Load each data, familiarize dfs, record basic infos.

merged_file <- readRDS("~/Documents/BCB430/data2_Jan/TCGA_data/TCGA Data/9246rnaSEQfiles.rds")



#==== D2 ===========================================================
all_pcgs <- readRDS("~/Documents/BCB430/data2_Jan/TCGA_data/TCGA Data/19438_lncRNAs_tcga_all_cancers_Dec20_wclinical_data.rds")
# 7387 x  19442
# TCGA-ID x Ensembl-ID (pcg)
#========================
# 19439-> tumor types,
# 19440 ->?
# 19441 -> survival
# 19442 -> pateint's sex
#=== Determine all types of cancers in df ======
#> levels(as.factor(d2[,19439]))
#[1] ""                                                                
#[2] "Bladder Urothelial Carcinoma"                                    
#[3] "Brain Lower Grade Glioma"                                        
#[4] "Breast invasive carcinoma"                                       
#[5] "Cervical squamous cell carcinoma and endocervical adenocarcinoma"
#[6] "Colon adenocarcinoma"                                            
#[7] "Head and Neck squamous cell carcinoma"                           
#[8] "Kidney renal clear cell carcinoma"                               
#[9] "Kidney renal papillary cell carcinoma"                           
#[10] "Liver hepatocellular carcinoma"                                  
#[11] "Lung adenocarcinoma"                                             
#[12] "Lung squamous cell carcinoma"                                    
#[13] "Ovarian serous cystadenocarcinoma"                               
#[14] "Prostate adenocarcinoma"                                         
#[15] "Sarcoma"                                                         
#[16] "Skin Cutaneous Melanoma"                                         
#[17] "Stomach adenocarcinoma"                                          
#[18] "Thyroid carcinoma"                                               
#[19] "Uterine Corpus Endometrial Carcinoma"     
pcgs_with_canctypes <- all_pcgs[,1:19439]
save(pcgs_with_canctypes, file="all_pcgs_with_canctypes_7387p.Rda")
#=========================================
# Compare the tumor types in PCAWG and TCGA,
#   1) "Pancreas Pancreatic dutal carcinoma" is missing in TCGA,
#   2) "Kidney Adenocarcinoma" - chromophobe type is missing in TCGA.


#===== D3 ==================================================================
all_lncRNAs <- readRDS("~/Documents/BCB430/data2_Jan/TCGA_data/TCGA Data/5919_lncRNAs_tcga_all_cancers_Dec20_wclinical_data.rds")
# 7387 x 5923
# TCGA-ID x Ensembl-ID (lncRNA)
#========================
# 5920-> tumor types,
# 5921 ->?
# 5922 -> survival
# 5923 -> pateint's sex
#=============================
# Same tumor type as D2. 



#==== GMB file ==============================================================
readRDS("~/Documents/BCB430/data2_Jan/TCGA_data/TCGA Data/167rnaSEQfiles.rds")

#==== PADD file ==============================================================
readRDS("~/Documents/BCB430/data2_Jan/TCGA_data/TCGA Data/182rnaSEQfiles_PAAD.rds")



#========================================
# 1. find Ensembl IDs for NEAT1 and MALAT1, get this 2 lncRNA gene expressions for 7387.
# 2. check whether if patients that are in TCGA data NOT in PCAWG.
# 3. check how many patients are there for each cancer type??

# 4. Is the normal and canc tissue type matters? (source type)
# 5. what's TSS code for ?

# 6. filter patients with NULL cancertype, then log1p both datasets.(pcgs and lncRNAs)

# 7. for each cancer type, find significant correlated pcgs, ~20 list for each lncRNAs.
# 8. RA, combine them together, filter genes with RA-score <0.05

# 9. compare the result with PCAWGs.
#============================================
# Task1a: 
# Ensemble gene ID for NEAT1: ENSG00000245532
# Ensemble gene ID for MALAT1: ENSG00000251562 (this is the one work!)/ ENSG00000278217

# check if these 2 ensemble ID is in the lncRNA dataset. 
"ENSG00000245532" %in% colnames(all_lncRNAs)
#[1] TRUE
"ENSG00000251562" %in% colnames(all_lncRNAs)
#[1] TRUE
#=======
# TAsk 1b:
# get this 2 lncRNA list from whole dataset. 
NEAT1_all_p <- all_lncRNAs[which(colnames(all_lncRNAs) == "ENSG00000245532")]
MALAT1_all_p <-all_lncRNAs[which(colnames(all_lncRNAs) == "ENSG00000251562")]
# which(colnames(all_lncRNAs) == "ENSG00000245532") -->> NEAT1
# [1] 3548
# which(colnames(all_lncRNAs) == "ENSG00000251562") -->> MALAT1
# [1] 2154

# save this 2 lncRNA data, unfiltered, unprocessed. for later.
save(NEAT1_all_p, file="NEAT1_all_patients7387.Rda")
save(MALAT1_all_p, file= "MALAT1_all_patients7387.Rda")
#=============================================
# Task2:
# Check whether if patients that are in TCGA data NOT in PCAWG.

# read text file : TCGA_IDs_usedinPCAWG.txt
# check if there's TCGA patients ID used in PCAWG
# how many of them among 7387?

#"/home/xwang/Documents/BCB430/data2_Jan/TCGA_data/TCGA Data/"

check_patient_datafile <- read.table(file = "/home/xwang/Documents/BCB430/data2_Jan/TCGA_data/TCGA Data/TCGA_IDs_usedinPCAWG.txt", 
                header = TRUE,sep = " ", stringsAsFactors = FALSE)

length(which(check_patient_datafile[,1]%in%rownames(MALAT1_all_p)))
#checked, all TCGA patient not used in PCAWG
length(which(check_patient_datafile[,2] %in% colnames(merged_file[2:9247])))
#============================================
# Task3:
# Check # patients  for each cancer type.
n_typesofcancer <- length(levels(as.factor(pcgs_with_canctypes[,19439])))
typesofcancer <- levels(as.factor(pcgs_with_canctypes[,19439]))

c1 <- 0
c2 <- 0
c3 <- 0
c4 <- 0
c5 <- 0
c6 <- 0
c7 <- 0
c8 <- 0
c9 <- 0
c10 <- 0
c11 <- 0
c12 <- 0
c13 <- 0
c14 <- 0
c15 <- 0
c16 <- 0
c17 <- 0
c18 <- 0
c19 <- 0


for(i in 1:nrow(pcgs_with_canctypes)){
  if (as.factor(pcgs_with_canctypes[i,19439]) == typesofcancer[1]){
    c1 <- c1 + 1
  }
  if (as.factor(pcgs_with_canctypes[i,19439]) == typesofcancer[2]){
    c2 <- c2 + 1
  }
  if (as.factor(pcgs_with_canctypes[i,19439]) == typesofcancer[3]){
    c3 <- c3 + 1
  }
  if (as.factor(pcgs_with_canctypes[i,19439]) == typesofcancer[4]){
    c4 <- c4 + 1
  }
  if (as.factor(pcgs_with_canctypes[i,19439]) == typesofcancer[5]){
    c5 <- c5 + 1
  }
  if (as.factor(pcgs_with_canctypes[i,19439]) == typesofcancer[6]){
    c6 <- c6 + 1
  }
  if (as.factor(pcgs_with_canctypes[i,19439]) == typesofcancer[7]){
    c7 <- c7 + 1
  }
  if (as.factor(pcgs_with_canctypes[i,19439]) == typesofcancer[8]){
    c8 <- c8 + 1
  }
  if (as.factor(pcgs_with_canctypes[i,19439]) == typesofcancer[9]){
    c9 <- c9 + 1
  }
  if (as.factor(pcgs_with_canctypes[i,19439]) == typesofcancer[10]){
    c10 <- c10 + 1
  }
  if (as.factor(pcgs_with_canctypes[i,19439]) == typesofcancer[11]){
    c11 <- c11 + 1
  }
  if (as.factor(pcgs_with_canctypes[i,19439]) == typesofcancer[12]){
    c12 <- c12 + 1
  }
  if (as.factor(pcgs_with_canctypes[i,19439]) == typesofcancer[13]){
    c13 <- c13 + 1
  }
  if (as.factor(pcgs_with_canctypes[i,19439]) == typesofcancer[14]){
    c14 <- c14 + 1
  }
  if (as.factor(pcgs_with_canctypes[i,19439]) == typesofcancer[15]){
    c15 <- c15 + 1
  }
  if (as.factor(pcgs_with_canctypes[i,19439]) == typesofcancer[16]){
    c16 <- c16 + 1
  }
  if (as.factor(pcgs_with_canctypes[i,19439]) == typesofcancer[17]){
    c17 <- c17 + 1
  }
  if (as.factor(pcgs_with_canctypes[i,19439]) == typesofcancer[18]){
    c18 <- c18 + 1
  }
  if (as.factor(pcgs_with_canctypes[i,19439]) == typesofcancer[19]){
    c19 <- c19 + 1
  }
  
}

# Number of patients in each canc type:
#> c1 = 63
#> c2 = 382
#> c3 = 488
#> c4 = 987
#> c5 = 284
#> c6 = 396
#> c7 = 456
#> c8 = 484
#> c9 = 255
#> c10 = 319
#> c11 = 432
#> c12 = 453
#> c13 = 305
#> c14 = 472
#> c15 = 224
#> c16 = 99
#> c17 = 345
#> c18 = 453
#> c19 = 490
#> c1+c2+c3+c4+c5+c6+c7+c8+c9+c10+c11+c12+c13+c14+c15+c16+c17+c18+c19 = 7387

n_each_canc <- c(c1,c2,c3,c4,c5,c6,c7,c8,c9,c10,c11,c12,c13,c14,c15,c16,c17,c18,c19)

# save types of cancer (names), and its numbers for later use. 
cancertypes_plus_numbers <- as.data.frame(cbind(typesofcancer, n_each_canc), stringsAsFactors = FALSE)
cancertypes_plus_numbers$n_each_canc <- as.integer(cancertypes_plus_numbers$n_each_canc)
save(cancertypes_plus_numbers, file="cancertypes_and_itsnumberofpatietns.Rda")



#===========
# later loading 
#load("~/Documents/BCB430/data2_Jan/TCGA_data/NEAT1_all_patients7387.Rda")
#load("~/Documents/BCB430/data2_Jan/TCGA_data/MALAT1_all_patients7387.Rda")
#load("~/Documents/BCB430/data2_Jan/TCGA_data/cancertypes_and_itsnumberofpatietns.Rda")
#load("~/Documents/BCB430/data2_Jan/TCGA_data/all_pcgs_with_canctypes_7387p.Rda")
#NEAT1_all_p
#MALAT1_all_p
#pcgs_with_canctypes
#cancertypes_plus_numbers


#=========
#END