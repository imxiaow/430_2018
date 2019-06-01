#====================================== Load all the data provided==========================================================================

lncRNAwithpatient <- readRDS("/home/xwang/Documents/BCB430/data_oct2nd/data/lnc_RNAseq_wPatientData_RankAg.rds")
View(lncRNAwithpatient)

proteincodingwithpatient <- readRDS("/home/xwang/Documents/BCB430/data_oct2nd/data/pcg_RNAseq_wPatientData_RankAg.rds")
View(proteincodingwithpatient)

highexpresslncRNAwithcancertypes <-  read.table("/home/xwang/Documents/BCB430/data_oct2nd/data/high_lncsmed4top5cancersPCAWG.txt", header=TRUE, sep=";", dec=".",stringsAsFactors = FALSE)
View(highexpresslncRNAwithcancertypes)

allpatientclinical = readRDS("data/Jan26_PCAWG_clinical")
View(allpatientclinical) #total 1228 patients with 45variables(personal health info)



#======================================== Explore the values of datafiles =======================================================================

#====  high_lncsmed4top5cancersPCAWG.txt (highexpresslncRNAwithcancertypes)  =======================================
str(highexpresslncRNAwithcancertypes)

#levels(highexpresslncRNAwithcancertypes$canc) #to view all the types of returns in cancer category
#[1] "Breast Infiltrating duct carcinoma"      "Kidney Adenocarcinoma, chromophobe type" "Kidney Adenocarcinoma, clear cell type" 
#[4] "Kidney Adenocarcinoma, papillary type"   "Liver Hepatocellular carcinoma"          "Ovary Serous cystadenocarcinoma"        
#[7] "Pancreas Pancreatic ductal carcinoma"   
# Total 7 cancer types contained in  high_lncsmed4top5cancersPCAWG.txt file. 

#levels(highexpresslncRNAwithcancertypes$gene) #to view all the genes in this file
#[1] "AC004447.2"      "AC005083.1"      "AC006126.4"      "AC009336.24"     "AC010127.5"      "AC011523.2"      "AC013394.2"      "AC017048.3"     
#[9] "AC017060.1"      "AC018730.1"      "AC019117.1"      "AC019181.2"      "AC034220.3"      "AC073046.25"     "AC073218.3"      "AC079630.2"     
#[17] "AC079630.4"      "AC093323.3"      "AC093673.5"      "AC093850.2"      "AC103563.8"      "AC106869.2"      "AC108488.3"      "AC113189.5"     
#[25] "AC133528.2"      "AC147651.5"      "ADAMTS9-AS1"     "ADORA2A-AS1"     "AFAP1-AS1"       "AP000355.2"      "AP001626.1"      "ASB16-AS1"      
#[33] "BAIAP2-AS1"      "C11orf95"        "C1RL-AS1"        "CASC7"           "CD27-AS1"        "CITF22-24E5.1"   "CRNDE"           "CTA-29F11.1"    
#[41] "CTB-131K11.1"    "CTB-25B13.12"    "CTB-27N1.1"      "CTB-43E15.4"     "CTC-503J8.6"     "CTD-2015H6.3"    "CTD-2037K23.2"   "CTD-2081K17.2"  
#[49] "CTD-2228K2.7"    "CTD-2369P2.2"    "CTD-2377D24.6"   "CTD-2377O17.1"   "CTD-2540F13.2"   "CTD-3184A7.4"    "DANCR"           "DBH-AS1"        
#[57] "DGCR5"           "DYNLL1-AS1"      "EGOT"            "EMX2OS"          "EPB41L4A-AS1"    "FAM83H-AS1"      "FGD5-AS1"        "FGF14-AS2"      
#[65] "FOXN3-AS1"       "GAS6-AS2"        "GS1-251I9.4"     "HOTAIRM1"        "HOXA-AS2"        "HOXB-AS3"        "HOXD-AS1"        "ILF3-AS1"       
#[73] "KMT2E-AS1"       "LBX2-AS1"        "LINC00035"       "LINC00087"       "LINC00094"       "LINC00152"       "LINC00261"       "LINC00462"      
#[81] "LINC00493"       "LINC00657"       "LINC00665"       "LINC00887"       "LINC00958"       "LINC00963"       "LINC00982"       "LINC00984"      
#[89] "LINC01003"       "LINC01018"       "MALAT1"          "MAPKAPK5-AS1"    "MCF2L-AS1"       "MFI2-AS1"        "MIR210HG"        "MIR663A"        
#[97] "MMP24-AS1"       "NCBP2-AS2"       "NEAT1"           "OIP5-AS1"        "OSER1-AS1"       "PAXIP1-AS1"      "PCAT6"           "PITPNA-AS1"     
#[105] "RAB11B-AS1"      "RBPMS-AS1"       "RMRP"            "RN7SL1"          "RNU12"           "RP11-1094M14.11" "RP11-111M22.3"   "RP11-1151B14.3" 
#[113] "RP11-118K6.3"    "RP11-122G18.5"   "RP11-139H15.1"   "RP11-145M9.4"    "RP11-156K13.1"   "RP11-157P1.4"    "RP11-166P13.3"   "RP11-18I14.10"  
#[121] "RP11-20I20.4"    "RP11-218M22.1"   "RP11-220I1.1"    "RP11-23P13.6"    "RP11-252A24.7"   "RP11-254F7.2"    "RP11-283G6.4"    "RP11-284F21.9"  
#[129] "RP11-285F7.2"    "RP11-290D2.6"    "RP11-295G20.2"   "RP11-296I10.3"   "RP11-2N1.2"      "RP11-304L19.5"   "RP11-305K5.1"    "RP11-309L24.4"  
#[137] "RP11-317J10.2"   "RP11-325F22.5"   "RP11-347C12.10"  "RP11-349A22.5"   "RP11-350J20.12"  "RP11-354P11.2"   "RP11-357H14.17"  "RP11-395B7.2"   
#[145] "RP11-395G23.3"   "RP11-395P17.3"   "RP11-421F16.3"   "RP11-422J8.1"    "RP11-427H3.3"    "RP11-440D17.3"   "RP11-469A15.2"   "RP11-47A8.5"    
#[153] "RP11-48O20.4"    "RP11-490M8.1"    "RP11-496I9.1"    "RP11-49I11.1"    "RP11-532F12.5"   "RP11-553L6.5"    "RP11-572C15.6"   "RP11-620J15.3"  
#[161] "RP11-622A1.2"    "RP11-622K12.1"   "RP11-644F5.11"   "RP11-65J3.1"     "RP11-713M15.2"   "RP11-728F11.3"   "RP11-728F11.4"   "RP11-732M18.3"  
#[169] "RP11-798K3.2"    "RP11-807H17.1"   "RP11-834C11.4"   "RP11-923I11.6"   "RP1-257A7.5"     "RP13-582O9.5"    "RP1-39G22.7"     "RP1-60O19.1"    
#[177] "RP4-635E18.8"    "RP4-639F20.1"    "RP4-758J18.10"   "RP4-797C5.2"     "RP5-1085F17.3"   "RP5-1103G7.4"    "RP5-1112D6.4"    "RP5-1120P11.1"  
#[185] "RP5-887A10.1"    "RP6-109B7.3"     "RPPH1"           "SBF2-AS1"        "SCARNA2"         "SCARNA9"         "SEMA3B-AS1"      "SERTAD4-AS1"    
#[193] "SMIM2-AS1"       "SNHG1"           "SNHG11"          "SNHG12"          "SNHG15"          "SNHG16"          "SNHG18"          "SNHG6"          
#[201] "SNHG7"           "SNHG8"           "SNHG9"           "SNORA76"         "SNORD3A"         "TAPSAR1"         "THAP9-AS1"       "TRIM52-AS1"     
#[209] "U47924.27"       "UNC5B-AS1"       "WAC-AS1"         "ZBED5-AS1"       "ZFAS1"           "ZNF503-AS2"      "ZNF667-AS1"    
# Total 215 unique highly expressed lncRNAs contained in  high_lncsmed4top5cancersPCAWG.txt file. 
#but total 493 entries(duplicated cancer types for each gene)
#*****Question1: What's the median in this data file for?? what is the median value represent?
# the medium of FPKM of the RNA-Seq data


#=====  lnc_RNAseq_wPatientData_RankAg.rds  (lncRNAwithpatient)  ============================================
PatientsID_Inc <- rownames(lncRNAwithpatient) #Patient IDs, total 485

length(colnames(lncRNAwithpatient))  
#total 5607 lncRNAs, with last 2 column for cancaer types (canc) [5608], and patient IDs[5609]
#[1] 5609

lncRNAwithpatient$canc 
#485 patients with its cancer type



#==== pcg_RNAseq_wPatientData_RankAg.rds (proteincodingwithpatient)  ==========================================
PatientsID_pcg <- rownames(proteincodingwithpatient)#Patient IDs, total 485

length(colnames(lncRNAwithpatient))  
#total 19953 genes with last two column [19954] for canc, [19955] for patient IDs
#[1] 19955


#*****Question2: RNA-Seq data for each pcg genes and lncRNAs, is it normalized? or raw data? P-value?
#Answer: the FPKM of RNA-Seq data


#==================== Jan26_PCAWG_clinical =================================================
#allpatientclinical = readRDS("data/Jan26_PCAWG_clinical")
#View(allpatientclinical) #total 1228 patients with 45variables(personal health info)



#==================== Preparing cancer specific  data =========================================================================================
canctype1 <- "Breast Infiltrating duct carcinoma" 
canctype2 <-"Kidney Adenocarcinoma, chromophobe type"
canctype3 <- "Kidney Adenocarcinoma, clear cell type" 
canctype4 <- "Kidney Adenocarcinoma, papillary type"
canctype5 <- "Liver Hepatocellular carcinoma"
canctype6 <- "Ovary Serous cystadenocarcinoma" 
canctype7 <- "Pancreas Pancreatic ductal carcinoma" 




#===================  Task1: for each types of cancer/tumor , get its correspond petient IDs 
#(cancertype -> pateint IDs (both lncRNAs and pcg))

#lncRNAwithpatient[,5608:5609] == proteincodingwithpatient[,19954:19955]
# use either one is fine
patientIDswithcanctypes <- proteincodingwithpatient[,19954:19955] #create variable contain patient Ids with canctypes
save(patientIDswithcanctypes,file="patientIDswithcanctypes.Rda") 

# can load "patientIDswithcanctypes.Rda" directly when needed
load("patientIDswithcanctypes.Rda")
#now munipulate patient ids for each cancer type
#   How to avoid using for loops..... apply()
n1 <-0
n2 <-0
n3 <-0
n4 <-0
n5 <-0
n6 <-0
n7 <-0
for (i in 1:length(patientIDswithcanctypes$canc)){
  if(patientIDswithcanctypes$canc[i] == canctype1){
    n1 <- n1 +1
  }
  if(patientIDswithcanctypes$canc[i] == canctype2){
    n2 <- n2+1
  }
  if(patientIDswithcanctypes$canc[i] == canctype3){
    n3 <- n3+1
  }
  if(patientIDswithcanctypes$canc[i] == canctype4){
    n4 <- n4 + 1
  }
  if(patientIDswithcanctypes$canc[i] == canctype5){
    n5 <- n5 +1
  }
  if(patientIDswithcanctypes$canc[i] == canctype6){
    n6 <-n6+1
  }
  if(patientIDswithcanctypes$canc[i] == canctype7){
    n7 <- n7 +1
  }
}

# write a general function for each cancer type
get_patient_cancer <- function(dataset, cancertypeX){
  cancertypeX.patientIDs <- vector("list")
  i = 1
  for (cancertype in 1:length(dataset$canc)){
    if (dataset$canc[cancertype] == cancertypeX){
      element <- dataset$patient[cancertype]
     # print(element)
    #  print(0)
      cancertypeX.patientIDs[i] <- element
      i = i + 1
    }
  }
  return (cancertypeX.patientIDs)
}
#test function, 1
canctype1.patientIDs <- get_patient_cancer(patientIDswithcanctypes, canctype1) #match the total number
canctype2.patientIDs <- get_patient_cancer(patientIDswithcanctypes, canctype2)
canctype3.patientIDs <- get_patient_cancer(patientIDswithcanctypes, canctype3)
canctype4.patientIDs <- get_patient_cancer(patientIDswithcanctypes, canctype4)
canctype5.patientIDs <- get_patient_cancer(patientIDswithcanctypes, canctype5)
canctype6.patientIDs <- get_patient_cancer(patientIDswithcanctypes, canctype6)
canctype7.patientIDs <- get_patient_cancer(patientIDswithcanctypes, canctype7)

# how to use apply or lapply function to solve this?
#all_types_canc <- c(canctype1, canctype2, canctype3, canctype4, canctype5, canctype6, canctype7)




#================     Task2: for each cancer/tumor type, get the higly expressed lncRNAs. 
highlncwithcanctypes <- highexpresslncRNAwithcancertypes[,2:3] #manipulate only the genename and cancertype
save(highlncwithcanctypes,file="highlncwithcanctypes.Rda") 
load("highlncwithcanctypes.Rda")

g1 <-0
g2 <-0
g3 <-0
g4 <-0
g5 <-0
g6 <-0
g7 <-0
#see if the number is correct
for (i in 1:length(highlncwithcanctypes$canc)){
  if(highlncwithcanctypes$canc[i] == canctype1){
    g1 <- g1 +1
  }
  if(highlncwithcanctypes$canc[i] == canctype2){
    g2 <- g2+1
  }
  if(highlncwithcanctypes$canc[i] == canctype3){
    g3 <- g3+1
  }
  if(highlncwithcanctypes$canc[i] == canctype4){
    g4 <- g4 + 1
  }
  if(highlncwithcanctypes$canc[i] == canctype5){
    g5 <- g5 +1
  }
  if(highlncwithcanctypes$canc[i] == canctype6){
    g6 <-g6+1
  }
  if(highlncwithcanctypes$canc[i] == canctype7){
    g7 <- g7 +1
  }
}

# write a general function for each cancer type
get_lncRNA_cancer <- function(dataset, cancertypeX){
  cancertypeX.geneIDs<- vector("list")
  i = 1
  for (cancertype in 1:length(dataset$canc)){
    if (dataset$canc[cancertype] == cancertypeX){
      element <- dataset$gene[cancertype]
      #print(elemecanctype1_lncRNAgeneIDs.Rdant)
      #print(0)
      cancertypeX.geneIDs[i] <- element
      i = i + 1
    }
  }
  return (cancertypeX.geneIDs)
}

#test function 2
canctype1.geneIDs <- get_lncRNA_cancer(highlncwithcanctypes, canctype1) #match the total number
canctype2.geneIDs <- get_lncRNA_cancer(highlncwithcanctypes, canctype2)
canctype3.geneIDs <- get_lncRNA_cancer(highlncwithcanctypes, canctype3)
canctype4.geneIDs <- get_lncRNA_cancer(highlncwithcanctypes, canctype4)
canctype5.geneIDs <- get_lncRNA_cancer(highlncwithcanctypes, canctype5)
canctype6.geneIDs <- get_lncRNA_cancer(highlncwithcanctypes, canctype6)
canctype7.geneIDs <- get_lncRNA_cancer(highlncwithcanctypes, canctype7)

# how to use apply or lapply function to solve this?
#all_types_canc <- c(canctype1, canctype2, canctype3, canctype4, canctype5, canctype6, canctype7)






#========================= Save the geneID lists and the patientID lists for all cancer types ========================================
#load them when needed in the futher
save(canctype1.geneIDs, file= "canctype1_lncRNAgeneIDs.Rda")
save(canctype2.geneIDs, file= "canctype2_lncRNAgeneIDs.Rda")
save(canctype3.geneIDs, file= "canctype3_lncRNAgeneIDs.Rda")
save(canctype4.geneIDs, file= "canctype4_lncRNAgeneIDs.Rda")
save(canctype5.geneIDs, file= "canctype5_lncRNAgeneIDs.Rda")
save(canctype6.geneIDs, file= "canctype6_lncRNAgeneIDs.Rda")
save(canctype7.geneIDs, file= "canctype7_lncRNAgeneIDs.Rda")

save(canctype1.patientIDs, file= "canctype1_patientIDs.Rda")
save(canctype2.patientIDs, file= "canctype2_patientIDs.Rda")
save(canctype3.patientIDs, file= "canctype3_patientIDs.Rda")
save(canctype4.patientIDs, file= "canctype4_patientIDs.Rda")
save(canctype5.patientIDs, file= "canctype5_patientIDs.Rda")
save(canctype6.patientIDs, file= "canctype6_patientIDs.Rda")
save(canctype7.patientIDs, file= "canctype7_patientIDs.Rda")





#=========================  Task3: For the same cancer type, first use one lncRNA generate pipline  =============================
# Use cancertype1 as test
load("canctype1_lncRNAgeneIDs.Rda")
load("canctype1_patientIDs.Rda")
# use one of the high expressed lncRNA in this cancertype1, for all patients(80), find the correlation of pcgs(filter) with this lncRNA.

first_lncRNA <- canctype1.geneIDs[1]
first_lncRNA <- as.character(first_lncRNA) # change the type of variable from list to character

# manipulate 80 patients data for both lncRNA and pcg RNA-Seq, this can be further used for other high_expressed lncRNA

# 1) lncRNAs RNA-Seq data
lncRNA_with_patientIdsfor_cancertype1 <- subset(lncRNAwithpatient, is.element(lncRNAwithpatient$patient, canctype1.patientIDs))
# 80 patients with 5607 lncRNAs RNA-Seq data

highexpr_lncRNA_with_patientIds_cancertype1 <-subset(lncRNA_with_patientIdsfor_cancertype1, select = is.element(colnames(lncRNA_with_patientIdsfor_cancertype1),canctype1.geneIDs))
# 80 patients with their 47 highly expressed lncRNA RNA-Seq data

save(highexpr_lncRNA_with_patientIds_cancertype1, file = "cancertype1_highexpr_lncRNAs_wpatients.Rda")
# can directly load later
#We have to log1p the lncRNA patient dataset.
log1p_highexpr_lncRNA_with_patientIds_canctype1 <- as.data.frame(apply(highexpr_lncRNA_with_patientIds_cancertype1,2, log1p))
save(log1p_highexpr_lncRNA_with_patientIds_canctype1, file= "log1p_cancertype1_highexpr_lncRNAs_wpatient.Rda")

# 2) Protein-coding genes RNA-Seq data
pcg_with_patientIdsfor_cancertype1 <- subset(proteincodingwithpatient, is.element(lncRNAwithpatient$patient, canctype1.patientIDs))
# 80 patients with 19955 pcgs RNA-Seq data
save(pcg_with_patientIdsfor_cancertype1, file = "cancertype1_pcg_wpatients.Rda")
# load it directly later


#================
#task continue in next script. 

#END
