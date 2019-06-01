#================================================================================================
# As Rscript 003 already compute significant pcgs for all lncRNAs associated with cancertype1, 
# This script aim to apply the same pipeline to all cancer types.
#
# Input file would be: 2 types
# 1) processed lncRNAs RNA-seq data from patients FOR each cancertype, (7 files in total)
# for example, "log1p_cancertype1_highexpr_lncRNAs_wpatient.Rda" -> 47 lncRNA with patients df
# 2) 17808 pcgs with patietns df (this is the same for all canctype)
# 
# Output: 7 files for its cancertypes
# Each highly expressed lncRNAs would have 3 files generated:
# 1) result (all pcg, need this for fdr), contain coeffient, p-value,intercept, genename (raw-result)
# 2) sig_result (after fdr), contain coeffient, p-value, fdr, intercept, genename   
# 3) patient data with only significant pcgs. (convenient for pllotting)
#
# Task 1: preprocess the data for rest of the cancer types. e.g. log1p, filter/subset ...
# Task 2: Consider parrallel/multiprocess by running this script with different cancer types.
#         The run time would be extremly high if compute one by one. 
#
#=========================  Previously defined functions  =============================================================
#helper function 1
col_mean_is_greater_than_zero <- function(col_from_dataset){
  m <- median(col_from_dataset) #e.g.: median(pcg_with_patientIdsfor_cancertype1[,1])
  if(m > 0){
    return(TRUE)
  }
  else{
    return(FALSE)
  }
}

#helper function 2
check_all_col <- function(dataset){
  col_name_greater_zero <- vector("list")
  ii = 1
  for(i in 1:19953){
    if(col_mean_is_greater_than_zero(dataset[,i])){
      col_name_greater_zero[ii] <- colnames(dataset)[i]
      ii = ii +1
    }
  }
  return(col_name_greater_zero)
}

#======================================================================================================================
setwd("/home/xwang/Documents/BCB430/data_oct2nd/")
#set working directory and load raw data 
lncRNAwithpatient <- readRDS("/home/xwang/Documents/BCB430/data_oct2nd/data/lnc_RNAseq_wPatientData_RankAg.rds")
proteincodingwithpatient <- readRDS("/home/xwang/Documents/BCB430/data_oct2nd/data/pcg_RNAseq_wPatientData_RankAg.rds")

#===========================  Task 1: Preprocess data for cancer type 2-7  =================================

for(i in 1:7){
  temp_n <- paste("canctype",i,sep="")
  filename <- paste(temp_n, "_lncRNAgeneIDs.Rda", sep="")
  load(filename)
  filename2 <- paste(temp_n,"_patientIDs.Rda", sep="")
  load(filename2)
  genelist_name <- paste(temp_n,".geneIDs",sep="")
  patlist_name <- paste(temp_n,".patientIDs",sep="")
  # 1) lncRNAs RNA-Seq data
  lncRNA_with_patientIdsfor_cancertype <- subset(lncRNAwithpatient, is.element(lncRNAwithpatient$patient, get(patlist_name)))
  highexpr_lncRNA_with_patientIds_cancertype <-subset(lncRNA_with_patientIdsfor_cancertype, select = is.element(colnames(lncRNA_with_patientIdsfor_cancertype),get(genelist_name)))
  #We have to log1p the lncRNA patient dataset.
  log1p_highexpr_lncRNA_with_patientIds_canctype <- as.data.frame(apply(highexpr_lncRNA_with_patientIds_cancertype,2, log1p))
  temp_save_f_name <- paste(temp_n,"_log1p_highexpr_lncRNAs_wpatient.Rda", sep="")
  save(log1p_highexpr_lncRNA_with_patientIds_canctype, file= temp_save_f_name)
  # 2) Protein-coding genes RNA-Seq data
  pcg_with_patientIdsfor_cancertype <- subset(proteincodingwithpatient, is.element(lncRNAwithpatient$patient,  get(patlist_name)))
  colnames_greater_zero <- check_all_col(pcg_with_patientIdsfor_cancertype)
  all_colnames <- colnames(pcg_with_patientIdsfor_cancertype)
  filter_pcg_cancertype <- subset(pcg_with_patientIdsfor_cancertype, select = is.element(all_colnames,colnames_greater_zero))
  # HAVE TO LOG1P the whole matrix.
  log1p_filtered_pcg_canctype <- as.data.frame(apply(filter_pcg_cancertype,2, log1p))
  temp_save_f_name2 <- paste(temp_n,"_log1p_filtered_pcg_wpatients.Rda",sep="")
  save(log1p_filtered_pcg_canctype, file = temp_save_f_name2)  
}

#====================================
#END