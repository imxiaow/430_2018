#===================  Defined functions =========================
# defined function, compute linear regression, anova, return file#1 (all pcg result)
find_cor_all_pcg_df <- function(lncRNA_col, pcg_df, store_signame, store_sigpv, store_sigcf){
  sig <- 0
  for (i in 1:length(pcg_df)){
    print(i)
    pcgname <- colnames(pcg_df)[i]
    pcgn_column <- subset(pcg_df, select = pcgname)
    lncRNA_pcg_data <- cbind(lncRNA_col, pcgn_column)
    check <- colnames(lncRNA_pcg_data)[2] #check pcg name
    check_value <- grepl("-",check)
    if(check_value){
      colnames(lncRNA_pcg_data)[2] <- gsub("-", "_", colnames(lncRNA_pcg_data)[2])
    }
    check2 <- colnames(lncRNA_pcg_data)[1] #check lncRNA name
    check_value2 <- grepl("-", check2)
    if(check_value2){ 
      colnames(lncRNA_pcg_data)[1]<- gsub("-","_",colnames(lncRNA_pcg_data)[1])
    }
    
    corre_test <- cor.test(unlist(lncRNA_pcg_data[1]),unlist(lncRNA_pcg_data[2]),method = "spearman")
    p_v <- corre_test$p.value
    corre <- cor(unlist(lncRNA_pcg_data[1]),unlist(lncRNA_pcg_data[2]), method="spearman")
    
    store_signame <- rbind(store_signame, pcgname)
    store_sigpv <- rbind(store_sigpv,p_v)
    store_sigcf <- rbind(store_sigcf, corre)
  }
  combined_cols <- cbind(store_signame, store_sigpv, store_sigcf)
  return(combined_cols)
}

#defined function, compute fdr by input data
fdr_get_sig <- function(result_data){
  raw_p <- result_data[,2]
  adjust_p <- p.adjust(raw_p, "fdr")
  comp_raw_adjust_p <- data.frame(raw=raw_p, fdr=adjust_p)
  adjust_sig_result_wrawfdr <- subset(comp_raw_adjust_p, comp_raw_adjust_p[,2] <= 0.05)
  adjust_sig_result <- subset(result_data, is.element(rownames(result_data),rownames(adjust_sig_result_wrawfdr)))
  adjust_sig_result <- cbind(adjust_sig_result, adjust_sig_result_wrawfdr[,2])
  colnames(adjust_sig_result)[4] <- "fdr"
  return(adjust_sig_result)
}

#=====================================================================
setwd("/home/xwang/Documents/BCB430/data_oct2nd/")


#===========================  Task 2: Compute  lncRNAs for input cancertype =================================
# run this script with command line arguments for each cancer type
args <- commandArgs(trailingOnly = TRUE)
cancertype <- args[1] # e.g. canctype1, canctype2,...(this is same as the filenames)

temp_filename1 <- paste(cancertype, "_log1p_highexpr_lncRNAs_wpatient.Rda",sep="")
temp_filename2 <- paste(cancertype, "_log1p_filtered_pcg_wpatients.Rda",sep="")

load(temp_filename1)
load(temp_filename2)

# all lncRNAs' name for this cancer type
lncRNA_names <- colnames(log1p_highexpr_lncRNA_with_patientIds_canctype) 

# save data in this directory
temp_dir_name <- paste("/home/xwang/Documents/BCB430/data_oct2nd/",cancertype,"/",sep ="")
setwd(temp_dir_name)

# get each lncRNA column one at a time
for(i in 1:length(lncRNA_names)){
  
  lncRNA_column <- subset(log1p_highexpr_lncRNA_with_patientIds_canctype, select = lncRNA_names[i])
  sig_pcg_names <- c()
  sig_pcg_pv <- c()
  sig_pcg_cf <- c()
 # sig_pcg_itcpt <- c()
  results <- find_cor_all_pcg_df(lncRNA_column, log1p_filtered_pcg_canctype, sig_pcg_names, sig_pcg_pv, sig_pcg_cf)
  colnames(results) <- cbind("pcg_names", "p-values", "coeffients")
  rownames(results) <- results[,1]
  temp_f1 <-paste("_raw_result_",cancertype,".Rda",sep="")
  filename1 <- paste(lncRNA_names[i], temp_f1,sep="") 
  save(results, file = filename1) # p-value without FDR correction 
  
  adjust_sig_result <- fdr_get_sig(results)
  temp_f2 <-paste("_fdrsig_result_", cancertype,".Rda",sep="")
  filename2 <-paste(lncRNA_names[i], temp_f2, sep="")
  save(adjust_sig_result, file = filename2)
  
  sig_pcg_patients <- subset(log1p_filtered_pcg_canctype,select=rownames(adjust_sig_result))
  temp_f3 <-paste("_fdrsig_wpatients_", cancertype, ".Rda", sep="")
  filename3 <- paste(lncRNA_names[i],temp_f3,sep="")
  save(sig_pcg_patients,file = filename3)
}

#==========
#END
