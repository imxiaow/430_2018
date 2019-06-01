setwd("/home/xwang/Documents/BCB430/data_oct2nd/")

load("log1p_cancertype1_highexpr_lncRNAs_wpatient.Rda") #47 lncRNA with patients df
load("cancertype1_log1p_filtered_pcg_wpatients.Rda") # 17808 pcgs with patietns df

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

#====================== Test for cancer type 1 =======================================================
lncRNA_names <- colnames(log1p_highexpr_lncRNA_with_patientIds_canctype1) 

# save data in this directory
setwd("/home/xwang/Documents/BCB430/data_oct2nd/canctype1/")

# get each lncRNA column one at a time
for(i in 1:length(lncRNA_names)){
  
  lncRNA_column <- subset(log1p_highexpr_lncRNA_with_patientIds_canctype1, select = lncRNA_names[i])
  sig_pcg_names <- c()
  sig_pcg_pv <- c()
  sig_pcg_cf <- c()
  sig_pcg_itcpt <- c()
  results <- find_cor_all_pcg_df(lncRNA_column, log1p_filtered_pcg_canctype1, sig_pcg_names, sig_pcg_pv, sig_pcg_cf)
  colnames(results) <- cbind("pcg_names", "p-values", "coeffients")
  rownames(results) <- results[,1]
  filename1 <- paste(lncRNA_names[i],"_raw_result_canctype1.Rda" ,sep="") 
  save(results, file = filename1) # p-value without FDR correction 
  
  adjust_sig_result <- fdr_get_sig(results)
  filename2 <-paste(lncRNA_names[i],"_fdrsig_result_canctype1.Rda" ,sep="")
  save(adjust_sig_result, file = filename2)
  
  sig_pcg_patients <- subset(log1p_filtered_pcg_canctype1,select=rownames(adjust_sig_result))
  filename3 <- paste(lncRNA_names[i],"_fdrsig_wpatients_canctype1.Rda",sep="")
  save(sig_pcg_patients,file = filename3)
}
