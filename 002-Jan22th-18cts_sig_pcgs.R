#=====================================================
# 3. for each cancer type, find significant correlated pcgs(fdr <0.05), output: 18 pcg-lists for each lncRNAs.
# 4. RA, combine them together, filter genes with RA-score <0.05
# 5. compare the result with PCAWGs.
#========================================
# load needed dataset, saved from previous script.
# c1_pcgs ~ c18_pcgs          18
# c1_MALAT1 ~ c18_MALAT1       18
# c1_NEAT1 ~ c18_NEAT1        18 


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
#===============================================================



#====================================
#END