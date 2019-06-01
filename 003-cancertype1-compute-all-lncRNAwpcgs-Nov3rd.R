#=============================================================================================================
# Compute significantly correlated protein coding genes for 47 lncRNAs that are associated with cancertype1.
# Automate by running this script with loaded data, 3 files generated for each lncRNA.
# 1) result (all pcg, need this for fdr), contain coeffient, p-value,intercept, genename 
# 2) sig_result (after fdr), contain coeffient, p-value, fdr, intercept, genename
# 3) patient data with only significant pcgs. 
#
# Look at the results for this 47 lncRNAs,  What's unusual about the result?
# Make sure everything is "correct" before apply it to rest of the cancer types.
#
#==============================================================
setwd("/home/xwang/Documents/BCB430/data_oct2nd/")

load("log1p_cancertype1_highexpr_lncRNAs_wpatient.Rda") #47 lncRNA with patients df
load("cancertype1_log1p_filtered_pcg_wpatients.Rda") # 17808 pcgs with patietns df

#===================  Defined functions =========================

# defined function, compute linear regression, anova, return file#1 (all pcg result)
find_cor_all_pcg_df <- function(lncRNA_col, pcg_df, store_signame, store_sigpv, store_sigcf, store_itcpt){
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
    m1 <- lm(as.formula(paste(colnames(lncRNA_pcg_data)[2],"~ 1")), data = lncRNA_pcg_data)
    temp <- paste(" ~ ",colnames(lncRNA_pcg_data)[1])
    m2 <- lm(as.formula(paste(colnames(lncRNA_pcg_data)[2], temp)), data = lncRNA_pcg_data)
    cmpr <- anova(m1,m2)
    p_v <- cmpr$`Pr(>F)`[2]
    store_signame <- rbind(store_signame, pcgname)
    store_sigpv <- rbind(store_sigpv,p_v)
    store_sigcf <- rbind(store_sigcf, m2$coefficients[2])
    store_itcpt <- rbind(store_itcpt, m2$coefficients[1])
  }
  combined_cols <- cbind(store_signame, store_sigpv, store_sigcf, store_itcpt)
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
  colnames(adjust_sig_result)[5] <- "fdr"
  return(adjust_sig_result)
}

#===========================  Compute 47 lncRNAs  ============================================

t0 <- Sys.time() # count time 
# all lncRNAs' name
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
  results <- find_cor_all_pcg_df(lncRNA_column, log1p_filtered_pcg_canctype1, sig_pcg_names, sig_pcg_pv, sig_pcg_cf,sig_pcg_itcpt)
  colnames(results) <- cbind("pcg_names", "p-values", "m2-coeffients", "m2-intercept")
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

tf <- Sys.time() - t0  #get the time of this script from actuall computation/calculation
print(tf) # print at the end, estimate 2.5hrs, actual Time difference of 2.237182 hours

#========
#END