load("canctype1_lncRNAgeneIDs.Rda")
load("canctype1_patientIDs.Rda")
load("cancertype1_pcg_wpatients.Rda")
load("log1p_cancertype1_highexpr_lncRNAs_wpatient.Rda")


#===============  continue from 001Prepdata_Oct16-23.R   ===========================================================================
# For cancertype1, select 1 highly expressed lncRNA (47 genes total), filter pcg (80 x 19956)
# Find the correlation of this lncRNA with filtered pcgs, plot it (80 points) for (# of pcg) graphs, 
# Identity potential target genes that are higly positive or negative. 
# familiar with the pipline, so that could apply to all other lncRNAs in this cancertype and other cancertypes.

#=============================  Task 1: filter pcg dataset  ========================================================================
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

colnames_greater_zero <- check_all_col(pcg_with_patientIdsfor_cancertype1)
all_colnames <- colnames(pcg_with_patientIdsfor_cancertype1)

filter_pcg_cancertype1 <- subset(pcg_with_patientIdsfor_cancertype1, select = is.element(all_colnames,colnames_greater_zero))
#is.element function -> put less number of items at last
save(filter_pcg_cancertype1, file = "cancertype1_filtered_pcg_wpatients.Rda")

# HAVE TO LOG1P the whole matrix.
log1p_filtered_pcg_canctype1 <- as.data.frame(apply(filter_pcg_cancertype1,2, log1p))
save(log1p_filtered_pcg_canctype1, file = "cancertype1_log1p_filtered_pcg_wpatients.Rda")

  
#============================  Task 2: find correlation of one lncRNA with all pcg expression =======================================
load("cancertype1_log1p_filtered_pcg_wpatients.Rda")

first_lncRNA <- as.character(canctype1.geneIDs[5])


# linear regression in R
# need 1) coefficient; 2) the result of comparsion to between 2 models (model1:null model(PCG1~1), model2:lm(PCG1~lncRNAE))
# write function and apply to each column
# 1) lm(pcg1~1) #null model
# 2) lm(pcg1 ~lncRNAE) -> coefficient
# anova(model1, model2), if p-value < 0.05 => significant correlated

# step 1: take the lncRNA coloumn (1x80) for all patients relate to the same cancer type e.g. cancer type 1
# step 2: for each pcg columns, 2 models -> one for null, one is the correlation of lncRNA and pcg. get the significance. (p value)
# step 3: apply this to all of the columns in log1p_pcg dataset
# step 4: filter the non-significant correlations, plot the significant correlated lncRNA with pcg graphs for all patients.


# step 1
lncRNA_column <- subset(log1p_highexpr_lncRNA_with_patientIds_canctype1, select = first_lncRNA)

# step 2
pcg1 <- colnames(log1p_filtered_pcg_canctype1)[1]
pcg1_column <- subset(log1p_filtered_pcg_canctype1, select = pcg1)
#lm(as.formula(paste(pcg1, "~1")), pcg1_column) # the way to use as.formula 

#combine lncRNA column with pcg colum in order to use lm to find the correlation between them
lncRNA_pcg_dat <- cbind(lncRNA_column, pcg1_column) # values are filtered

m1 <- lm(as.formula(paste(colnames(lncRNA_pcg_dat)[2], "~1")), lncRNA_pcg_dat) #null model
temp <- paste("~", colnames(lncRNA_pcg_dat)[1])
m2 <- lm(as.formula(paste(colnames(lncRNA_pcg_dat)[2], temp)), lncRNA_pcg_dat) #linear regression of 2 genes.

cmpr <- anova(m1,m2)
p_v <- cmpr$`Pr(>F)`[2] # the P-value, If P-value is <0.05, then these 2 genes are significantly correlated.


#step 3: apply this to all columns in pcg df
find_cor_all_pcg_df <- function(lncRNA_col, pcg_df, store_signame, store_sigpv, store_sigcf, store_itcpt){
  sig <- 0
  for (i in 1:length(pcg_df)){
    print(i)
    pcgname <- colnames(pcg_df)[i]
    pcgn_column <- subset(pcg_df, select = pcgname)
    lncRNA_pcg_data <- cbind(lncRNA_col, pcgn_column)
    check <- colnames(lncRNA_pcg_data)[2]
    check_value <- grepl("-",check)
    if(check_value){ # check pcg names
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
  
   # sig <- sig + 1  #count number of significant
    store_signame <- rbind(store_signame, pcgname)
    store_sigpv <- rbind(store_sigpv,p_v)
    store_sigcf <- rbind(store_sigcf, m2$coefficients[2])
    store_itcpt <- rbind(store_itcpt, m2$coefficients[1])
  }
 # print(sig) #1037
  combined_cols <- cbind(store_signame, store_sigpv, store_sigcf, store_itcpt)
  return(combined_cols)
}

sig_pcg_names <- c()
sig_pcg_pv <- c()
sig_pcg_cf <- c()
sig_pcg_itcpt <- c()

t0 <- Sys.time() # count time used for one lncRNA 
results <- find_cor_all_pcg_df(lncRNA_column, log1p_filtered_pcg_canctype1, sig_pcg_names, sig_pcg_pv, sig_pcg_cf,sig_pcg_itcpt)
tf <- Sys.time() - t0 # time go over the pcg df: 3.45777 mins

colnames(results) <- cbind("pcg_names", "p-values", "m2-coeffients", "m2-intercept")
rownames(results) <- results[,1]

save(results, file="first_lncRNA_sig_results_allpcg.Rda") # p-value without FDR correction 
#load this data
#load("first_lncRNA_sig_results_allpcg.Rda")

raw_p <- results[,2]
adjust_p <- p.adjust(raw_p, "fdr")
comp_raw_adjust_p <- data.frame(raw=raw_p, fdr=adjust_p)
adjust_sig_result_wrawfdr <- subset(comp_raw_adjust_p, comp_raw_adjust_p[,2] <= 0.05)
adjust_sig_result <- subset(results, is.element(rownames(results),rownames(adjust_sig_result_wrawfdr)))
adjust_sig_result <- cbind(adjust_sig_result, adjust_sig_result_wrawfdr[,2])
colnames(adjust_sig_result)[5] <- "fdr"

save(adjust_sig_result, file = "first_lncRNA_sig_results_adjusted_wfdr.Rda")
#load this data
#load("first_lncRNA_sig_results_adjusted_wfdr.Rda")

#in order to plot the significant lncRNA with PCGs, we need the coeffient and intercept from lm(linear regression)
# also, we need the patients data (point) for the plot,
# we have to get the patient data for significant pcg (by subset it)
sig_pcg_patients <- subset(log1p_filtered_pcg_canctype1,select=rownames(adjust_sig_result))

save(sig_pcg_patients,file="first_lncRNA_sig_result_wPatients_data.Rda")
#load this after by
#load("first_lncRNA_sig_result_wPatients_data.Rda")


#===============================   Task 3: Plotting  ==================================================================================
# ggplot 2 ...
if(!require(ggplot2)) {
  install.packages("ggplot2")
  library(ggplot2)
}

#plot()
#abline() add straight Lines to plot
#geom_line 

#========================================================================================================
# pipline for one lncRNA and all pcgs, apply this pipeline to all lncRNAs (highly expressed) associated
# with cancer type 1. 
# Question to think about: how to store the sig data and plots for each lncRNAs separetely while running the program?
# How to name the files...

# Then, apply to all cacer types. Manage the data generated. 
#END