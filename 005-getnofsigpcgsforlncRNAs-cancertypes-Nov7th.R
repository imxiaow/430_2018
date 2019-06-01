#=====================================================================================================
# For each cancer type, get the number of PCGs for each lncRNA from that cancertype.
# So for each cancer type, there is a table with columns of "num_of_sig_PCGs", "lncRNA","cancer_type"
# Output would be seven tablee, we can combine them together if needed. 
#
# Next step would be generate table of lncRNA vs. PCGs of the same patients, for plotting. (10 for each)
#
# Task 1: read/load the file, count number of rows or columns, store this value.
# Task 2: create/ add this value, lncRNA name, cancertype to a table/df.
# Task 3: after go over all the files, save this df with corresponded name.
#
#======================================================================
setwd("/home/xwang/Documents/BCB430/data_oct2nd") # current directory

#=========================================
# Task 1, task2, task3
for(i in 1:7){
  setwd("/home/xwang/Documents/BCB430/data_oct2nd") # current directory
  foldername <- paste("canctype",i,sep="")
  pathnames_temp <- getwd()
  pathname <- paste(pathnames_temp,"/", foldername,sep="")
  setwd(pathname)
  
  files <- list.files()
  
  lncRNAname_c <- c()
  num_pcg_c <- c()
  canctype_c <- c()
  
  for(f in 1:length(files)){
    check_value <- grepl("fdrsig_result", files[f])
    if(check_value){
      canctype <- foldername
      lncRNA_name <- strsplit(files[f],"_")[[1]][1]
      load(files[f])
      num_pcg <- nrow(adjust_sig_result)
      
      lncRNAname_c <- rbind(lncRNAname_c, lncRNA_name)
      num_pcg_c <- rbind(num_pcg_c, num_pcg)
      canctype_c <- rbind(canctype_c, canctype)
    }
  }
  combined_c <- cbind(lncRNAname_c, num_pcg_c,canctype_c)
  colnames(combined_c) <- cbind("lncRNA", "num_PCGs", "Canctype")
  rownames(combined_c) <- combined_c[,1]
  filename_table <- paste(foldername,"_Summary_lncRNAs_num_pcg.Rda")
  setwd("/home/xwang/Documents/BCB430/data_oct2nd/Summary/") 
  save(combined_c, file=filename_table)
}


#=====================================================================
setwd("/home/xwang/Documents/BCB430/data_oct2nd/Summary/") 
#combine all the df together 
fs <- list.files()

all_summary <- c()
for(fl in 1:length(fs)){
  check_f <- grepl("_Summary_", fs[fl])
  if(check_f){
    load(fs[fl])
    all_summary <- rbind(all_summary,combined_c)
  }
}
save(all_summary, file="all_canctype_Summary_lncRNAs_num_pcg.Rda")

#===========================
#END