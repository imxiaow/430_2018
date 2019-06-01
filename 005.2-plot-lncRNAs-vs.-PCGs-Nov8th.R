#==================================================================================================================
# Continue from previous script, the purpose of this script would be:
# 1) for each lncRNAs, get the table (~10 for each lncRNA) of lncRNA vs. pcg of same number of patients.(e.g.80)
# 2) Once we have this table, we can plot this table.
# 3) Save plots as PDF. 
# 4) Also, use plot to visilize the distribution of num_sig_pcgs for all the lncRNAs under different canctype. 
# 
#==============================================================================================
#setwd("/home/xwang/Documents/BCB430/data_oct2nd/")
# 1)
# have go to each cancer type folders by reset the directory everytime. 
# need "_fdrsig_wpatients_"canctypeX.Rda  -> load this file, variable name: sig_pcg_patients,
#                                             pcgs are column wise, get it each time,
# 2)
# lncRNA RNA-seq data are located in the original directory, 
# "canctypeX_log1p_highexpr_lncRNAs_wpatient.Rda" -> load this file, 
# variable name: log1p_highexpr_lncRNA_with_patientIds_canctype
# lncRNas are column wise, get it each time, match with sig-pcg list filename
# 
# 3)
# num_sig_pcg summary for each cancer type, and the summary of all cancer types are located in Summary folder.
# change the directory to Summary folder, and plot 8 tables. Save them.
#
#==================================================================================
if (!require(ggplot2)) {
  install.packages("ggplot2")
  library(ggplot2) }
#======================================================
# Task 1: plot summary tables
setwd("/home/xwang/Documents/BCB430/data_oct2nd/Summary/")

fs <- list.files()
pdf("summary_lncRNAs_numpcgs.pdf")
for(i in 1:length(fs)){
  check_f <- grepl("canctype", fs[i])
  if(check_f){
    load(fs[i])
    data <- as.data.frame(combined_c)
    data$lncRNA <- as.character(data$lncRNA)
    data$num_PCGs <- as.numeric(as.character(data$num_PCGs))
    gtilename <- strsplit(fs[i], "_")[[1]][1]
    print(ggplot(data, aes(x=lncRNA, y=num_PCGs,fill= Canctype)) + geom_bar(stat='identity') + geom_text(aes(label=num_PCGs), vjust=1.6,color="Black",size=1.5, position = position_identity()) + theme(legend.position = "right") + ggtitle(gtilename) + theme(axis.text.x = element_text(size =3, angle = 90, hjust =1)))
  }
  check_f2 <- grepl("cancertype",fs[i])
  if(check_f2){
    load(fs[i])
    data <- as.data.frame(all_summary)
    data$lncRNA <- as.character(data$lncRNA)
    data$num_PCGs <- as.numeric(as.character(data$num_PCGs))
    gtilename <- strsplit(fs[i], "_")[[1]][1]
    print(ggplot(data, aes(x=lncRNA, y=num_PCGs,fill= Canctype)) + geom_bar(stat='identity') + theme(legend.position = "right") + ggtitle(gtilename) + theme(axis.text.x = element_text(size =2,angle = 90, hjust =1)))
  }
}
dev.off()
#======================================================


#==================
#END