#=================================================================
# get the genelist for each cancer type for NEAT1, MALAT1.
# Use as comparsion from simple cor result, and also use for RankAggreg.
# Can eailsy find out the overlap pcgs from 2 methods?? 
# compare the plotting from both methods, some of the lncRNAs have more sig pcgs, some have less. 
# rank the list with rankAggreg method, Question: the gene list for aggregateRanks() only need gene names?
#             - How did they generate the 

# Task:
# For each foler, go into its folder directory, find same "gene_name_fdrsig_result_canctypeX.Rda" for 7 canc.
# the variable name called "adjust_sig_result", store it, go back to home dir, 1-> 7.
#====================================================
if (!require(RobustRankAggreg)) {
  install.packages("RobustRankAggreg")
  library(RobustRankAggreg) }

#====================== NEAT1 ====================================
for(i in 1:7){
  setwd("/home/xwang/Documents/BCB430/data_oct2nd/cor_results") # current directory
  foldername <- paste("cor_canctype",i,sep="")
  pathnames_temp <- getwd()
  pathname <- paste(pathnames_temp,"/", foldername,sep="")
  setwd(pathname)
  
  files <- list.files()

  for(f in 1:length(files)){
    #check_value <- grepl("fdrsig_result", files[f])
    check_gene_name <- grepl("NEAT1_fdrsig_result",files[f])
    if(check_gene_name){
      canctype <- foldername
      load(files[f])
      pcgs_lst <- names(sort(adjust_sig_result[,4]))
      fn <- paste("NEAT1_", foldername, "_pcglst.Rda")
      setwd("/home/xwang/Documents/BCB430/data_oct2nd/cor_results/NEAT1")
      save(pcgs_lst, file=fn)
    }
  }
}

#==== for NEAT1 sig_result check ==== 
for(i in 1:7){
  setwd("~/Documents/BCB430/data_oct2nd/cor_results") # current directory
  foldername <- paste("cor_canctype",i,sep="")
  pathnames_temp <- getwd()
  pathname <- paste(pathnames_temp,"/", foldername,sep="")
  setwd(pathname)
  
  files <- list.files()
  
  for(f in 1:length(files)){
    #check_value <- grepl("fdrsig_result", files[f])
    check_gene_name <- grepl("NEAT1_fdrsig_result",files[f])
    if(check_gene_name){
      canctype <- foldername
      load(files[f])
      pcgs_lst <- names(adjust_sig_result[,3])
      cor_result <- (as.character(adjust_sig_result[,3])) # the coeffient, cor result
      comb_list <- cbind(pcgs_lst, cor_result)
      # fn <- paste("NEAT1_", foldername, "_pcglst.Rda")
      fn <- paste("NEAT1_", foldername, "_pcgwithcorsig.Rda")
      setwd("~/Documents/BCB430/data_oct2nd/cor_results/NEAT1")
      save(comb_list, file=fn)
    }
  }
}
# combine them into a list 
setwd("/home/xwang/Documents/BCB430/data_oct2nd/cor_results/NEAT1")

lncRNAname <-"NEAT1"
files <- list.files()

for(f in 1:length(files)){
  if(files[f] == "NEAT1_ cor_canctype1 _pcgwithcorsig.Rda"){
    NEAT1_pcg_cor_lst1 <- get(load(files[f]))
  }

  if(files[f] == "NEAT1_ cor_canctype2 _pcgwithcorsig.Rda"){
    NEAT1_pcg_cor_lst2 <- get(load(files[f]))
  }
  
  if(files[f] == "NEAT1_ cor_canctype3 _pcgwithcorsig.Rda"){
    NEAT1_pcg_cor_lst3 <- get(load(files[f]))
  }
  
  if(files[f] == "NEAT1_ cor_canctype4 _pcgwithcorsig.Rda"){
    NEAT1_pcg_cor_lst4 <- get(load(files[f]))
  }
  
  if(files[f] == "NEAT1_ cor_canctype5 _pcgwithcorsig.Rda"){
    NEAT1_pcg_cor_lst5 <- get(load(files[f]))
  }
  
  if(files[f] == "NEAT1_ cor_canctype6 _pcgwithcorsig.Rda"){
    NEAT1_pcg_cor_lst6 <- get(load(files[f]))
  }
  
  if(files[f] == "NEAT1_ cor_canctype7 _pcgwithcorsig.Rda"){
    NEAT1_pcg_cor_lst7 <- get(load(files[f]))
  }
}
#check if one list have positive value while other list have negative value. y

check_gene_inotherlist <- function(ls1,ls2){
  #if there's gene also appear in other list, 
  # check the value
  if (length(intersect(ls1,ls2)) != 0){
   # print(intersect(ls1,ls2))
    common_g <- intersect(ls1,ls2)
    cat("Total common gene #: ", length(common_g), "\n")
    test1 <- as.data.frame(ls1, stringsAsFactors = FALSE)
    test1[,2] <- as.numeric(test1[,2])
    
    test2 <- as.data.frame(ls2,stringsAsFactors = FALSE)
    test2[,2] <- as.numeric(test2[,2])
    
    WRONG <- c()
    for(i in 1:length(common_g)){
      ls1_value <- test1$cor_result[which(test1$pcgs_lst == common_g[i])]
      ls2_value <- test2$cor_result[which(test2$pcgs_lst == common_g[i])]
      if(is.null(ls1_value) | is.na(ls1_value)| ls1_value==""){
        print(common_g[i])
        WRONG <- c(WRONG, common_g[i])
      }
      else{
        if(ls1_value >=0){
          if(ls2_value >=0){
            #pass
          }
          else{
            print(common_g[i])
            WRONG <- c(WRONG, common_g[i])
          }
        }
        
        if(ls1_value <0){
          if(ls2_value <0){
            #pass
          }
          else{
            print(common_g[i])
            WRONG <- c(WRONG, common_g[i])
          }
        }
      }

    }
      cat("Wrong #: ",length(WRONG))
  }
}







#====================== MALAT1 ====================================
for(i in 1:7){
  setwd("/home/xwang/Documents/BCB430/data_oct2nd/cor_results/") # current directory
  foldername <- paste("cor_canctype",i,sep="")
  pathnames_temp <- getwd()
  pathname <- paste(pathnames_temp,"/", foldername,sep="")
  setwd(pathname)
  
  files <- list.files()
  
  for(f in 1:length(files)){
    #check_value <- grepl("fdrsig_result", files[f])
    check_gene_name <- grepl("MALAT1_fdrsig_result",files[f])
    if(check_gene_name){
      canctype <- foldername
      load(files[f])
      pcgs_lst <- names(sort(adjust_sig_result[,4]))
      fn <- paste("MALAT1_", foldername, "_pcglst.Rda")
      setwd("/home/xwang/Documents/BCB430/data_oct2nd/cor_results/MALAT1")
      save(pcgs_lst, file=fn)
    }
  }
}

#==================== MALAT1 ===================
#set variable name for each Rda, make them into one glist
setwd("/home/xwang/Documents/BCB430/data_oct2nd/cor_results/MALAT1")

lncRNAname <-"MALAT1"
files <- list.files()

for(f in 1:length(files)){
  cantype <- strsplit(files[f],"_")[[1]][3]
  if(cantype == "canctype1 "){
    MALAT1_pcg_lst1 <- get(load(files[f]))

  }
  if(cantype == "canctype2 "){
    MALAT1_pcg_lst2 <- get(load(files[f]))

  }
  if(cantype == "canctype3 "){
    MALAT1_pcg_lst3 <- get(load(files[f]))

  }
  if(cantype == "canctype4 "){
    MALAT1_pcg_lst4 <- get(load(files[f]))

  }
  if(cantype == "canctype5 "){
    MALAT1_pcg_lst5 <- get(load(files[f]))
  }
  if(cantype == "canctype6 "){
    MALAT1_pcg_lst6 <- get(load(files[f]))
  }
  if(cantype == "canctype7 "){
    MALAT1_pcg_lst7 <- get(load(files[f]))
  }
}

MALAT1_glist <- list(MALAT1_pcg_lst1, MALAT1_pcg_lst2, MALAT1_pcg_lst3, MALAT1_pcg_lst4, MALAT1_pcg_lst5, MALAT1_pcg_lst6,MALAT1_pcg_lst7)
ra_MALAT1 <- aggregateRanks(glist=MALAT1_glist)
save(ra_MALAT1, file="MALAT1_RA_result.Rda")

#================ NEAT1 ========================
setwd("/home/xwang/Documents/BCB430/data_oct2nd/cor_results/NEAT1")

lncRNAname <-"NEAT1"
files <- list.files()

for(f in 1:length(files)){
  cantype <- strsplit(files[f],"_")[[1]][3]
  if(cantype == "canctype1 "){
    NEAT1_pcg_lst1 <- get(load(files[f]))
    
  }
  if(cantype == "canctype2 "){
    NEAT1_pcg_lst2 <- get(load(files[f]))
    
  }
  if(cantype == "canctype3 "){
    NEAT1_pcg_lst3 <- get(load(files[f]))
    
  }
  if(cantype == "canctype4 "){
    NEAT1_pcg_lst4 <- get(load(files[f]))
    
  }
  if(cantype == "canctype5 "){
    NEAT1_pcg_lst5 <- get(load(files[f]))
  }
  if(cantype == "canctype6 "){
    NEAT1_pcg_lst6 <- get(load(files[f]))
  }
  if(cantype == "canctype7 "){
    NEAT1_pcg_lst7 <- get(load(files[f]))
  }
}
NEAT1_glist <- list(NEAT1_pcg_lst1, NEAT1_pcg_lst2, NEAT1_pcg_lst3, NEAT1_pcg_lst4, NEAT1_pcg_lst5, NEAT1_pcg_lst6, NEAT1_pcg_lst7)
ra_NEAT1 <- aggregateRanks(glist=NEAT1_glist)
save(ra_NEAT1, file="NEAT1_RA_result.Rda")

#===========================
# question: is the glist already ranked? rank by p-value? 

# After meeting with juri, 1) rank the gene list by fdr,  before put into rank Aggreg, 
#                          2) for 7 lists, compute genes by group, e.g. 7/7, 6/7.,..1/7. (function)
                    #      3) R.A. them, note the number of returned genelist. 
                        #  4) what's the R.A. score of each group, 
                        #  5) boxplot them to better visualize. 
# 1) group lists
# 2) get R.A. score for each group's gene. 
# 3) box-plot 


#==============  Defined function  =======================================================

group_genes_by_apearance <- function(lst1, lst2, lst3, lst4, lst5, lst6, lst7, count_number){
  # combine 7 cancer type list together to count times of appearance.
  combine_lst <- c(lst1, lst2, lst3, lst4, lst5, lst6, lst7)
  print(combine_lst)
  # create new list by times of appearance
  seen1 <- c()
  seen2 <- c()
  seen3 <- c()
  seen4 <- c()
  seen5 <- c()
  seen6 <- c()
  seen7 <- c()
  
  for(item in 1:length(combine_lst)){
    # if this item is never seen before, 
    # add it to its list
    if (! combine_lst[item]%in%seen1){
      seen1 <- c(seen1, combine_lst[item])
    }
    else{
      
      if(! combine_lst[item]%in%seen2){
        seen2 <- c(seen2, combine_lst[item])
      }
      else{
        
        if(! combine_lst[item]%in%seen3){
          seen3 <- c(seen3, combine_lst[item])
        }
        else{
          
          if(! combine_lst[item]%in%seen4){
            seen4 <- c(seen4, combine_lst[item])
          }
          else{
            
            if(! combine_lst[item]%in%seen5){
              seen5 <- c(seen5, combine_lst[item])
            }
            else{
              
              if(! combine_lst[item]%in%seen6){
                seen6 <- c(seen6, combine_lst[item])
              }
              else{
                
                if(! combine_lst[item]%in%seen7){
                  seen7 <- c(seen7, combine_lst[item])
                }
              }
            }
          }
        }
      }
    }
  }
  # after this for loop, item that appear in 7 times will in list seen6.
  print(length(seen7))
  print(length(seen6))
  print(length(seen5))
  print(length(seen4))
  print(length(seen3))
  print(length(seen2))
  print(length(seen1))
  
  if (count_number == 7){
    return(seen7)
  }
  if (count_number == 6){
    # for seen6, get item that not appear in seen7, 
    new_lst <- c()
    for (i in 1:length(seen6)){
      if(! seen6[i]%in%seen7){
        new_lst <- c(new_lst, seen6[i])
      }
    }
    return(new_lst)
  }
  # for seen5, get item that not appear in seen6,
  if (count_number == 5){
    new_lst <- c()
    for (i in 1:length(seen5)){
      if(! seen5[i]%in%seen6){
        new_lst <- c(new_lst, seen5[i])
      }
    }
    return(new_lst)
  }
  # for seen4, get item that not appear in seen5, 
  if (count_number == 4){
    new_lst <- c()
    for (i in 1:length(seen4)){
      if(! seen4[i]%in%seen5){
        new_lst <- c(new_lst, seen4[i])
      }
    }
    return(new_lst)
  }
  # for seen3, get item that not appear in seen4, 
  if (count_number == 3){
    new_lst <- c()
    for (i in 1:length(seen3)){
      if(! seen3[i]%in%seen4){
        new_lst <- c(new_lst, seen3[i])
      }
    }
    return(new_lst)
  }
  # for seen2, get item that not appear in seen3, 
  if (count_number == 2){
    new_lst <- c()
    for (i in 1:length(seen2)){
      if(! seen2[i]%in%seen3){
        new_lst <- c(new_lst, seen2[i])
      }
    }
    return(new_lst)
  }
  # for seen1, get item that not appear in seen2. 
  if (count_number == 1){
    new_lst <- c()
    for (i in 1:length(seen1)){
      if(! seen1[i]%in%seen2){
        new_lst <- c(new_lst, seen1[i])
      }
    }
    return(new_lst)
  }
  
}

#========================================================
# use previously defined function, to group genes by appearance.
# MALAT1
 
# 0 for 7/7
MALAT1_appear7 <- group_genes_by_apearance(MALAT1_pcg_lst1, MALAT1_pcg_lst2, MALAT1_pcg_lst3, 
                                           MALAT1_pcg_lst4, MALAT1_pcg_lst5, MALAT1_pcg_lst6, 
                                          MALAT1_pcg_lst7, 7)
# MALAT1_appear7 list is NULL

# printout message:
# 0       -> appear 7
# 34      -> appear 6
# 270     -> appear 5
# 1386    -> appear 4
# 3434    -> appear 3
# 7080    -> appear 2
# 12748   -> appear 1

# 34 for 6/7
MALAT1_appear6 <- group_genes_by_apearance(MALAT1_pcg_lst1, MALAT1_pcg_lst2, MALAT1_pcg_lst3, 
                                           MALAT1_pcg_lst4, MALAT1_pcg_lst5, MALAT1_pcg_lst6, 
                                           MALAT1_pcg_lst7, 6)
# 236 for 5/7
MALAT1_appear5 <- group_genes_by_apearance(MALAT1_pcg_lst1, MALAT1_pcg_lst2, MALAT1_pcg_lst3, 
                                           MALAT1_pcg_lst4, MALAT1_pcg_lst5, MALAT1_pcg_lst6, 
                                           MALAT1_pcg_lst7, 5)
# 1116 for 4/7
MALAT1_appear4 <- group_genes_by_apearance(MALAT1_pcg_lst1, MALAT1_pcg_lst2, MALAT1_pcg_lst3, 
                                           MALAT1_pcg_lst4, MALAT1_pcg_lst5, MALAT1_pcg_lst6, 
                                           MALAT1_pcg_lst7, 4)
# 2048 for 3/7
MALAT1_appear3 <- group_genes_by_apearance(MALAT1_pcg_lst1, MALAT1_pcg_lst2, MALAT1_pcg_lst3, 
                                           MALAT1_pcg_lst4, MALAT1_pcg_lst5, MALAT1_pcg_lst6, 
                                           MALAT1_pcg_lst7, 3)
# 3646 for 2/7
MALAT1_appear2 <- group_genes_by_apearance(MALAT1_pcg_lst1, MALAT1_pcg_lst2, MALAT1_pcg_lst3, 
                                           MALAT1_pcg_lst4, MALAT1_pcg_lst5, MALAT1_pcg_lst6, 
                                           MALAT1_pcg_lst7, 2)
# 5668 for 1/7
MALAT1_appear1 <- group_genes_by_apearance(MALAT1_pcg_lst1, MALAT1_pcg_lst2, MALAT1_pcg_lst3, 
                                           MALAT1_pcg_lst4, MALAT1_pcg_lst5, MALAT1_pcg_lst6, 
                                           MALAT1_pcg_lst7, 1)

# save the data into MALAT1 directory
setwd("/home/xwang/Documents/BCB430/data_oct2nd/cor_results/MALAT1")
save(MALAT1_appear1, file="MALAT1_appear1.Rda")
save(MALAT1_appear2, file="MALAT1_appear2.Rda")
save(MALAT1_appear3, file="MALAT1_appear3.Rda")
save(MALAT1_appear4, file="MALAT1_appear4.Rda")
save(MALAT1_appear5, file="MALAT1_appear5.Rda")
save(MALAT1_appear6, file="MALAT1_appear6.Rda")

#==============================================================
# NEAT1

# 15 for 7/7
NEAT1_appear7 <- group_genes_by_apearance(NEAT1_pcg_lst1, NEAT1_pcg_lst2, NEAT1_pcg_lst3, 
                                          NEAT1_pcg_lst4, NEAT1_pcg_lst5, NEAT1_pcg_lst6, 
                                          NEAT1_pcg_lst7, 7)
# printout message:
#[1] 15
#[1] 106
#[1] 425
#[1] 1174
#[1] 2570
#[1] 5148
#[1] 10599

# 91 for 6/7
NEAT1_appear6 <- group_genes_by_apearance(NEAT1_pcg_lst1, NEAT1_pcg_lst2, NEAT1_pcg_lst3, 
                                          NEAT1_pcg_lst4, NEAT1_pcg_lst5, NEAT1_pcg_lst6, 
                                          NEAT1_pcg_lst7, 6)
# 319 for 5/7
NEAT1_appear5 <- group_genes_by_apearance(NEAT1_pcg_lst1, NEAT1_pcg_lst2, NEAT1_pcg_lst3, 
                                          NEAT1_pcg_lst4, NEAT1_pcg_lst5, NEAT1_pcg_lst6, 
                                          NEAT1_pcg_lst7, 5)
# 749 for 4/7
NEAT1_appear4 <- group_genes_by_apearance(NEAT1_pcg_lst1, NEAT1_pcg_lst2, NEAT1_pcg_lst3, 
                                          NEAT1_pcg_lst4, NEAT1_pcg_lst5, NEAT1_pcg_lst6, 
                                          NEAT1_pcg_lst7, 4)
# 1396 for 3/7
NEAT1_appear3 <- group_genes_by_apearance(NEAT1_pcg_lst1, NEAT1_pcg_lst2, NEAT1_pcg_lst3, 
                                          NEAT1_pcg_lst4, NEAT1_pcg_lst5, NEAT1_pcg_lst6, 
                                          NEAT1_pcg_lst7, 3)
# 2578 for 2/7
NEAT1_appear2 <- group_genes_by_apearance(NEAT1_pcg_lst1, NEAT1_pcg_lst2, NEAT1_pcg_lst3, 
                                          NEAT1_pcg_lst4, NEAT1_pcg_lst5, NEAT1_pcg_lst6, 
                                          NEAT1_pcg_lst7, 2)
# 5451 for 1/7
NEAT1_appear1 <- group_genes_by_apearance(NEAT1_pcg_lst1, NEAT1_pcg_lst2, NEAT1_pcg_lst3, 
                                          NEAT1_pcg_lst4, NEAT1_pcg_lst5, NEAT1_pcg_lst6, 
                                          NEAT1_pcg_lst7, 1)

# save the data into MALAT1 directory
setwd("/home/xwang/Documents/BCB430/data_oct2nd/cor_results/NEAT1")
save(NEAT1_appear1, file="NEAT1_appear1.Rda")
save(NEAT1_appear2, file="NEAT1_appear2.Rda")
save(NEAT1_appear3, file="NEAT1_appear3.Rda")
save(NEAT1_appear4, file="NEAT1_appear4.Rda")
save(NEAT1_appear5, file="NEAT1_appear5.Rda")
save(NEAT1_appear6, file="NEAT1_appear6.Rda")
save(NEAT1_appear7, file="NEAT1_appear7.Rda")


#====================================================================
# right now, I have grouped genes, RA scores, 
# create new data frame, for each group, get RA score for each gene in that group.
# have 7 data frame, each contain genename and its RA score. 
# plot them into one box-plot. 


#==============
#END