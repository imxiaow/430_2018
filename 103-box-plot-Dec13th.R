#============================================================================
# load data for box-plot 
setwd("/home/xwang/Documents/BCB430/data_oct2nd/cor_results/NEAT1")
load("NEAT1_appear1.Rda")
load("NEAT1_appear2.Rda")
load("NEAT1_appear3.Rda")
load("NEAT1_appear4.Rda")
load("NEAT1_appear5.Rda")
load("NEAT1_appear6.Rda")
load("NEAT1_appear7.Rda")

load("NEAT1_RA_result.Rda")

setwd("/home/xwang/Documents/BCB430/data_oct2nd/cor_results/MALAT1")
load("MALAT1_appear1.Rda")
load("MALAT1_appear2.Rda")
load("MALAT1_appear3.Rda")
load("MALAT1_appear4.Rda")
load("MALAT1_appear5.Rda")
load("MALAT1_appear6.Rda")

load("MALAT1_RA_result.Rda")

#================================================

if (!require(ggplot2)) {
  install.packages("ggplot2")
  library(ggplot) }

#====================  Defined function =======================================

get_group_ID <- function(col_lst_character, lst1, lst2, lst3, lst4, lst5, lst6, lst7){
  group_ID <- c()
  query_lst <- as.character(col_lst_character)
  
  for(i in 1:length(query_lst)){
    if(query_lst[i]%in% lst1){
      group_ID <- c(group_ID, "1/7")
    }
    if(query_lst[i]%in% lst2){
      group_ID <- c(group_ID, "2/7")
    }
    if(query_lst[i]%in% lst3){
      group_ID <- c(group_ID, "3/7")
    }
    if(query_lst[i]%in% lst4){
      group_ID <- c(group_ID, "4/7")
    }
    if(query_lst[i]%in% lst5){
      group_ID <- c(group_ID, "5/7")
    }
    if(query_lst[i]%in% lst6){
      group_ID <- c(group_ID, "6/7")
    }
    if(query_lst[i]%in% lst7){
      group_ID <- c(group_ID, "7/7")
    }
  }
  print(length(group_ID))
  return(group_ID)
}

#==================================================
#MALAT1

MALAT1_groupID <- get_group_ID(ra_MALAT1[,1],MALAT1_appear1, MALAT1_appear2, MALAT1_appear3, 
                               MALAT1_appear4, MALAT1_appear5, MALAT1_appear6, NULL)

ra_MALAT1 <- cbind(ra_MALAT1, MALAT1_groupID)
ra_MALAT1[,1] <- as.character(ra_MALAT1[,1])
ra_MALAT1[,3] <- as.character(ra_MALAT1[,3])

#NEAT1
NEAT1_groupID <- get_group_ID(ra_NEAT1[,1], NEAT1_appear1, NEAT1_appear2, NEAT1_appear3, 
                              NEAT1_appear4, NEAT1_appear5, NEAT1_appear6, NEAT1_appear7)

ra_NEAT1 <- cbind(ra_NEAT1, NEAT1_groupID)
ra_NEAT1[,1] <- as.character(ra_NEAT1[,1])
ra_NEAT1[,3] <- as.character(ra_NEAT1[,3])


#====================================================
# box-plot

# log RA-score to better visualize 
MALAT1_bp <- ggplot(ra_MALAT1, aes(x = ra_MALAT1$MALAT1_groupID, y = log(ra_MALAT1$Score),fill=MALAT1_groupID)) + labs(x="MALAT1 Gene Groups",y="RA Score") + geom_boxplot(outlier.size=0.3)

# log RA-score to better visualize 
NEAT1_bp <- ggplot(ra_NEAT1, aes(x = ra_NEAT1$NEAT1_groupID, y = log(ra_NEAT1$Score),fill=NEAT1_groupID)) + labs(x="NEAT1 Gene Groups",y="RA Score") + geom_boxplot(outlier.size=0.3)

#===================
#END