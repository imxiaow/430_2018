#========================================================================================
# 1. filter patients with NULL cancertype,
# 2. log1p both datasets.(pcgs and lncRNAs). Subset dataset into cancer-type specific.
# 3. for each cancer type, find significant correlated pcgs, 18 pcg-lists for each lncRNAs.
# 4. RA, combine them together, filter genes with RA-score <0.05
# 5. compare the result with PCAWGs.
#========================================
# load needed dataset, saved from previous script.

#pcgs_with_canctypes
load("~/Documents/BCB430/data2_Jan/TCGA_data/all_pcgs_with_canctypes_7387p.Rda")
#NEAT1_all_p
load("~/Documents/BCB430/data2_Jan/TCGA_data/NEAT1_all_patients7387.Rda")
#MALAT1_all_p
load("~/Documents/BCB430/data2_Jan/TCGA_data/MALAT1_all_patients7387.Rda")

#cancertypes_plus_numbers
load("~/Documents/BCB430/data2_Jan/TCGA_data/cancertypes_and_itsnumberofpatietns.Rda")

#=================================================================================
# 1. filter out patients with cancer type 1: ""
#    in 1) pcgs_with_canctypes
#       2) NEAT1_all_p
#       3) MALAT1_all_p
# delete patients with ""

p <- rownames(pcgs_with_canctypes)
c <- pcgs_with_canctypes[,19439]
pID_with_canctype <- cbind(p,c)
pID_with_canctype <- as.data.frame(pID_with_canctype, stringsAsFactors = FALSE)
save(pID_with_canctype, file="patientID_with_canctypes_unfilter.Rda")

need_filter_ctype <- cancertypes_plus_numbers$typesofcancer[1]
#length(pID_with_canctype$p[which(need_filter_ctype == pID_with_canctype$c)])
#[1] 63  -> patients

#(0)
# get 63 unique patient IDs
need_filter_pID <- pID_with_canctype$p[which(need_filter_ctype == pID_with_canctype$c)]


# (1)
# subset NEAT1_all_p with pID not in need_filter_pID
c1_filtered_NEAT1_all_p <- subset(NEAT1_all_p, (rownames(NEAT1_all_p)%in%need_filter_pID)==FALSE)
# 7324 patient left, 7387-63 = 7324.
save(c1_filtered_NEAT1_all_p, file="c1_filtered_NEAT1_patients7324.Rda")
# check again to see if filtered data contain any PID that need to be filtered.
# which(need_filter_pID%in%rownames(c1_filtered_NEAT1_all_p))
#integer(0)


# (2)
# subset MALAT1_all_p with pID not in need_filter_pID
c1_filtered_MALAT1_all_p <- subset(MALAT1_all_p, (rownames(MALAT1_all_p)%in%need_filter_pID)==FALSE)
# 7324 left.
save(c1_filtered_MALAT1_all_p, file="c1_filtered_MALAT1_patients7324.Rda")
# check again
# which(need_filter_pID%in%rownames(c1_filtered_MALAT1_all_p))
#integer(0)


# (3)
c1_filtered_pcgs_with_canctypes <- subset(pcgs_with_canctypes, (rownames(pcgs_with_canctypes)%in%need_filter_pID)==FALSE)
#7324 left
save(c1_filtered_pcgs_with_canctypes, file="c1_filtered_pcgs_with_canctypes_patients7324.Rda")

#=================================================================================


#============================
# Output Rda files from this script, that can load later:
#   
#  1) filtered cancertype1 for NEAT1, 7324 patients left. 
#     "~/Documents/BCB430/data2_Jan/TCGA_data/c1_filtered_NEAT1_patients7324.Rda"
#     #c1_filtered_NEAT1_all_p       
#
#  2) filtered cancertype1 for MALAT1, 7324 patients left.
#     "~/Documents/BCB430/data2_Jan/TCGA_data/c1_filtered_MALAT1_patients7324.Rda"
#     #c1_filtered_MALAT1_all_p
#     
#  3) filtered cancertype1 for all_pcgs,  7324 patients left.
#     "~/Documents/BCB430/data2_Jan/TCGA_data/c1_filtered_pcgs_with_canctypes_patients7324.Rda"
#     #c1_filtered_pcgs_with_canctypes
#
#
#============
#END