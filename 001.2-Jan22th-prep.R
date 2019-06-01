#====================================================================================
# 2. log1p both datasets.(pcgs and lncRNAs). Subset dataset into cancer-type specific.
# 3. for each cancer type, find significant correlated pcgs, 18 pcg-lists for each lncRNAs.
# 4. RA, combine them together, filter genes with RA-score <0.05
# 5. compare the result with PCAWGs.

#========================================
# load needed dataset, saved from previous script.

# 1) filtered cancertype1 for NEAT1, 7324 patients left. 
load("~/Documents/BCB430/data2_Jan/TCGA_data/c1_filtered_NEAT1_patients7324.Rda")
#c1_filtered_NEAT1_all_p       

# 2) filtered cancertype1 for MALAT1, 7324 patients left.
load("~/Documents/BCB430/data2_Jan/TCGA_data/c1_filtered_MALAT1_patients7324.Rda")
#c1_filtered_MALAT1_all_p
     
# 3) filtered cancertype1 for all_pcgs,  7324 patients left.
load("~/Documents/BCB430/data2_Jan/TCGA_data/c1_filtered_pcgs_with_canctypes_patients7324.Rda")
#c1_filtered_pcgs_with_canctypes

# 4) cancertypes_plus_numbers
load("~/Documents/BCB430/data2_Jan/TCGA_data/cancertypes_and_itsnumberofpatietns.Rda")
#cancertypes_plus_numbers

# 5) all patientIDs with canctypes unfilter
load("~/Documents/BCB430/data2_Jan/TCGA_data/patientID_with_canctypes_unfilter.Rda")
#pID_with_canctype

#======================================================================
# 2a. log1p both datasets.(pcgs and lncRNAs).
#           in 1) c1_filtered_NEAT1_all_p
#              2) c1_filtered_MALAT1_all_p
#              3) c1_filtered_pcgs_with_canctypes
# log1p to the whole matrix

log1p_c1_filtered_NEAT1_all_p <- as.data.frame(apply(c1_filtered_NEAT1_all_p,2, log1p))

log1p_c1_filtered_MALAT1_all_p <- as.data.frame(apply(c1_filtered_MALAT1_all_p,2, log1p))

log1p_c1_filtered_pcgs_with_canctypes <- as.data.frame(apply(c1_filtered_pcgs_with_canctypes[,1:19438],2, log1p))

save(log1p_c1_filtered_NEAT1_all_p, file="log1p_c1_filtered_NEAT1_all_patients7324.Rda")
save(log1p_c1_filtered_MALAT1_all_p, file="log1p_c1_filtered_MALAT1_all_patients7324.Rda")
save(log1p_c1_filtered_pcgs_with_canctypes, file = "log1p_c1_filtered_pcgs_with_canctypes_patients7324.Rda")


#======================================================================
# 2b.Subset dataset into cancer-type specific.

# get each can type string
c1_type <- cancertypes_plus_numbers$typesofcancer[2]
c2_type <- cancertypes_plus_numbers$typesofcancer[3]
c3_type <- cancertypes_plus_numbers$typesofcancer[4]
c4_type <- cancertypes_plus_numbers$typesofcancer[5]
c5_type <- cancertypes_plus_numbers$typesofcancer[6]
c6_type <- cancertypes_plus_numbers$typesofcancer[7]
c7_type <- cancertypes_plus_numbers$typesofcancer[8]
c8_type <- cancertypes_plus_numbers$typesofcancer[9]
c9_type <- cancertypes_plus_numbers$typesofcancer[10]
c10_type <- cancertypes_plus_numbers$typesofcancer[11]

c11_type <- cancertypes_plus_numbers$typesofcancer[12]
c12_type <- cancertypes_plus_numbers$typesofcancer[13]
c13_type <- cancertypes_plus_numbers$typesofcancer[14]
c14_type <- cancertypes_plus_numbers$typesofcancer[15]
c15_type <- cancertypes_plus_numbers$typesofcancer[16]
c16_type <- cancertypes_plus_numbers$typesofcancer[17]
c17_type <- cancertypes_plus_numbers$typesofcancer[18]
c18_type <- cancertypes_plus_numbers$typesofcancer[19]

# get correspond pids
c1_PID <- pID_with_canctype$p[which(pID_with_canctype$c == c1_type)]
c2_PID <- pID_with_canctype$p[which(pID_with_canctype$c == c2_type)]
c3_PID <- pID_with_canctype$p[which(pID_with_canctype$c == c3_type)]
c4_PID <- pID_with_canctype$p[which(pID_with_canctype$c == c4_type)]
c5_PID <- pID_with_canctype$p[which(pID_with_canctype$c == c5_type)]
c6_PID <- pID_with_canctype$p[which(pID_with_canctype$c == c6_type)]
c7_PID <- pID_with_canctype$p[which(pID_with_canctype$c == c7_type)]
c8_PID <- pID_with_canctype$p[which(pID_with_canctype$c == c8_type)]
c9_PID <- pID_with_canctype$p[which(pID_with_canctype$c == c9_type)]
c10_PID <- pID_with_canctype$p[which(pID_with_canctype$c == c10_type)]

c11_PID <- pID_with_canctype$p[which(pID_with_canctype$c == c11_type)]
c12_PID <- pID_with_canctype$p[which(pID_with_canctype$c == c12_type)]
c13_PID <- pID_with_canctype$p[which(pID_with_canctype$c == c13_type)]
c14_PID <- pID_with_canctype$p[which(pID_with_canctype$c == c14_type)]
c15_PID <- pID_with_canctype$p[which(pID_with_canctype$c == c15_type)]
c16_PID <- pID_with_canctype$p[which(pID_with_canctype$c == c16_type)]
c17_PID <- pID_with_canctype$p[which(pID_with_canctype$c == c17_type)]
c18_PID <- pID_with_canctype$p[which(pID_with_canctype$c == c18_type)]

#subset NEAT1 into 18 cantyep list.
c1_NEAT1 <- subset(log1p_c1_filtered_NEAT1_all_p, (rownames(log1p_c1_filtered_NEAT1_all_p)%in% c1_PID)==TRUE)
c2_NEAT1 <- subset(log1p_c1_filtered_NEAT1_all_p, (rownames(log1p_c1_filtered_NEAT1_all_p)%in% c2_PID)==TRUE)
c3_NEAT1 <- subset(log1p_c1_filtered_NEAT1_all_p, (rownames(log1p_c1_filtered_NEAT1_all_p)%in% c3_PID)==TRUE)
c4_NEAT1 <- subset(log1p_c1_filtered_NEAT1_all_p, (rownames(log1p_c1_filtered_NEAT1_all_p)%in% c4_PID)==TRUE)
c5_NEAT1 <- subset(log1p_c1_filtered_NEAT1_all_p, (rownames(log1p_c1_filtered_NEAT1_all_p)%in% c5_PID)==TRUE)
c6_NEAT1 <- subset(log1p_c1_filtered_NEAT1_all_p, (rownames(log1p_c1_filtered_NEAT1_all_p)%in% c6_PID)==TRUE)
c7_NEAT1 <- subset(log1p_c1_filtered_NEAT1_all_p, (rownames(log1p_c1_filtered_NEAT1_all_p)%in% c7_PID)==TRUE)
c8_NEAT1 <- subset(log1p_c1_filtered_NEAT1_all_p, (rownames(log1p_c1_filtered_NEAT1_all_p)%in% c8_PID)==TRUE)
c9_NEAT1 <- subset(log1p_c1_filtered_NEAT1_all_p, (rownames(log1p_c1_filtered_NEAT1_all_p)%in% c9_PID)==TRUE)
c10_NEAT1 <- subset(log1p_c1_filtered_NEAT1_all_p, (rownames(log1p_c1_filtered_NEAT1_all_p)%in% c10_PID)==TRUE)

c11_NEAT1 <- subset(log1p_c1_filtered_NEAT1_all_p, (rownames(log1p_c1_filtered_NEAT1_all_p)%in% c11_PID)==TRUE)
c12_NEAT1 <- subset(log1p_c1_filtered_NEAT1_all_p, (rownames(log1p_c1_filtered_NEAT1_all_p)%in% c12_PID)==TRUE)
c13_NEAT1 <- subset(log1p_c1_filtered_NEAT1_all_p, (rownames(log1p_c1_filtered_NEAT1_all_p)%in% c13_PID)==TRUE)
c14_NEAT1 <- subset(log1p_c1_filtered_NEAT1_all_p, (rownames(log1p_c1_filtered_NEAT1_all_p)%in% c14_PID)==TRUE)
c15_NEAT1 <- subset(log1p_c1_filtered_NEAT1_all_p, (rownames(log1p_c1_filtered_NEAT1_all_p)%in% c15_PID)==TRUE)
c16_NEAT1 <- subset(log1p_c1_filtered_NEAT1_all_p, (rownames(log1p_c1_filtered_NEAT1_all_p)%in% c16_PID)==TRUE)
c17_NEAT1 <- subset(log1p_c1_filtered_NEAT1_all_p, (rownames(log1p_c1_filtered_NEAT1_all_p)%in% c17_PID)==TRUE)
c18_NEAT1 <- subset(log1p_c1_filtered_NEAT1_all_p, (rownames(log1p_c1_filtered_NEAT1_all_p)%in% c18_PID)==TRUE)

save(c1_NEAT1, file= "c1_NEAT1_382p.Rda")
save(c2_NEAT1, file= "c2_NEAT1_488p.Rda")
save(c3_NEAT1, file= "c3_NEAT1_987p.Rda")
save(c4_NEAT1, file= "c4_NEAT1_284p.Rda")
save(c5_NEAT1, file= "c5_NEAT1_396p.Rda")
save(c6_NEAT1, file= "c6_NEAT1_456p.Rda")
save(c7_NEAT1, file= "c7_NEAT1_484p.Rda")
save(c8_NEAT1, file= "c8_NEAT1_255p.Rda")
save(c9_NEAT1, file= "c9_NEAT1_319p.Rda")
save(c10_NEAT1, file= "c10_NEAT1_432p.Rda")
save(c11_NEAT1, file= "c11_NEAT1_453p.Rda")
save(c12_NEAT1, file= "c12_NEAT1_305p.Rda")
save(c13_NEAT1, file= "c13_NEAT1_472p.Rda")
save(c14_NEAT1, file= "c14_NEAT1_224p.Rda")
save(c15_NEAT1, file= "c15_NEAT1_99p.Rda")
save(c16_NEAT1, file= "c16_NEAT1_345p.Rda")
save(c17_NEAT1, file= "c17_NEAT1_453p.Rda")
save(c18_NEAT1, file= "c18_NEAT1_490p.Rda")

# subset MALAT1 into 18 cantype list.
c1_MALAT1 <- subset(log1p_c1_filtered_MALAT1_all_p, (rownames(log1p_c1_filtered_MALAT1_all_p)%in% c1_PID)==TRUE)
c2_MALAT1 <- subset(log1p_c1_filtered_MALAT1_all_p, (rownames(log1p_c1_filtered_MALAT1_all_p)%in% c2_PID)==TRUE)
c3_MALAT1 <- subset(log1p_c1_filtered_MALAT1_all_p, (rownames(log1p_c1_filtered_MALAT1_all_p)%in% c3_PID)==TRUE)
c4_MALAT1 <- subset(log1p_c1_filtered_MALAT1_all_p, (rownames(log1p_c1_filtered_MALAT1_all_p)%in% c4_PID)==TRUE)
c5_MALAT1 <- subset(log1p_c1_filtered_MALAT1_all_p, (rownames(log1p_c1_filtered_MALAT1_all_p)%in% c5_PID)==TRUE)
c6_MALAT1 <- subset(log1p_c1_filtered_MALAT1_all_p, (rownames(log1p_c1_filtered_MALAT1_all_p)%in% c6_PID)==TRUE)
c7_MALAT1 <- subset(log1p_c1_filtered_MALAT1_all_p, (rownames(log1p_c1_filtered_MALAT1_all_p)%in% c7_PID)==TRUE)
c8_MALAT1 <- subset(log1p_c1_filtered_MALAT1_all_p, (rownames(log1p_c1_filtered_MALAT1_all_p)%in% c8_PID)==TRUE)
c9_MALAT1 <- subset(log1p_c1_filtered_MALAT1_all_p, (rownames(log1p_c1_filtered_MALAT1_all_p)%in% c9_PID)==TRUE)
c10_MALAT1 <- subset(log1p_c1_filtered_MALAT1_all_p, (rownames(log1p_c1_filtered_MALAT1_all_p)%in% c10_PID)==TRUE)

c11_MALAT1 <- subset(log1p_c1_filtered_MALAT1_all_p, (rownames(log1p_c1_filtered_MALAT1_all_p)%in% c11_PID)==TRUE)
c12_MALAT1 <- subset(log1p_c1_filtered_MALAT1_all_p, (rownames(log1p_c1_filtered_MALAT1_all_p)%in% c12_PID)==TRUE)
c13_MALAT1 <- subset(log1p_c1_filtered_MALAT1_all_p, (rownames(log1p_c1_filtered_MALAT1_all_p)%in% c13_PID)==TRUE)
c14_MALAT1 <- subset(log1p_c1_filtered_MALAT1_all_p, (rownames(log1p_c1_filtered_MALAT1_all_p)%in% c14_PID)==TRUE)
c15_MALAT1 <- subset(log1p_c1_filtered_MALAT1_all_p, (rownames(log1p_c1_filtered_MALAT1_all_p)%in% c15_PID)==TRUE)
c16_MALAT1 <- subset(log1p_c1_filtered_MALAT1_all_p, (rownames(log1p_c1_filtered_MALAT1_all_p)%in% c16_PID)==TRUE)
c17_MALAT1 <- subset(log1p_c1_filtered_MALAT1_all_p, (rownames(log1p_c1_filtered_MALAT1_all_p)%in% c17_PID)==TRUE)
c18_MALAT1 <- subset(log1p_c1_filtered_MALAT1_all_p, (rownames(log1p_c1_filtered_MALAT1_all_p)%in% c18_PID)==TRUE)

save(c1_MALAT1,file="c1_MALAT1_382p.Rda")
save(c2_MALAT1,file="c2_MALAT1_488p.Rda")
save(c3_MALAT1,file="c3_MALAT1_987p.Rda")
save(c4_MALAT1,file="c4_MALAT1_284p.Rda")
save(c5_MALAT1,file="c5_MALAT1_396p.Rda")
save(c6_MALAT1,file="c6_MALAT1_456p.Rda")
save(c7_MALAT1,file="c7_MALAT1_484p.Rda")
save(c8_MALAT1,file="c8_MALAT1_255p.Rda")
save(c9_MALAT1,file="c9_MALAT1_319p.Rda")
save(c10_MALAT1,file="c10_MALAT1_432p.Rda")

save(c11_MALAT1,file="c11_MALAT1_453p.Rda")
save(c12_MALAT1,file="c12_MALAT1_305p.Rda")
save(c13_MALAT1,file="c13_MALAT1_472p.Rda")
save(c14_MALAT1,file="c14_MALAT1_224p.Rda")
save(c15_MALAT1,file="c15_MALAT1_99p.Rda")
save(c16_MALAT1,file="c16_MALAT1_345p.Rda")
save(c17_MALAT1,file="c17_MALAT1_453p.Rda")
save(c18_MALAT1,file="c18_MALAT1_490p.Rda")

#subset pcg into 18 cantype list.
c1_pcgs <- subset(log1p_c1_filtered_pcgs_with_canctypes, (rownames(log1p_c1_filtered_pcgs_with_canctypes)%in% c1_PID)==TRUE)
c2_pcgs <- subset(log1p_c1_filtered_pcgs_with_canctypes, (rownames(log1p_c1_filtered_pcgs_with_canctypes)%in% c2_PID)==TRUE)
c3_pcgs <- subset(log1p_c1_filtered_pcgs_with_canctypes, (rownames(log1p_c1_filtered_pcgs_with_canctypes)%in% c3_PID)==TRUE)
c4_pcgs <- subset(log1p_c1_filtered_pcgs_with_canctypes, (rownames(log1p_c1_filtered_pcgs_with_canctypes)%in% c4_PID)==TRUE)
c5_pcgs <- subset(log1p_c1_filtered_pcgs_with_canctypes, (rownames(log1p_c1_filtered_pcgs_with_canctypes)%in% c5_PID)==TRUE)
c6_pcgs <- subset(log1p_c1_filtered_pcgs_with_canctypes, (rownames(log1p_c1_filtered_pcgs_with_canctypes)%in% c6_PID)==TRUE)
c7_pcgs <- subset(log1p_c1_filtered_pcgs_with_canctypes, (rownames(log1p_c1_filtered_pcgs_with_canctypes)%in% c7_PID)==TRUE)
c8_pcgs <- subset(log1p_c1_filtered_pcgs_with_canctypes, (rownames(log1p_c1_filtered_pcgs_with_canctypes)%in% c8_PID)==TRUE)
c9_pcgs <- subset(log1p_c1_filtered_pcgs_with_canctypes, (rownames(log1p_c1_filtered_pcgs_with_canctypes)%in% c9_PID)==TRUE)
c10_pcgs <- subset(log1p_c1_filtered_pcgs_with_canctypes, (rownames(log1p_c1_filtered_pcgs_with_canctypes)%in% c10_PID)==TRUE)

c11_pcgs <- subset(log1p_c1_filtered_pcgs_with_canctypes, (rownames(log1p_c1_filtered_pcgs_with_canctypes)%in% c11_PID)==TRUE)
c12_pcgs <- subset(log1p_c1_filtered_pcgs_with_canctypes, (rownames(log1p_c1_filtered_pcgs_with_canctypes)%in% c12_PID)==TRUE)
c13_pcgs <- subset(log1p_c1_filtered_pcgs_with_canctypes, (rownames(log1p_c1_filtered_pcgs_with_canctypes)%in% c13_PID)==TRUE)
c14_pcgs <- subset(log1p_c1_filtered_pcgs_with_canctypes, (rownames(log1p_c1_filtered_pcgs_with_canctypes)%in% c14_PID)==TRUE)
c15_pcgs <- subset(log1p_c1_filtered_pcgs_with_canctypes, (rownames(log1p_c1_filtered_pcgs_with_canctypes)%in% c15_PID)==TRUE)
c16_pcgs <- subset(log1p_c1_filtered_pcgs_with_canctypes, (rownames(log1p_c1_filtered_pcgs_with_canctypes)%in% c16_PID)==TRUE)
c17_pcgs <- subset(log1p_c1_filtered_pcgs_with_canctypes, (rownames(log1p_c1_filtered_pcgs_with_canctypes)%in% c17_PID)==TRUE)
c18_pcgs <- subset(log1p_c1_filtered_pcgs_with_canctypes, (rownames(log1p_c1_filtered_pcgs_with_canctypes)%in% c18_PID)==TRUE)

save(c1_pcgs, file = "c1_pcgs_382p.Rda")
save(c2_pcgs, file="c2_pcgs_488p.Rda")
save(c3_pcgs, file="c3_pcgs_987p.Rda")
save(c4_pcgs, file="c4_pcgs_284p.Rda")
save(c5_pcgs, file="c5_pcgs_396p.Rda")
save(c6_pcgs, file="c6_pcgs_456p.Rda")
save(c7_pcgs, file="c7_pcgs_484p.Rda")
save(c8_pcgs, file="c8_pcgs_255p.Rda")
save(c9_pcgs, file="c9_pcgs_319p.Rda")
save(c10_pcgs, file="c10_pcgs_432p.Rda")
save(c11_pcgs, file="c11_pcgs_453p.Rda")

save(c12_pcgs, file="c12_pcgs_305p.Rda")
save(c13_pcgs, file="c13_pcgs_472p.Rda")
save(c14_pcgs, file="c14_pcgs_224p.Rda")
save(c15_pcgs, file="c15_pcgs_99p.Rda")
save(c16_pcgs, file="c16_pcgs_345p.Rda")
save(c17_pcgs, file="c17_pcgs_453p.Rda")
save(c18_pcgs, file="c18_pcgs_490p.Rda")

#============================
# Output Rda files from this script, that can load later:
# 1) log1p  c1_filtered_NEAT1_all_p
#           "~/Documents/BCB430/data2_Jan/TCGA_data/log1p_c1_filtered_NEAT1_all_patients7324.Rda"
#    log1p_c1_filtered_NEAT1_all_p     
#
#
# 2) log1p  log1p_c1_filtered_MALAT1_all_p
#           "~/Documents/BCB430/data2_Jan/TCGA_data/log1p_c1_filtered_MALAT1_all_patients7324.Rda"
#    log1p_c1_filtered_MALAT1_all_p
#
#
# 3) log1p  c1_filtered_pcgs_with_canctypes
#           ~/Documents/BCB430/data2_Jan/TCGA_data/log1p_c1_filtered_pcgs_with_canctypes_patients7324.Rda"
#     log1p_c1_filtered_pcgs_with_canctypes
#
#
#c1_pcgs ~ c18_pcgs          18
#c1_MALAT1 ~ c18_MALAT1       18
# c1_NEAT1 ~ c18_NEAT1        18 


#============
#END