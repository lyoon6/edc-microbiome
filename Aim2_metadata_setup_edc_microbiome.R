#################################################
#    lara yoon                                  #
#     lyoon6@ucla.edu                           #
#     EDC- MICROBIOME METADATA SETUP            #
#     last updated 07/27/2021                   #
#################################################

library(pacman)
pacman::p_load("devtools", "tidyverse", "reshape2", "ggplot2", "gtools", 
               "naniar", "mice", "lattice", "readxl",
               "scales", "grid", "ggsci", "ggpubr", "viridis", 
               "picante", "knitr", "gridExtra")  

pacman::p_loaded()

#############################################
##                                         ##
##       Loading the Data Sets             ##
##                                         ##
#############################################

# Cleaned & imputed data from "CHEAR_Readin" file
dat <- read_csv("Data/FinalData/AIM2_EDC_MICRO.csv")
summary(dat)

# Stool sample manifest 
manifest <- read.delim("Microbiome/Input/metadata.txt")
manifest.sub <- manifest %>%  
  filter(PID!="Control") %>% 
  filter(QC.flag!="Duplicate= MICRO") %>%  
  select(c("sampleid", "Sample.Name", "PID")) %>%  
  dplyr::rename("chear_pid"="PID")


# WQS index from CHEAR microbiome analysis
wqs <- read.csv("C:/Users/laray/Box/Karin's Lab/GOCS/CHEAR/Final Reports and Code/Rerun/Output/WQS/wqsindexdata_pos_pm1.csv")
wqs2 <- wqs %>%  
  rename(chear_pid=pid)  
wqs2$chear_pid <- as.character(wqs2$chear_pid)

#############################################
##                                         ##
##       Tidy                              ##
##                                         ##
#############################################

names(dat)

dat1 <- dat %>%  
  select(-c("X1")) %>%  
  mutate(menarche=ifelse(age_men<age, 1, 0)) %>%  
  mutate(log_dehp = log10(mol_dehp), 
         log_hiphth = log10(mol_hi_phth), 
         log_lophth = log10(mol_lo_phth), 
         log_phenf = log10(mol_phenf), 
         log_parbf = log10(mol_parbf)) %>%  
  mutate(chear_pid=as.character(chear_pid))


myfactvars <- c("lab", "moments", "bfeedcat", "dp_mapuche", "fatcat", "fatcat_2pm","birth_mode", "menarche", "medu", "urine_b1", "urine_b4", "urine_1pm", "bd_2pm", "stool", "abx")
dat1[,myfactvars] <- lapply(dat1[,myfactvars], factor)


#############################################
##                                         ##
##       Merge                             ##
##                                         ##
#############################################

metaall <- left_join(dat1, manifest.sub, by="chear_pid")
metaall <- metaall %>%
  select( sampleid, Sample.Name, chear_pid, everything())

metawqs <- left_join(metaall, wqs2, by="chear_pid")
metaall <- metawqs

n_distinct(metaall$chear_pid) #261


#############################################
##                                         ##
##       Cut by Tanner                     ##
##                                         ##
#############################################

metab1 <- metaall %>%  
  filter(moments=="b1" & urine_b1==1) %>%  
  filter(!is.na(lab))

metab4 <- metaall %>%  
  filter(moments=="b4" & urine_b4==1)%>%  
  filter(!is.na(lab))

metapm1 <- metaall %>%  
  filter(moments=="pm1" & urine_1pm==1)%>%  
  filter(!is.na(lab))

meta <- metaall %>%  
  select(sampleid, Sample.Name, chear_pid, urine_b1, urine_b4, urine_1pm, stool) %>%  
  distinct()


#############################################
##                                         ##
##       Create tertiles                   ##
##                                         ##
#############################################

## Should probably loop this... 

# B1 
metab1$bp3T <- ntile(metab1$log_bp3,3)
metab1$bpaT <- ntile(metab1$log_bpa,3)
metab1$bpsT <- ntile(metab1$log_bps,3)
metab1$tcsT <- ntile(metab1$log_tcs,3)
metab1$phenT <- ntile(metab1$log_phenf,3)
metab1$mepbT <- ntile(metab1$log_mepb,3)
metab1$etpbT <- ntile(metab1$log_etpb,3)
metab1$prpbT <- ntile(metab1$log_prpb,3)
metab1$parbT <- ntile(metab1$log_parbf,3)
metab1$mbpT <- ntile(metab1$log_mbp,3)
metab1$mbzpT <- ntile(metab1$log_mbzp,3)
metab1$mcppT <- ntile(metab1$log_mcpp,3)
metab1$mecppT <- ntile(metab1$log_mecpp,3)
metab1$mehhpT <- ntile(metab1$log_mehhp,3)
metab1$mehpT <- ntile(metab1$log_mehp,3)
metab1$meohpT <- ntile(metab1$log_meohp,3)
metab1$mepT <- ntile(metab1$log_mep,3)
metab1$mibpT <- ntile(metab1$log_mibp,3)
metab1$dehpT <- ntile(metab1$log_dehp,3)
metab1$hiphthT <- ntile(metab1$log_hiphth,3)
metab1$lophthT <- ntile(metab1$log_lophth,3)
# B4
metab4$bp3T <- ntile(metab4$log_bp3,3)
metab4$bpaT <- ntile(metab4$log_bpa,3)
metab4$bpsT <- ntile(metab4$log_bps,3)
metab4$tcsT <- ntile(metab4$log_tcs,3)
metab4$phenT <- ntile(metab4$log_phenf,3)
metab4$mepbT <- ntile(metab4$log_mepb,3)
metab4$etpbT <- ntile(metab4$log_etpb,3)
metab4$prpbT <- ntile(metab4$log_prpb,3)
metab4$parbT <- ntile(metab4$log_parbf,3)
metab4$mbpT <- ntile(metab4$log_mbp,3)
metab4$mbzpT <- ntile(metab4$log_mbzp,3)
metab4$mcppT <- ntile(metab4$log_mcpp,3)
metab4$mecppT <- ntile(metab4$log_mecpp,3)
metab4$mehhpT <- ntile(metab4$log_mehhp,3)
metab4$mehpT <- ntile(metab4$log_mehp,3)
metab4$meohpT <- ntile(metab4$log_meohp,3)
metab4$mepT <- ntile(metab4$log_mep,3)
metab4$mibpT <- ntile(metab4$log_mibp,3)
metab4$dehpT <- ntile(metab4$log_dehp,3)
metab4$hiphthT <- ntile(metab4$log_hiphth,3)
metab4$lophthT <- ntile(metab4$log_lophth,3)
# PM1
metapm1$bp3T <- ntile(metapm1$log_bp3,3)
metapm1$bpaT <- ntile(metapm1$log_bpa,3)
metapm1$bpsT <- ntile(metapm1$log_bps,3)
metapm1$tcsT <- ntile(metapm1$log_tcs,3)
metapm1$phenT <- ntile(metapm1$log_phenf,3)
metapm1$mepbT <- ntile(metapm1$log_mepb,3)
metapm1$etpbT <- ntile(metapm1$log_etpb,3)
metapm1$prpbT <- ntile(metapm1$log_prpb,3)
metapm1$parbT <- ntile(metapm1$log_parbf,3)
metapm1$mbpT <- ntile(metapm1$log_mbp,3)
metapm1$mbzpT <- ntile(metapm1$log_mbzp,3)
metapm1$mcppT <- ntile(metapm1$log_mcpp,3)
metapm1$mecppT <- ntile(metapm1$log_mecpp,3)
metapm1$mehhpT <- ntile(metapm1$log_mehhp,3)
metapm1$mehpT <- ntile(metapm1$log_mehp,3)
metapm1$meohpT <- ntile(metapm1$log_meohp,3)
metapm1$mepT <- ntile(metapm1$log_mep,3)
metapm1$mibpT <- ntile(metapm1$log_mibp,3)
metapm1$dehpT <- ntile(metapm1$log_dehp,3)
metapm1$hiphthT <- ntile(metapm1$log_hiphth,3)
metapm1$lophthT <- ntile(metapm1$log_lophth,3)

# Make factors 
# B1 
metab1$bp3T <- as_factor(metab1$bp3T)
metab1$bpaT <- as_factor(metab1$bpaT)
metab1$bpsT <- as_factor(metab1$bpsT)
metab1$tcsT <- as_factor(metab1$tcsT)
metab1$phenT <- as_factor(metab1$phenT)
metab1$mepbT <- as_factor(metab1$mepbT)
metab1$etpbT <- as_factor(metab1$etpbT)
metab1$prpbT <- as_factor(metab1$prpbT)
metab1$parbT <- as_factor(metab1$parbT)
metab1$mbpT <- as_factor(metab1$mbpT)
metab1$mbzpT <- as_factor(metab1$mbzpT)
metab1$mcppT <- as_factor(metab1$mcppT)
metab1$mecppT <- as_factor(metab1$mecppT)
metab1$mehhpT <- as_factor(metab1$mehhpT)
metab1$mehpT <- as_factor(metab1$mehpT)
metab1$meohpT <- as_factor(metab1$meohpT)
metab1$mepT <- as_factor(metab1$mepT)
metab1$mibpT <- as_factor(metab1$mibpT)
metab1$dehpT <- as_factor(metab1$dehpT)
metab1$hiphthT <- as_factor(metab1$hiphthT)
metab1$lophthT <- as_factor(metab1$lophthT)
# B4
metab4$bp3T <- as_factor(metab4$bp3T)
metab4$bpaT <- as_factor(metab4$bpaT)
metab4$bpsT <- as_factor(metab4$bpsT)
metab4$tcsT <- as_factor(metab4$tcsT)
metab4$phenT <- as_factor(metab4$phenT)
metab4$mepbT <- as_factor(metab4$mepbT)
metab4$etpbT <- as_factor(metab4$etpbT)
metab4$prpbT <- as_factor(metab4$prpbT)
metab4$parbT <- as_factor(metab4$parbT)
metab4$mbpT <- as_factor(metab4$mbpT)
metab4$mbzpT <- as_factor(metab4$mbzpT)
metab4$mcppT <- as_factor(metab4$mcppT)
metab4$mecppT <- as_factor(metab4$mecppT)
metab4$mehhpT <- as_factor(metab4$mehhpT)
metab4$mehpT <- as_factor(metab4$mehpT)
metab4$meohpT <- as_factor(metab4$meohpT)
metab4$mepT <- as_factor(metab4$mepT)
metab4$mibpT <- as_factor(metab4$mibpT)
metab4$dehpT <- as_factor(metab4$dehpT)
metab4$hiphthT <- as_factor(metab4$hiphthT)
metab4$lophthT <- as_factor(metab4$lophthT)
# PM1
metapm1$bp3T <- as_factor(metapm1$bp3T)
metapm1$bpaT <- as_factor(metapm1$bpaT)
metapm1$bpsT <- as_factor(metapm1$bpsT)
metapm1$tcsT <- as_factor(metapm1$tcsT)
metapm1$phenT <- as_factor(metapm1$phenT)
metapm1$mepbT <- as_factor(metapm1$mepbT)
metapm1$etpbT <- as_factor(metapm1$etpbT)
metapm1$prpbT <- as_factor(metapm1$prpbT)
metapm1$parbT <- as_factor(metapm1$parbT)
metapm1$mbpT <- as_factor(metapm1$mbpT)
metapm1$mbzpT <- as_factor(metapm1$mbzpT)
metapm1$mcppT <- as_factor(metapm1$mcppT)
metapm1$mecppT <- as_factor(metapm1$mecppT)
metapm1$mehhpT <- as_factor(metapm1$mehhpT)
metapm1$mehpT <- as_factor(metapm1$mehpT)
metapm1$meohpT <- as_factor(metapm1$meohpT)
metapm1$mepT <- as_factor(metapm1$mepT)
metapm1$mibpT <- as_factor(metapm1$mibpT)
metapm1$dehpT <- as_factor(metapm1$dehpT)
metapm1$hiphthT <- as_factor(metapm1$hiphthT)
metapm1$lophthT <- as_factor(metapm1$lophthT)


#############################################
##                                         ##
##       Save                              ##
##                                         ##
#############################################

# full long
write.csv(metaall, "Data/FinalData/Aim2MetadataFinal.csv", row.names = FALSE)
write.table(metaall, "Data/FinalData/Aim2MetadataFinal.txt", row.names = FALSE, sep="\t", quote = FALSE, col.names=TRUE)
write.table(meta,"Data/FinalData/Aim2MetadataFinalDistinct.txt", row.names = FALSE, sep="\t", quote = FALSE, col.names=TRUE)

# strat by tanner 
write.table(metab1, "Data/FinalData/Aim2MetadataFinalB1.txt", row.names = FALSE, sep="\t", quote = FALSE, col.names=TRUE)
write.table(metab4, "Data/FinalData/Aim2MetadataFinalB4.txt", row.names = FALSE, sep="\t", quote = FALSE, col.names=TRUE)
write.table(metapm1, "Data/FinalData/Aim2MetadataFinalPM1.txt", row.names = FALSE, sep="\t", quote = FALSE, col.names=TRUE)

# save as csv
write.csv(metab1, "Data/FinalData/Aim2MetadataFinalB1.csv", row.names = FALSE)
write.csv(metab4, "Data/FinalData/Aim2MetadataFinalB4.csv", row.names = FALSE)
write.csv(metapm1, "Data/FinalData/Aim2MetadataFinalPM1.csv", row.names = FALSE)




