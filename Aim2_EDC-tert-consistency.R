#################################################
#    lara yoon                                  #
#     lyoon6@ucla.edu                           #
#     EDC- TERTILE CONSISTENCY (AIM2)           #
#     last updated 08/08/2021                   #
#################################################

library(pacman)
pacman::p_load("devtools", "tidyverse", "reshape2", "ggplot2", "gtools", 
               "naniar", "mice", "lattice", "readxl",
               "scales", "grid", "ggsci", "ggpubr", "viridis", 
               "picante", "knitr", "gridExtra", "tableone", "gmodels")  

pacman::p_loaded()


#############################################
##                                         ##
##       Load data                         ##
##                                         ##
#############################################

# From metadata setup file 
###############################################

metab1 <- read_csv("Data/FinalData/Aim2MetadataFinalB1.csv")
metab4 <- read_csv("Data/FinalData/Aim2MetadataFinalB4.csv")
metapm1 <- read_csv("Data/FinalData/Aim2MetadataFinalPM1.csv")


#############################################
##                                         ##
##       Clean up                          ##
##                                         ##
#############################################

# Make factor variables 
myfacts <- c("bp3T", "bpaT", "bpsT", "tcsT", "phenT", 
             "mepbT", "etpbT", "prpbT", "parbT", 
             "mbpT", "mbzpT", "mcppT", "mecppT", "mehhpT", 
             "mehpT", "meohpT", "mepT", "mibpT", 
             "dehpT", "hiphthT", "lophthT")
metab1[,myfacts] <- lapply(metab1[,myfacts] , factor)
metab4[,myfacts] <- lapply(metab4[,myfacts] , factor)
metapm1[,myfacts] <- lapply(metapm1[,myfacts] , factor)


# Keep only relevant information 
myvars <- c("chear_pid", "bp3T", "bpaT", "bpsT", "tcsT", "phenT", 
            "mepbT", "etpbT", "prpbT", "parbT", 
            "mbpT", "mbzpT", "mcppT", "mecppT", "mehhpT", 
            "mehpT", "meohpT", "mepT", "mibpT", 
            "dehpT", "hiphthT", "lophthT", 
            "log_bp3", "log_bpa", "log_bps", "log_tcs", "log_phenf",
            "log_mepb", "log_etpb", "log_prpb", "log_parbf",
            "log_mbp", "log_mbzp", "log_mcpp", "log_mecpp","log_mehhp",  
            "log_mehp", "log_meohp", "log_mep", "log_mibp", 
            "log_dehp", "log_hiphth", "log_lophth") 
metab1 <- metab1[,myvars]
metab4 <- metab4[,myvars]
metapm1 <- metapm1[,myvars]

# Add prefix 
colnames(metab1)[2:43] <- paste(colnames(metab1[,c(2:43)]), "_b1", sep = "")
colnames(metab4)[2:43] <- paste(colnames(metab4[,c(2:43)]), "_b4", sep = "")
colnames(metapm1)[2:43] <- paste(colnames(metapm1[,c(2:43)]), "_pm1", sep = "")

# Merge together 
edcall <- metab1 %>%  
  full_join(metab4, by="chear_pid") %>%  
  full_join(metapm1, by="chear_pid")


#############################################
##                                         ##
##       Look at frequency of tertiles     ##
##                                         ##
#############################################

myedcs <- ("bp3", "bpa", "bps", "tcs", "phen", 
"mepb", "etpb", "prpb", "parb", 
"mbp", "mbzp", "mcpp", "mecpp", "mehhp", 
"mehp", "meohp", "mep", "mibp", "dehp", "hiphth", "lophth")

freqfun <- function(edc){}


  
  # Get frequencies
  lophth_b1_b4 <- CrossTable(edcall$lophthT_b1, edcall$lophthT_b4)
  lophth_b1_pm1 <- CrossTable(edcall$lophthT_b1, edcall$lophthT_pm1)
  lophth_b4_pm1 <- CrossTable(edcall$lophthT_b4, edcall$lophthT_pm1)
  
  # Convert row frequencies to table 
  lophth_b1_b4.df <- as.data.frame(lophth_b1_b4[["prop.row"]])
  lophth_b1_b4.df <- lophth_b1_b4.df %>%  rename(B1T=x, B4T=y) %>%  arrange(B1T)
  lophth_b1_pm1.df <- as.data.frame(lophth_b1_pm1[["prop.row"]])
  lophth_b1_pm1.df <- lophth_b1_pm1.df %>%  rename(B1T=x, PM1T=y) %>%  arrange(B1T)
  lophth_b4_pm1.df <- as.data.frame(lophth_b4_pm1[["prop.row"]])
  lophth_b4_pm1.df <- lophth_b4_pm1.df %>%  rename(B4T=x, PM1T=y) %>%  arrange(B4T)
  
  # Summary table plot
  lophth_b1_b4.p <- ggtexttable(lophth_b1_b4.df, rows = NULL)
  lophth_b1_pm1.p <- ggtexttable(lophth_b1_pm1.df, rows = NULL)
  lophth_b4_pm1.p <- ggtexttable(lophth_b4_pm1.df, rows = NULL)
  
  # Arrange together 
  lophth.p <- ggarrange(lophth_b1_b4.p, lophth_b1_pm1.p, lophth_b4_pm1.p, 
            ncol = 3, nrow = 1)
  annotate_figure(lophth.p, top = text_grob("lophth Tertiles", face = "bold"))

