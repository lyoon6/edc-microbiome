#################################################
#    lara yoon                                  #
#     lyoon6@ucla.edu                           #
#     EDC- MICROBIOME MAASLIN  (AIM2)           #
#     last updated 09/08/2021                   #
#################################################

library(pacman)
pacman::p_load(tidyr, gtools, tidyverse, broom, data.table, ggplot2, sjPlot, gWQS, kableExtra, phyloseq, qiime2R, Maaslin2, here, microbiome)
pacman::p_loaded()


#############################################
##                                         ##
##       Loading the Data Sets             ##
##                                         ##
#############################################

setwd("C:/Users/laray/Box/Karin's Lab/GOCS/CHEAR/Final Reports and Code/Rerun/Statistical_Analysis_Datasets")

# Metadata
#############################################
# Stool sample manifest 
manifest <- read.delim("microbiome/metadata.txt")
manifest.sub <- manifest %>%  
  filter(PID!="Control") %>% 
  filter(QC.flag!="Duplicate= MICRO") %>%  
  select(c("sampleid", "Sample.Name", "PID")) %>%  
  dplyr::rename("chear_pid"="PID")

# CHEAR's samples 
microchear <- read.csv("micheals_epi_table_filtered_2.csv")
microchear.sub <- microchear %>%  
  dplyr::rename("chear_pid"="pid") %>%  
  select(c("chear_pid"))
microchear.sub$chear_pid <- as.character(microchear.sub$chear_pid )
  
# CHEAR's metadata
epichear <- read.csv("1977_Stool_Analytic_Data_082421.csv")
epichear.sub <- epichear %>%  
  dplyr::rename("chear_pid"="pid")
epichear.sub$chear_pid <- as.character(epichear.sub$chear_pid)

# Merge 
metadata <- microchear.sub %>%  
  left_join(manifest.sub, by="chear_pid") %>%  
  left_join(epichear.sub, by="chear_pid") %>%  
  relocate(sampleid, .before = chear_pid)

# Split by timepoint 
aim2b1 <- filter(metadata, timepoint=='B1' & !is.na(log_bp3))
aim2b4 <- filter(metadata, timepoint=='B4' & !is.na(log_bp3))
aim2pm1 <- filter(metadata, timepoint=='1PM' & !is.na(log_bp3))

# Save 
write.table(aim2b1, "C:/Users/laray/Box/Karin's Lab/GOCS/CHEAR/Final Reports and Code/Rerun/Statistical_Analysis_Datasets/microbiome/Aim2MetadataFinalB1.txt", row.names = FALSE, sep="\t", quote = FALSE, col.names=TRUE)
write.table(aim2b4, "C:/Users/laray/Box/Karin's Lab/GOCS/CHEAR/Final Reports and Code/Rerun/Statistical_Analysis_Datasets/microbiome/Aim2MetadataFinalB4.txt", row.names = FALSE, sep="\t", quote = FALSE, col.names=TRUE)
write.table(aim2pm1, "C:/Users/laray/Box/Karin's Lab/GOCS/CHEAR/Final Reports and Code/Rerun/Statistical_Analysis_Datasets/microbiome/Aim2MetadataFinalPM1.txt", row.names = FALSE, sep="\t", quote = FALSE, col.names=TRUE)

# Refer to locations of files  
metab1 <- here("C:/Users/laray/Box/Karin's Lab/GOCS/CHEAR/Final Reports and Code/Rerun/Statistical_Analysis_Datasets/microbiome/Aim2MetadataFinalB1.txt")
metab4 <- here("C:/Users/laray/Box/Karin's Lab/GOCS/CHEAR/Final Reports and Code/Rerun/Statistical_Analysis_Datasets/microbiome/Aim2MetadataFinalB4.txt")
metapm1 <- here("C:/Users/laray/Box/Karin's Lab/GOCS/CHEAR/Final Reports and Code/Rerun/Statistical_Analysis_Datasets/microbiome/Aim2MetadataFinalPM1.txt")

# ASV/Tax tables (refer to locations of files) 
#############################################
feat <- here("C:/Users/laray/Box/Karin's Lab/GOCS/CHEAR/Final Reports and Code/Rerun/Statistical_Analysis_Datasets/microbiome", "tablesilva.qza")
taxo <- here("C:/Users/laray/Box/Karin's Lab/GOCS/CHEAR/Final Reports and Code/Rerun/Statistical_Analysis_Datasets/microbiome", "silva_97_taxonomy.qza" )

#############################################
##                                         ##
##       Creating phyloseq objects         ##
##                                         ##
#############################################

# convert imported CHEAR files to phyloseq object 
pob1<-qza_to_phyloseq(features=feat,
                      taxonomy=taxo,
                      metadata=metab1)
pob4<-qza_to_phyloseq(features=feat,
                      taxonomy=taxo,
                      metadata=metab4)
popm1<-qza_to_phyloseq(features=feat,
                       taxonomy=taxo,
                       metadata=metapm1)

pob1
# otu_table()   OTU Table:         [ 12032 taxa and 201 samples ]
# sample_data() Sample Data:       [ 201 samples by 72 sample variables ]
# tax_table()   Taxonomy Table:    [ 12032 taxa by 7 taxonomic ranks ]

pob4
# otu_table()   OTU Table:         [ 12032 taxa and 239 samples ]
# sample_data() Sample Data:       [ 239 samples by 72 sample variables ]
# tax_table()   Taxonomy Table:    [ 12032 taxa by 7 taxonomic ranks ]

popm1
# otu_table()   OTU Table:         [ 12032 taxa and 162 samples ]
# sample_data() Sample Data:       [ 162 samples by 72 sample variables ]
# tax_table()   Taxonomy Table:    [ 12032 taxa by 7 taxonomic ranks ]


#############################################
##                                         ##
##       Filtering                         ##
##                                         ##
#############################################

### TAXANOMIC FILTERING ###
# prune ASVs that are not present in at least 2 samples [JJacobs suggestion]
pob1.a <- prune_taxa(taxa_sums(pob1) >= 2, pob1) # taxa_sums returns the total number of individuals observed from each species/taxa/ASV
pob1.a # [ 9209 taxa and 201 samples ]
pob4.a <- prune_taxa(taxa_sums(pob4) >= 2, pob4) 
pob4.a #  [ 10246 taxa and 239 samples ]
popm1.a <- prune_taxa(taxa_sums(popm1) >= 2, popm1) 
popm1.a # [ 8018 taxa and 162 samples ]


# What phyla are represented? 
get_taxa_unique(pob1.a, "Phylum")
# [1] "Proteobacteria"      "Actinobacteria"      "Firmicutes"          "Verrucomicrobia"     "Bacteroidetes"       "Euryarchaeota"       "Tenericutes"        
# [8] "Patescibacteria"     "Cyanobacteria"       NA                    "Fusobacteria"        "Spirochaetes"        "Elusimicrobia"       "Acidobacteria"      
# [15] "Synergistetes"       "Lentisphaerae"       "Epsilonbacteraeota"  "Chloroflexi"         "Gemmatimonadetes"    "Planctomycetes"      "Kiritimatiellaeota" 
# [22] "Deinococcus-Thermus"  
get_taxa_unique(pob4.a, "Phylum") #same
get_taxa_unique(popm1.a, "Phylum") #same


# READ COUNTS
# Create table, number of features for each phyla 
tableb1.features_phyla <- as.data.frame(table(tax_table(pob1.a)[,"Phylum"], exclude=NULL))
view(tableb1.features_phyla)
tableb4.features_phyla <- as.data.frame(table(tax_table(pob4.a)[,"Phylum"], exclude=NULL))
view(tableb4.features_phyla) 
tablepm1.features_phyla <- as.data.frame(table(tax_table(popm1.a)[,"Phylum"], exclude=NULL))
view(tablepm1.features_phyla) 
#Majority are Firmicutes & Bacteroidetes; there are 700-951 with missing phylum


# Remove NA and ambiguous phylum annotation 
pob1.b <- subset_taxa(pob1.a, !is.na(Phylum) & !Phylum %in% c("", "uncharacterized"))
pob4.b <- subset_taxa(pob4.a, !is.na(Phylum) & !Phylum %in% c("", "uncharacterized"))
popm1.b <- subset_taxa(popm1.a, !is.na(Phylum) & !Phylum %in% c("", "uncharacterized"))

tableb1.features_phyla2 <- as.data.frame(table(tax_table(pob1.b)[,"Phylum"], exclude=NULL))
tableb4.features_phyla2 <- as.data.frame(table(tax_table(pob4.b)[,"Phylum"], exclude=NULL))
tablepm1.features_phyla2 <- as.data.frame(table(tax_table(popm1.b)[,"Phylum"], exclude=NULL))

pob1.b #[ 8359 taxa and 201 samples ]
pob4.b #[ 9295 taxa and 239 samples ]
popm1.b #[ 7318 taxa and 162 samples ]


# drop samples below a read count of 10,000 [JJ suggestion]; shouldn't change anything here 
pob1.c <- prune_samples(sample_sums(pob1.b) >= 10000, pob1.b)
pob4.c <- prune_samples(sample_sums(pob4.b) >= 10000, pob4.b)
popm1.c <- prune_samples(sample_sums(popm1.b) >= 10000, popm1.b)


#############################################
##                                         ##
##       Data clean up                     ##
##                                         ##
#############################################

# EDCs: No changes necessary 

# Age: Scale 
sample_data(pob1.c)$age_c <- scale(sample_data(pob1.c)$Age, scale=FALSE)
sample_data(pob4.c)$age_c <- scale(sample_data(pob4.c)$Age, scale=FALSE)
sample_data(popm1.c)$age_c <- scale(sample_data(popm1.c)$Age, scale=FALSE)

# Creatinine: Scale 
sample_data(pob1.c)$cre_c <- scale(sample_data(pob1.c)$log_cre, scale=FALSE)
sample_data(pob4.c)$cre_c <- scale(sample_data(pob4.c)$log_cre, scale=FALSE)
sample_data(popm1.c)$cre_c <- scale(sample_data(popm1.c)$log_cre, scale=FALSE)

# Maternal education
sample_data(pob1.c)$medu = as_factor(sample_data(pob1.c)$medu)
sample_data(pob4.c)$medu = as_factor(sample_data(pob4.c)$medu)
sample_data(popm1.c)$medu = as_factor(sample_data(popm1.c)$medu)

sample_data(pob1.c)$medu <- factor(sample_data(pob1.c)$medu, labels=c("Secondary school or less", "More than secondary school"))
sample_data(pob4.c)$medu <- factor(sample_data(pob4.c)$medu, labels=c("Secondary school or less", "More than secondary school"))
sample_data(popm1.c)$medu <- factor(sample_data(popm1.c)$medu, labels=c("Secondary school or less", "More than secondary school"))

# Birth mode 
sample_data(pob1.c)$birth_mode = as_factor(sample_data(pob1.c)$birth_mode)
sample_data(pob4.c)$birth_mode = as_factor(sample_data(pob4.c)$birth_mode)
sample_data(popm1.c)$birth_mode = as_factor(sample_data(popm1.c)$birth_mode)

sample_data(pob1.c)$birth_mode <- factor(sample_data(pob1.c)$birth_mode, labels=c("Cesarean", "Vaginal"))
sample_data(pob4.c)$birth_mode <- factor(sample_data(pob4.c)$birth_mode, labels=c("Cesarean", "Vaginal"))
sample_data(popm1.c)$birth_mode <- factor(sample_data(popm1.c)$birth_mode, labels=c("Cesarean", "Vaginal"))

# Breast feeding
sample_data(pob1.c)$bfeedcat = as_factor(sample_data(pob1.c)$bfeedcat)
sample_data(pob4.c)$bfeedcat = as_factor(sample_data(pob4.c)$bfeedcat)
sample_data(popm1.c)$bfeedcat = as_factor(sample_data(popm1.c)$bfeedcat)

sample_data(pob1.c)$bfeedcat <- factor(sample_data(pob1.c)$bfeedcat, labels=c("<3 months", "3-6 months", ">6 months"))
sample_data(pob4.c)$bfeedcat <- factor(sample_data(pob4.c)$bfeedcat, labels=c("<3 months", "3-6 months", ">6 months"))
sample_data(popm1.c)$bfeedcat <- factor(sample_data(popm1.c)$bfeedcat, labels=c("<3 months", "3-6 months", ">6 months"))

# Average calories: Scale
sample_data(pob1.c)$avgcal_c <- scale(sample_data(pob1.c)$AvgCal, scale=FALSE)
sample_data(pob4.c)$avgcal_c <- scale(sample_data(pob4.c)$AvgCal, scale=FALSE)
sample_data(popm1.c)$avgcal_c <- scale(sample_data(popm1.c)$AvgCal, scale=FALSE)

# Antibiotic use 
sample_data(pob1.c)$patb_5_depo = as_factor(sample_data(pob1.c)$patb_5_depo)
sample_data(pob4.c)$patb_5_depo = as_factor(sample_data(pob4.c)$patb_5_depo)
sample_data(popm1.c)$patb_5_depo = as_factor(sample_data(popm1.c)$patb_5_depo)

sample_data(pob1.c)$patb_5_depo <- factor(sample_data(pob1.c)$patb_5_depo, labels=c("No", "Yes", "Missing/unknown"))
sample_data(pob4.c)$patb_5_depo <- factor(sample_data(pob4.c)$patb_5_depo, labels=c("No", "Yes", "Missing/unknown"))
sample_data(popm1.c)$patb_5_depo <- factor(sample_data(popm1.c)$patb_5_depo, labels=c("No", "Yes", "Missing/unknown"))

# Lab 
sample_data(pob1.c)$Lab = as_factor(sample_data(pob1.c)$Lab)
sample_data(pob4.c)$Lab = as_factor(sample_data(pob4.c)$Lab)
sample_data(popm1.c)$Lab = as_factor(sample_data(popm1.c)$Lab)

#############################################
##                                         ##
##      PREP FOR Maaslin2                  ##
##                                         ##
#############################################

# agglomerate to species level to control for duplicate taxa
# pob1.d <- tax_glom(pob1.c, taxrank="Species", NArm=TRUE)
# pob4.d <- tax_glom(pob4.c, taxrank="Species", NArm=TRUE)
# popm1.d <- tax_glom(popm1.c, taxrank="Species", NArm=TRUE)
# 
# ntaxa(pob1.c); ntaxa(pob1.d) #8359 taxa; 622
# ntaxa(pob4.c); ntaxa(pob4.d) #9295; 632
# ntaxa(popm1.c); ntaxa(popm1.d) #7318; 590

# MODIFICATION: agglomerate to genus
pob1.d <- tax_glom(pob1.c, taxrank="Genus", NArm=TRUE)
pob4.d <- tax_glom(pob4.c, taxrank="Genus", NArm=TRUE)
popm1.d <- tax_glom(popm1.c, taxrank="Genus", NArm=TRUE)

ntaxa(pob1.c); ntaxa(pob1.d) #8359 taxa; 354
ntaxa(pob4.c); ntaxa(pob4.d) #9295; 358
ntaxa(popm1.c); ntaxa(popm1.d) #7318; 347


# define individual objects as data frames
taxtableb1<-as.data.frame(tax_table(pob1.d)) 
asvtableb1<-as.data.frame(otu_table(pob1.d))
metatableb1<-meta(pob1.d) 

taxtableb4<-as.data.frame(tax_table(pob4.d)) 
asvtableb4<-as.data.frame(otu_table(pob4.d))
metatableb4<-meta(pob4.d) 

taxtablepm1<-as.data.frame(tax_table(popm1.d)) 
asvtablepm1<-as.data.frame(otu_table(popm1.d))
metatablepm1<-meta(popm1.d) 

# Create new label
taxtableb1$taxa <- paste0(taxtableb1$Phylum, "; ", taxtableb1$Class, "; ", taxtableb1$Order, "; ", taxtableb1$Family, "; ", taxtableb1$Genus)
taxtableb4$taxa <- paste0(taxtableb4$Phylum, "; ", taxtableb4$Class, "; ", taxtableb4$Order, "; ", taxtableb4$Family, "; ", taxtableb4$Genus)
taxtablepm1$taxa <- paste0(taxtablepm1$Phylum, "; ", taxtablepm1$Class, "; ", taxtablepm1$Order, "; ", taxtablepm1$Family, "; ", taxtablepm1$Genus)

taxtableb1sub <- select(taxtableb1, c("taxa"))
taxtableb4sub <- select(taxtableb4, c("taxa"))
taxtablepm1sub <- select(taxtablepm1, c("taxa"))

# Transpose asv table (not required here)
# asvtableb1 <- t(asvtableb1)
# asvtableb4 <- t(asvtableb4)
# asvtablepm1 <- t(asvtablepm1)

# Merge tax and asv tables
asvnamesb1 <-  merge(taxtableb1sub, asvtableb1, by="row.names")
asvnamesb4 <-  merge(taxtableb4sub, asvtableb4, by="row.names")
asvnamespm1 <-  merge(taxtablepm1sub, asvtablepm1, by="row.names")

# Get rownames 
rownames(asvnamesb1) <- asvnamesb1$taxa
rownames(asvnamesb4) <- asvnamesb4$taxa
rownames(asvnamespm1) <- asvnamespm1$taxa

# Remove row.names and taxa columns
asvnamesb1 <- select(asvnamesb1, -c("Row.names", "taxa"))
asvnamesb4 <- select(asvnamesb4, -c("Row.names", "taxa"))
asvnamespm1 <- select(asvnamespm1, -c("Row.names", "taxa"))

# Transpose to samples as rows
asvnamesb1t <- as.data.frame(t(asvnamesb1)) 
asvnamesb4t <- as.data.frame(t(asvnamesb4)) 
asvnamespm1t <- as.data.frame(t(asvnamespm1)) 


#############################################
##                                         ##
##       MaAsLin2                          ##
##                                         ##
#############################################

# Running maaslin based on significant associations in prior analyses 
getwd()
setwd("C:/Users/laray/Box/Projects/2020-Dissertation/Analysis")

# Check 
names(metatableb1)
levels(metatableb1$patb_5_depo)

# B1:
############################################
edcb1 <- names(metatableb1[15:33])

for(i in edcb1){
  Maaslin2(input_data=asvnamesb1t, input_metadata=metatableb1, 
           output= paste0("Output/Aim2/09082021/maasNArm/maas_b1_",i), 
           fixed_effects=c(paste0(i),"age_c", "cre_c", "avgcal_c", "bfeedcat", "birth_mode", "medu", "patb_5_depo", "Lab"), 
           reference=c("bfeedcat, 3-6 months", "birth_mode, Vaginal", "medu, Secondary school or less", "patb_5_depo, No", "Lab, CHEAR"))
}

# B4:  
############################################
edcb4 <- names(metatableb4[15:33])

for(i in edcb4){
  Maaslin2(input_data=asvnamesb4t, input_metadata=metatableb4, 
           output= paste0("Output/Aim2/09082021/maasNArm/maas_b4_",i), 
           fixed_effects=c(paste0(i),"age_c", "cre_c", "avgcal_c", "bfeedcat", "birth_mode", "medu", "patb_5_depo", "Lab"), 
           reference=c("bfeedcat, 3-6 months", "birth_mode, Vaginal", "medu, Secondary school or less", "patb_5_depo, No", "Lab, CHEAR"))
}

# PM1: 
############################################
edcpm1 <- names(metatablepm1[15:33])

for(i in edcpm1){
  Maaslin2(input_data=asvnamespm1t, input_metadata=metatablepm1, 
           output= paste0("Output/Aim2/09082021/maasNArm/maas_pm1_",i), 
           fixed_effects=c(paste0(i),"age_c", "cre_c", "avgcal_c", "bfeedcat", "birth_mode", "medu", "patb_5_depo"), 
           reference=c("bfeedcat, 3-6 months", "birth_mode, Vaginal", "medu, Secondary school or less", "patb_5_depo, No"))
}

