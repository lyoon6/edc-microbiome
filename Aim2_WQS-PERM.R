#################################################
#   CHEAR 1977 Microbiome Aim (Aim 2)           #
#                                               #
#    Writted by Lara Yoon                       #
#     September 09, 2021                        #
#                                               #
#    Sensitivity analysis:                      #
#      Regress WQS on other things              #
#################################################


## Load libraries
library(pacman)
pacman::p_load(here, tidyr, gtools, tidyverse, broom, data.table, ggplot2, sjPlot, gWQS, kableExtra, vegan, phyloseq, qiime2R)
pacman::p_loaded()


#############################################
##                                         ##
##       Loading the Data Sets             ##
##                                         ##
#############################################

setwd('C:/Users/laray/Box/Projects/2020-Dissertation/Analysis')

# Microbiome files for import 
feat <- here("Microbiome/Input", "tablesilva.qza")
taxo <- here("Microbiome/Input", "silva_97_taxonomy.qza" )
treeroot <- here("Microbiome/Input", "rooted_tree_masked_alignment.qza")

# Metadata
metapm1 <- here("Data/FinalData/Aim2MetadataFinalPM1.txt")



#############################################
##                                         ##
##       Loading the Data Sets             ##
##                                         ##
#############################################
# convert imported CHEAR files to phyloseq objects (separate by time point) 
# note that phyloseq does not allow duplicate row.names, so we will have to do the analysis separate by time point 
popm1<-qza_to_phyloseq(features=feat,
                       taxonomy=taxo,
                       tree=treeroot,
                       metadata=metapm1)


#############################################
##                                         ##
##       Filtering                         ##
##                                         ##
#############################################

popm1
# otu_table()   OTU Table:         [ 6270 taxa and 165 samples ]
# sample_data() Sample Data:       [ 165 samples by 124 sample variables ]
# tax_table()   Taxonomy Table:    [ 6270 taxa by 7 taxonomic ranks ]
# phy_tree()    Phylogenetic Tree: [ 6270 tips and 6120 internal nodes ]


# Remove NA and ambiguous phylum annotation from each of the phyloseq objects
pobpm1 <- subset_taxa(pob4, !is.na(Phylum) & !Phylum %in% c("", "uncharacterized"))

# Examine phyla frequency 
table.features_phyla2 <- as.data.frame(table(tax_table(poall)[,"Phylum"], exclude=NULL))

# Examine prevalence of each feature
prevdf = apply(X=otu_table(poall), MARGIN = ifelse(taxa_are_rows(poall), yes=1, no=2), FUN=function(x) {sum(x>0)})

# Add taxonomy and total read counts
prevdf = data.frame(Prevalence = prevdf, TotalAbundance = taxa_sums(poall), tax_table(poall))

# Compute the total and average prevalence of the features in each phylum
prevdfSUMS <- plyr::ddply(prevdf, "Phylum", function(df1){cbind(mean(df1$Prevalence),sum(df1$Prevalence))})
names(prevdfSUMS)[2] <- "meanprev"
names(prevdfSUMS)[3] <- "sumprev"

# Filter out low prevalence phyla from phyloseq objects (prevalence <=1)
PhyForFilt <- c("Chlamydiae", "Kiritimatiellaeota", "Deinococcus-Thermus", "Deferribacteres")
subpopm1 <- subset_taxa(popm1, !Phylum %in% PhyForFilt)

# Subset to remaining phyla in dataframes
prevdf.A = subset(prevdf, Phylum %in% get_taxa_unique(poall, "Phylum"))

# Define prevalence threshold as 2% of total samples (.02*261=5.22)
prevthresh = 0.02 * nsamples(subpopm1) 

# Execute prevalence filter, using `prune_taxa()` function
keepTaxa = rownames(prevdf.A)[(prevdf.A$Prevalence >= prevthresh)] # can use this for all 3 tanner stages
subpopm1.b = prune_taxa(keepTaxa, subpopm1)
subpopm1.b
# otu_table()   OTU Table:         [ 1975 taxa and 165 samples ]
# sample_data() Sample Data:       [ 165 samples by 124 sample variables ]
# tax_table()   Taxonomy Table:    [ 1975 taxa by 7 taxonomic ranks ]
# phy_tree()    Phylogenetic Tree: [ 1975 tips and 1948 internal nodes ]

# prune those less than 10,000
subpopm1.b <- prune_samples(sample_sums(subpopm1.b) >= 10000, subpopm1.b)
subpopm1.b
# otu_table()   OTU Table:         [ 1975 taxa and 164 samples ]
# sample_data() Sample Data:       [ 164 samples by 124 sample variables ]
# tax_table()   Taxonomy Table:    [ 1975 taxa by 7 taxonomic ranks ]
# phy_tree()    Phylogenetic Tree: [ 1975 tips and 1948 internal nodes ]


#############################################
##                                         ##
##       PERMANOVA                         ##
##                                         ##
#############################################

subpopm1.b.rf

# rarefY 
subpopm1.b.rf <- rarefy_even_depth(subpopm1.b, sample.size=min(sample_sums(subpopm1.b)), rngseed=801, replace=F)

# calculate distances 
bray_dist_pm1 <- phyloseq::distance(subpopm1.b.rf, method="bray")

#### [TEMP FIX: replace NA in avgcal and abx with mean/median]
sample_data(subpopm1.b.rf)$avgcal[is.na(sample_data(subpopm1.b.rf)$avgcal)] <- mean(sample_data(subpopm1.b.rf)$avgcal, na.rm = TRUE)
sample_data(subpopm1.b.rf)$abx[is.na(sample_data(subpopm1.b.rf)$abx)] <- median(sample_data(subpopm1.b.rf)$abx, na.rm = TRUE)

# Drop missing 
# 1. Omit NAs for EDCs (they will all have the same NAs)
sample_data(subpopm1.b.rf) <- sample_data(subpopm1.b.rf)[!(is.na(sample_data(subpopm1.b.rf)$wqs)), ]
view(sample_data(subpopm1.b.rf))

# Run model 
permwqs <- adonis2(bray_dist_pm1 ~ sample_data(subpopm1.b.rf)$wqs + sample_data(subpopm1.b.rf)$log_cre + sample_data(subpopm1.b.rf)$fatp + sample_data(subpopm1.b.rf)$age 
                           + sample_data(subpopm1.b.rf)$medu + sample_data(subpopm1.b.rf)$avgcal + sample_data(subpopm1.b.rf)$abx 
                           + sample_data(subpopm1.b.rf)$bfeedcat + sample_data(subpopm1.b.rf)$birth_mode)
permwqs


