# RNA-Seq Introduction Notebooks 
# Authors: Acosta, Joe
# Date: 07/07/2026

# ---------------------------------------

# Linking github repository
# usethis::create_from_github(
#   "https://github.com/JA-Acosta/rna_seq_intro.git",
#   destdir = "C:/Users/Drago/source/repos"
# )

# Data set info :
#   - NCBI GEO Data set found w/accession ID: 
#   - FPKM Normalized
#   - Don't forget gunzip
#   _ Added data set to .gitignore

# ---------------------------------------

# Read the gene expression data
# Retrieve the Metadata (or clinical data)
# Learn useful packages to join data, tidy data, or manipulate data.
# R Libraries Learned:
#   - dplyr
#   - tidyverse
#   - GEOquery





# script to manipulate gene expression data
# setwd("C:/Users/Drago/source/repos/rna_seq_intro/scripts")

# @Note: always place the setwd function at the top of the script. It allows
# future you to know where the R-Script directory is and easily set it


# Load Libraries
library(dplyr)
library(tidyverse)
library(GEOquery)

# @Note: If needed install the packages w/ install.packages("")
#   - First download RTools; will need to know the R version
#   - R.version.string; "R version 4.6.0 (2026-04-24 ucrt)"
#   - Downloaded Rtools45 (supports version 4.6.0) on CRAN


# Read in the Data
dat <- read.csv(file = "./data/GSE183947_fpkm.csv")
dim(dat)

# @Note: Since the data set is the raw data for each sample, we
# don't know which sample is cancerous and which is normal.The
# meta data provides us this information

# Get the Meta Data
gse <- getGEO(GEO = "GSE183947", GSEMatrix = TRUE)
metadata <- pData(phenoData(gse[[1]]))
head(metadata)

# Extract Relevant Meta Data columns
metadata.subset <- select(metadata, c(1, 10, 11, 17))

# rename is used to rename the columns
# mutate is used to create new cols, update cols, 
# Using the pipe operator in Tidyverse
metadata.modified <- metadata %>%
  select(1, 10, 11, 17) %>%
  rename(tissue = characteristics_ch1) %>%
  rename(metastasis = characteristics_ch1.1) %>%
  mutate(tissue = gsub("tissue: ", "", tissue)) %>%
  mutate(metastasis  = gsub("metastasis: ", "", metastasis))

# always covert the data to the long format since it is easier
# to work with

# reshaping data
# rename gene to X
# use gather() to convert from wide to long format
dat.long <- dat %>%
  rename(gene = X) %>%
  gather(key = 'samples', value = 'FPKM', -gene)


# add metadata to gene expression data
# join dataframes = dat.long + metadata.modified
dat.long <- dat.long %>%
  left_join(., metadata.modified, by = c("samples" = "description"))


head(dat.long)

# explore data
dat.long %>%
  filter(gene == 'BRCA1' | gene == 'BRCA2') %>%
  group_by(gene, tissue) %>%
  summarise(mean_FPKM = mean(FPKM),
            median_FPKM = median(FPKM)) %>%
  arrange(mean_FPKM)


# can calc additional iqr, upper, std_dev, more
# arrange will show in ascending but descending is -col_name


