# Practice RScript
# Demonstrates the capacity to apply core concepts for bioinformatics

## Authors: Acosta, Joe A
## Date: 07/09/2026

# Working Directory
setwd("C:/Users/Drago/source/repos/rna_seq_intro/scripts")

# Load Libraries
library(tidyverse)
library(GEOquery)

# Load Data
dat <- read.csv(file = "../data/GSE183947_fpkm.csv")

dim(dat)

# Load Metadata
gse <- getGEO(GEO = "GSE183947", GSEMatrix = TRUE) # didn't remember GSEMatrix
metadata <- pData(phenoData(gse[[1]]))

head(metadata)
colnames(metadata)

# Create Metadata Subset and Clean It
metadata.modified <- metadata %>%
  select(c(1, 10, 11, 17)) %>%
  rename("tissue" = characteristics_ch1) %>%
  rename("metastasis" = characteristics_ch1.1) %>%
  rename("sample" = description) %>%
  mutate(tissue = gsub("tissue: ", "", tissue)) %>%
  mutate(metastasis = gsub("metastasis: ", "", metastasis))

head(metadata.modified)
dim(metadata.modified)

# Create Long Dataset
dat.long <- dat %>%
  rename(gene = X) %>%
  gather(key = "sample", value = "FPKM", -gene) # didn't remember gather

# Merge Metadata and Long Dataset
dat.long <- dat.long %>%
  left_join(., metadata.modified, by = c("sample" = "sample"))

head(dat.long)
dim(dat.long)

# Data Exploration
dat.long %>%
  filter(gene == "BRCA1" | gene == "BRCA2") %>% # didn't remember filter
  group_by(gene, tissue) %>% # didn't remember groupby
  summarise( # didn't remembered summarise
    mean_FPKM = mean(FPKM),
    median_FPKM = median(FPKM)
  ) %>%
  arrange(mean_FPKM)

# First Collect the relevant genes then conduct data exploration on them.
# In this case, we are looking at the BRCA1 and BRCA2 genes 
