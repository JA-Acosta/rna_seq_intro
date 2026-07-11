# Practice RScript
# Demonstrates the capacity to apply core concepts for bioinformatics

## Authors: Acosta, Joe A
## Date: 07/10/2026

# get and set working directory
getwd()
setwd("C:/Users/Drago/source/repos/rna_seq_intro/scripts")

getwd()


# load libraries
library(tidyverse)
library(GEOquery)


# load dataset
dat <- read.csv("C:/Users/Drago/source/repos/rna_seq_intro/data/GSE183947_fpkm.csv")

head(dat)
dim(dat)


# load metadata
gse <- getGEO(GEO = "GSE183947", GSEMatrix = TRUE)
metadata <- pData(phenoData(gse[[1]]))

head(metadata)
dim(metadata)


# create a cleaned metadata subset
metadata.modified <- metadata %>%
  select(c(1, 10, 11, 17)) %>% # forgot at first; recalled without refering to notes
  rename(tissue = characteristics_ch1 ) %>%
  rename(metastasis = characteristics_ch1.1 ) %>%
  rename(sample = description) %>%
  mutate(tissue = gsub("tissue: ", "", tissue)) %>%
  mutate(metastasis = gsub("metastasis: ", "", metastasis))

head(metadata.modified)
dim(metadata.modified)


# convert dataset to long form
dat.long <- dat %>%
  rename(gene = X) %>%
  gather(key = "sample", value = "FPKM", -gene) # forgot at first; recalled without refering to notes

head(dat.long)
dim(dat.long)


# join metadata to dataset
dat.long <- dat.long %>%
  left_join(., metadata.modified, by = c("sample" = "sample"))

head(dat.long)
dim(dat.long)


# conduct basic exploratory data analysis on dataset
eda <- dat.long %>%
  filter(gene == "BRCA1" | gene == "BRCA2") %>%
  group_by(gene, tissue) %>%
  summarise(
    mean_FPKM = mean(FPKM),
    median_FPKM = median(FPKM)
  )

# arrange by median descending
arrange(eda, desc(median_FPKM))




