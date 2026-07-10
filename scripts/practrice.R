# Practice RScript
# Demonstrates the capacity to apply core concepts for bioinformatics

## Authors: Acosta, Joe A
## Date: 07/09/2026


# get and set current wd
getwd() # check current directory
setwd("C:/Users/Drago/source/repos/rna_seq_intro/scripts")

getwd() # make sure the directory changes


# load libraries
library(tidyverse)
library(GEOquery)

## loads all of the working data required for the project

# load data set
dat <- read.csv("../data/GSE183947_fpkm.csv")
head(dat)
dim(dat)


# load metadata
gse = getGEO(GEO = "GSE183947", GSEMatrix = TRUE)
metadata <- pData(phenoData(gse[[1]]))


# modified metadata set
# -  get the specific columns I need
# -  clean cols
# -  check the output
metadata.modified <- metadata %>%
  select(c(1, 10, 11, 17)) %>%
  rename(tissue = characteristics_ch1) %>%
  rename(metastasis = characteristics_ch1.1) %>%
  rename(sample = description) %>%
  mutate(tissue = gsub("tissue: ", "", tissue)) %>%
  mutate(metastasis = gsub("metastasis: ", "", metastasis))

head(metadata.modified)
dim(metadata.modified)


# data set conversion to long format
dat.long <- dat %>%
  rename(gene = X) %>%
  gather(key = "sample", value = "FPKM", -gene)

head(dat.long)
dim(dat.long)


# merge long data with modified metadata
dat.long <- dat.long %>%
  left_join(., metadata.modified, by = c("sample" = "sample"))

head(dat.long)
dim(dat.long)


## Exploratory data analysis

# gets the relevant genes
# groups them by gene and tissue
# gets mean, median, std
eda <- dat.long %>%
  filter(gene == "BRCA1" | gene == "BRCA2") %>%
  group_by(gene, tissue) %>%
  summarise(
    mean_FPKM = mean(FPKM),
    median_FPKM = median(FPKM),
    std_FPKM = sd(FPKM)
  )

# rearrange based on mean in ascending order
arrange(eda, mean_FPKM)
