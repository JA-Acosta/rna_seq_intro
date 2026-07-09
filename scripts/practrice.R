# Practice RScript
# Demonstrates the capacity to apply core concepts for bioinformatics

## Authors: Acosta, Joe A
## Date: 07/09/2026

## Sets Working Directory
setwd("C:/Users/Drago/source/repos/rna_seq_intro/scripts")

## Loads Librarires
library(tidyverse)
library(GEOquery)

## Load Dataset
dat <- read.csv("../data/GSE183947_fpkm.csv")
head(dat)

## Load Metadata
gse <- getGEO(GEO="GSE183947", GSEMatrix=TRUE) # didn't remember
metadata <- pData(phenoData(gse[[1]]))
head(metadata)

## Extract Vital Fields From Metadata
metadata.subset <- select(metadata, c("title", "characteristics_ch1", "characteristics_ch1.1", "description"))

## Create a Modified Metadata Subset 
metadata.modified <- metadata %>%
  select(c(1, 10, 11, 17)) %>%
  rename(metastatis = characteristics_ch1.1) %>% # didn't remember rename
  rename(tissue = characteristics_ch1) %>%
  mutate(tissue = gsub("tissue: ", "", tissue)) %>% # didn't remember mutate
  mutate(metastatis = gsub("metastasis: ", "", metastatis))

## Create a Long Daataset
dat.long <- dat %>%
  rename(gene = X) %>%
  gather(key="samples", value="FPKM", -gene) # didn't remember gather

## Merge Data with Metadata
dat.long <- dat.long %>%
  left_join(., metadata.modified, by=c("samples" = "description" )) # didn't remember join

## Final Dataset
head(dat.long)

# Questions Developed over course of practice:
# 1. What are the different types of joins
# 2. Why make the data a long dataset
# 3. Look up gsub


