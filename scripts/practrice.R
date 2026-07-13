# Practice RScript
# Demonstrates the capacity to apply core concepts for bioinformatics

## Authors: Acosta, Joe A
## Date: 07/13/2026

# set working directory
setwd("C:/Users/Drago/source/repos/rna_seq_intro/scripts")
getwd()

# load libraries
library(tidyverse)
library(GEOquery)
library(ggplot2)
library(dplyr)

# load dataset
dat <- read.csv("../data/GSE183947_fpkm.csv")

head(dat)
dim(dat)

# load metadata
gse <- getGEO(GEO = "GSE183947", GSEMatrix = TRUE)
metadata <- pData(phenoData(gse[[1]]))

head(metadata)
dim(metadata)

# merge dataset and relevant metadata
## create a long form of the dataset
dat.long <- dat %>%
  rename(gene = X) %>%
  gather(key = "sample", value = "FPKM", -gene)

head(dat.long)
dim(dat.long)

## create a modified metadata with all relevant data
metadata.modified <- metadata %>%
  select(1, 10, 11, 17) %>%
  rename(sample = description) %>%
  rename(tissue = characteristics_ch1 ) %>%
  rename(metastasis = characteristics_ch1.1 ) %>%
  mutate(tissue = gsub("tissue: ", "", tissue)) %>%
  mutate(metastasis = gsub("metastasis: ", "", metastasis))

head(metadata.modified)
dim(metadata.modified)

dat.long <- dat.long %>%
  left_join(., metadata.modified, by = c("sample" = "sample"))

head(dat.long)
dim(dat.long)

# conduct basic eda
## compare the mean, median, and std of BRCA1 and BRCA2 genes in tumor and normal
## samples

dat.long %>%
  filter(gene == "BRCA1" | gene == "BRCA2") %>%
  group_by(gene, tissue) %>%
  summarise(
    mean_FPKM = mean(FPKM),
    median_FPKM = median(FPKM),
    std_FPKM = sd(FPKM)
  ) %>%
  ## arrange by mean descending
  arrange(desc(mean_FPKM))

# create plots (Box, Bar, Density, Scatter, Heat)

## Box
## Compares BRCA1 FPKM among samples with metastasis and those without
dat.long %>%
  filter(gene == "BRCA1") %>%
  ggplot(., aes(x = metastasis, y = FPKM, fill = metastasis)) +
  geom_boxplot()

ggsave("BRCA1_metastasis_box.png", path = "../outputs", width = 10, height = 8)

## Bar
## Compares the gene expression across tumor and normal samples
dat.long %>%
  filter(gene == "BRCA1") %>%
  ggplot(., aes(x = sample, y = FPKM, fill = tissue)) +
  geom_col()

ggsave("BRCA1_FPKM_values.png", path = "../outputs", width = 10, height = 8)

## Density
## Sums expression of FPKM for BRCA1 gene across tumor and normal samples
## showing the ranges for FPRM given the sample type
dat.long %>%
  filter(gene == "BRCA1") %>%
  ggplot(., aes(x = FPKM, fill = tissue)) +
  geom_density(alpha = .6)

ggsave("BRCA1_FPKM_expression.png", path="../outputs", width = 10, height = 8)

## Scatter
## Compares the expression rates of BRCA1 and BRCA2 among tumor and normal cells
## to see how expression differs when a cell becomes cancerous
dat.long %>%
  filter(gene == "BRCA1" | gene == "BRCA2") %>%
  spread(key = gene, value = FPKM) %>% # still confused by how this works
  ggplot(., aes(x = BRCA1, y = BRCA2, colour = tissue)) +
  geom_point() +
  geom_smooth(method = "lm", se = FALSE)# didn't remember geom_smooth

ggsave("BRCA1_BRCA2_expression.png", path = "../outputs", width = 10, height = 8)

## Heat
genes.of.interest <- c('BRCA1', 'BRCA2', 'TP53', 'ALK', 'MYCN')
# used the same gene list as the tutorial since, at the present moment, I
# haven't developed the skills or expertise in breast cancer bioinformatics
# needed to extract this information sans a tutorial

dat.long %>%
  filter(gene %in% genes.of.interest) %>%
  ggplot(., aes(x = sample, y = gene, fill = FPKM)) +
  geom_tile() +
  scale_fill_gradient(low = "white", high = "red") # confused scale_fill with
  # scale_color gradient

ggsave("heat_genes_of_interest.png", path = "../outputs", height = 8, width = 10)

# conclusion
# managed to retain the fundamentals required to plot. Next step will be to learn
# the fundamentals of breast cancers and what these values mean. Started this by
# developing an understanding of FPKM and where it comes from and how it compares
# to other normalization strategies. Though still a bit confused, will review.