# Practice RScript
# Demonstrates the capacity to apply core concepts for bioinformatics

## Authors: Acosta, Joe A
## Date: 07/14/2026

# set working directory
setwd("C:/Users/Drago/source/repos/rna_seq_intro/scripts")
getwd()

# load libraries
library(tidyr)
library(tidyverse)
library(ggplot2)
library(GEOquery)

# load dataset
dat <- read.csv("../data/GSE183947_fpkm.csv")
head(dat)
dim(dat)

# load corresponding metadata
gse <- getGEO(GEO = "GSE183947", GSEMatrix = TRUE)
metadata <- pData(phenoData(gse[[1]]))
head(metadata)
dim(metadata)

# clean metadata
metadata.modified <- metadata %>%
  select(1, 10, 11, 17) %>%
  rename(sample = description) %>%
  rename(tissue = characteristics_ch1) %>%
  rename(metastasis = characteristics_ch1.1 ) %>%
  mutate(tissue = gsub("tissue: ", "", tissue)) %>%
  mutate(metastasis = gsub("metastasis: ", "", metastasis))
  
# merge cleaned metadata and dataset
## create long datset
dat.long <- dat %>%
  rename(gene = X) %>%
  gather(key = "sample", value = "FPKM", -gene)

## merge both datasets
dat.long <- dat.long %>%
left_join(., metadata.modified, by = c("sample" = "sample"))

head(dat.long)
dim(dat.long)

# conduct basic eda
dat.long %>%
  filter(gene == "BRCA1" | gene == "BRCA2") %>%
  group_by(gene, tissue) %>%
  summarise(
    mean_FPKM = mean(FPKM),
    median_FPKM = median(FPKM)
  ) %>% # arrange by median descending
  arrange(desc(mean_FPKM))

# conduct more in depth eda with plots
## Box plot: Compare the impact of metastasis on gene expression of BRCA1
dat.long %>%
  filter(gene == "BRCA1") %>%
  ggplot(., aes(x = metastasis, y = FPKM, fill = metastasis)) +
  geom_boxplot()

## save Box plot
ggsave(filename = "BRCA1_metastasis_box.png", path = "../outputs", width = 10, height = 8)

## Bar plot: Compares gene expression of BRCA1 across all samples collected
dat.long %>%
  filter(gene == "BRCA1") %>%
  ggplot(., aes(x = sample, y = FPKM, fill = tissue)) +
  geom_col()

## save Bar plot
ggsave(filename = "BRCA1_expression_bar.png", path = "../outputs", width = 10, height = 8)

## Density plot: Compares gene expression of BRCA1 among tumor and normal samples
dat.long %>%
  filter(gene == "BRCA1") %>%
  ggplot(., aes(x = FPKM, fill = tissue)) +
  geom_density(alpha = .6)

## save Density plot
ggsave(filename = "BRCA1_density.png", path = "../outputs", width = 10, height = 8)

## Scatter plot: Compares gene expression of BRCA1 and BRCA2 among normal and tumor samples
dat.long %>%
  filter(gene == "BRCA1" | gene == "BRCA2") %>%
  spread(key = gene, value = FPKM) %>%
  ggplot(., aes(x = BRCA1, y = BRCA2, colour = tissue)) +
  geom_point() +
  geom_smooth(method = "lm", se = FALSE)

# save Scatter plot
ggsave(filename = "BRCA1_BRCA2_scatter.png", path = "../outputs", width = 10, height = 8)

## Heat map

genes.of.interest <- c('BRCA1', 'BRCA2', 'TP53', 'ALK', 'MYCN')
# used the same gene list as the tutorial since, at the present moment, I
# haven't developed the skills or expertise in breast cancer bioinformatics
# needed to extract this information sans a tutorial

dat.long %>%
  filter(gene %in% genes.of.interest) %>%
  ggplot(., aes(x = sample, y = gene, fill = FPKM)) +
  geom_tile() +
  scale_fill_gradient(low = 'white', high = 'red')

ggsave(filename = "genes_heatmap.png", path = '../outputs', width = 10, height = 8)