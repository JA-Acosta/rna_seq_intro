# Practice RScript
# Demonstrates the capacity to apply core concepts for bioinformatics

## Authors: Acosta, Joe A
## Date: 07/12/2026

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
metadata = pData(phenoData(gse[[1]]))

head(metadata)
dim(metadata)


# create metadata subset
metadata.modified <- metadata %>%
  select(1, 10, 11, 17) %>%
  rename(tissue = characteristics_ch1 ) %>%
  rename(metastasis = characteristics_ch1.1 ) %>%
  rename(sample = description) %>%
  mutate(tissue = gsub("tissue: ", "", tissue)) %>%
  mutate(metastasis = gsub("metastasis: ", "", metastasis))

head(metadata.modified)
dim(metadata.modified)


# merge metadata with dataset
dat.long <- dat %>%
  rename(gene = X) %>%
  gather(key = "sample", value = "FPKM", -gene)

head(dat.long)
dim(dat.long)

dat.long <- dat.long %>%
  left_join(., metadata.modified, by = c("sample" = "sample"))

head(dat.long)
dim(dat.long)

# do basic eda
# extract the mean, median, and std for both BRCA1 and BRCA2
dat.long %>%
  filter(gene == "BRCA1" | gene == "BRCA2") %>%
  group_by(gene) %>%
  summarise(
    mean_FPKM = mean(FPKM),
    median_FPKM = median(FPKM),
    std_FPKM = sd(FPKM)
  )


# create various plots for eda (bar, box, density, scatter, heat)

# Bar Plot
# Compare the expression of BRCA1 among all samples collected
dat.long %>%
  filter(gene == "BRCA1") %>%
  ggplot(., aes(x = sample, y = FPKM, fill = tissue)) +
  geom_col() # didn't remember column plot; looked at notes

ggsave("brca1_sample.png", path = "../outputs", width = 10, height = 8)
# didn't remember how to save a figure; looked at notes


# Box Plot
# Compare the impact of metastasis on gene expression of BRCA1
dat.long %>%
  filter(gene == "BRCA1") %>%
  ggplot(., aes(x = metastasis, y = FPKM, fill = metastasis)) +
  geom_boxplot()


# Density Plot
# looks at the expression of BRCA1 in tumor and normal tissue
dat.long %>%
  filter(gene == "BRCA1") %>%
  ggplot(., aes(x = FPKM, fill = tissue)) +
  geom_density(alpha = .6) # didn't remember to include a density plot at first


# Scatter Plot
# comparing the gene expression between BRCA1 and BRCA2
dat.long %>%
  filter(gene == "BRCA1" | gene == "BRCA2") %>%
  spread(key = "gene", value = "FPKM") %>% # didn't remember spread
  ggplot(., aes(x = BRCA1, y = BRCA2, colour = tissue)) + # didn't remember colours
  geom_point() +
  geom_smooth(method = "lm", se = FALSE) # didn't remember smooth

# Heat Map
genes.of.interest <- c('BRCA1', 'BRCA2', 'TP53', 'ALK', 'MYCN') # copied the list
# of genes since I don't have familiarity with them yet.

dat.long %>%
  filter(gene %in% genes.of.interest) %>%
  ggplot(., aes(x = sample, y = gene, fill = FPKM)) + # didn't remember args
  geom_tile() + # didn't remember tile
  scale_fill_gradient(low = 'white', high = 'red')# didn't remember scale_fill_gradient


# Questions

# How does spread works?

# How do we add labels to the plots?
# Titles, axis, etc.

# In ggplot what is the difference between colour and colours?


