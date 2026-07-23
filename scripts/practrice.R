# Practice RScript
# Demonstrates the capacity to apply core concepts for bioinformatics

## Authors: Acosta, Joe A
## Date: 07/23/2026

## Prepare the Environment:
# set working directory
setwd("C:/Users/Drago/source/repos/rna_seq_intro/scripts")
getwd()

# load libraries
library(ggplot2)
library(ggthemes)
library(tidyverse)
library(GEOquery)

## Load and Prepare Dataset:
# load dataset
dat <- read.csv("../data/GSE183947_fpkm.csv")

head(dat)
dim(dat)

# load metadata
gse <- getGEO(GEO = "GSE183947", GSEMatrix = TRUE)
metadata <- pData(phenoData(gse[[1]]))

head(metadata)
dim(metadata)

# merge metadata and dataset
dat.long <- dat |>
  rename(gene = X) |>
  gather(key = "sample", value = "fpkm", -gene)

metadata.modified <- metadata |>
  select(1, 10, 11, 17) |>
  rename(sample = description) |>
  rename(tissue = characteristics_ch1) |>
  rename(metastasis = characteristics_ch1.1) |>
  mutate(tissue = gsub("tissue: ", "", tissue)) |>
  mutate(metastasis = gsub("metastasis: ", "", metastasis))

dat.long <- dat.long |>
  left_join(metadata.modified, by = c("sample" = "sample"))

glimpse(dat.long)
dim(dat.long)

## Conduct Exploratory Data Analysis:

# Scatter: How are BRCA1 and BRCA2 gene expressions related
dat.long |>
  filter(gene == "BRCA1" | gene == "BRCA2") |>
  spread(key = "gene", value = "fpkm") |>
  ggplot(aes(x = BRCA1, y = BRCA2, color = tissue)) +
  geom_point() +
  geom_smooth(method = "lm", se = FALSE) +
  labs(
    title = "BRCA1 vs BRCA2 Gene Expression",
    subtitle = "Comparing gene expression in tumor and normal cells",
    x = "BRCA1 (fpkm)",
    y = "BRCA2 (fpkm)",
  )

# ggsave("brca1_brca2_scatter.png", path = "../outputs", width = 10, height = 8)

# In normal tissue, BRCA1 and BRCA2 have a positive correlation. In tumor tissue
# the correlation is negative

# Evaluate the data given 

dat.long |>
  filter(gene == "BRCA1" | gene == "BRCA2" | gene == 'MYCN') |>
  spread(key = "gene", value = "fpkm") |>
  ggplot(aes(x = BRCA1, y = BRCA2, colour = MYCN)) +
  geom_point() +
  facet_wrap(~tissue)

# a bit inconclusive to look at the MYCN as well as the BRCA1 and BRCA2  


# Density and Histogram 
# explores how BRCA1 expression and BRCA2 expression differ among normal and
# tumor cells
dat.long |>
  filter(gene == "BRCA1") |>
  ggplot(aes(x = fpkm, fill = tissue)) +
  geom_density(alpha = .6)

dat.long |>
  filter(gene == "BRCA2") |>
  ggplot(aes(x = fpkm, fill = tissue)) +
  geom_histogram(binwidth = .5) +
  facet_wrap(~tissue)

# Box and Violin
# plots the gene expression for BRCA1
dat.long |>
  filter(gene == "BRCA1") |>
  ggplot(aes(y = fpkm, fill = tissue)) +
  geom_boxplot()

dat.long |>
  filter(gene == "BRCA1") |>
  ggplot(aes(x = tissue, y = fpkm, fill = tissue)) +
  geom_violin() +
  geom_boxplot(width = .15)

# Bar
dat.long |>
  filter(gene == "BRCA1") |>
  ggplot(aes(x = sample, y = fpkm, fill = tissue)) +
  geom_col()

# Heat

genes.of.interest <- c('BRCA1', 'BRCA2', 'TP53', 'ALK', 'MYCN')
# used the same gene list as the tutorial since, at the present moment, I
# haven't developed the skills or expertise in breast cancer bioinformatics
# needed to extract this information sans a tutorial

dat.long |>
  filter(gene %in% genes.of.interest) |>
  ggplot(aes(x = sample, y = gene, fill = fpkm)) +
  geom_tile() +
  scale_fill_gradient(low = "white", high = "red")
