# Practice RScript
# Demonstrates the capacity to apply core concepts for bioinformatics

## Authors: Acosta, Joe A
## Date: 07/22/2026

# Loading Environment:
# set the working directory
setwd("C:/Users/Drago/source/repos/rna_seq_intro/scripts")
getwd()

# import libraries
library(tidyverse)
library(GEOquery)
library(ggplot2)
library(ggthemes)

# Loading and Preparing Dataset:
# import data
dat <- read.csv('../data/GSE183947_fpkm.csv')

head(dat)
dim(dat)

# import metadata
gse <- getGEO(GEO = 'GSE183947', GSEMatrix = TRUE)
metadata <- pData(phenoData(gse[[1]]))

head(metadata)
dim(metadata)

# clean and merge data and metadata
# used R built in pipe instead of tidyverse; recently learned in r4ds
dat.long <- dat |>
  rename(gene = X) |>
  gather(key = "sample", value = "fpkm", -gene)

glimpse(dat.long) # used instead of head to practice; recently learned in r4ds
dim(dat.long)

metadata.modified <- metadata |>
  select(1, 10, 11, 17) |>
  rename(tissue = characteristics_ch1) |>
  rename(metastasis = characteristics_ch1.1) |>
  rename(sample = description) |>
  mutate(tissue = gsub("tissue: ", "", tissue)) |>
  mutate(metastasis = gsub("metastasis: ", "", metastasis))

dat.long <- dat.long %>%
  left_join(., metadata.modified, by = c("sample" = "sample"))

glimpse(dat.long)
dim(dat.long)

# do simple EDA
dat.long %>%
  filter(gene == 'BRCA1' | gene == 'BRCA2') %>%
  group_by(gene, tissue) %>%
  summarise(
    mean_fpkm = mean(fpkm),
    median_fpkm = median(fpkm)
  ) %>%
  arrange(., desc(median_fpkm))

# EDA Visualization:
# Scatter Plot: identify the relationships between BRCA1 and BRCA2
dat.long %>%
  filter(gene == "BRCA1" | gene == "BRCA2") %>%
  spread(., key = "gene", value = "fpkm") %>%
  ggplot(., aes(x = BRCA1, y = BRCA2, colour = tissue, shape = tissue)) +
  geom_point() +
  labs(
    title = "Relationship between BRCA1 and BRCA2",
    subtitle = "quantifies the relationship among both genes in cancerous and normal tissue",
    x = "BRCA1 (fpkm)",
    y = "BRCA2 (fpkm)"
  ) +
  scale_color_colorblind()

ggsave("brca1_brca2_relationship.png", path = "../outputs", width = 10, height = 8)

# Box Plot: identify the relationship between BRCA1 expression and metastatic tissue
dat.long %>%
  filter(gene == "BRCA1") %>%
  ggplot(., aes(y = fpkm, fill = metastasis)) +
  geom_boxplot() +
  labs(
    title = "BRCA1 Expression given cell Metastasis",
    x = "Metastasis (no/yes)",
    y = "BRCA1 (fpkm)"
  )

ggsave("brca1_metastasis.png", path = "../outputs", width = 10, height = 8)
  
# Bar Plot: visualize the expression rates for BRCA1 across all individuals
dat.long %>%
  filter(gene == "BRCA1") %>%
  ggplot(., aes(x = sample, y = fpkm, fill = tissue)) +
  geom_col() +
  labs(
    title = "FPKM Expressions Among Tissues",
    x = "Samples",
    y = "FPKM"
  )

ggsave(filename = "brca1_tissue.png", path = "../outputs", width = 10, height = 8)

# Distribution Plot
dat.long %>%
  filter(gene == "BRCA1") %>%
  ggplot(., aes(x = fpkm, fill = tissue)) +
  geom_density(alpha = .6) +
  labs(
    title = "BRCA1 FPKM Expression in Normal and Tumor Tissues",
    x = "FPKM",
    y = "Expression Quantity"
  )

ggsave(filename = "brca1_fpkm_expression.png", path = "../outputs", width = 10, height = 8)


# Heat Map
genes.of.interest <- c('BRCA1', 'BRCA2', 'TP53', 'ALK', 'MYCN')
# used the same gene list as the tutorial since, at the present moment, I
# haven't developed the skills or expertise in breast cancer bioinformatics
# needed to extract this information sans a tutorial

dat.long %>%
  filter(gene %in% genes.of.interest) %>%
  ggplot(., aes(x = sample, y = gene, fill = fpkm)) +
  geom_tile() +
  scale_fill_gradient(low = 'white', high = 'red') + 
  labs(
    title = "Sample FPKM Gene Expression",
    subtitle = "fpkm expression for select genes relevant to cancer",
    x = "Sample",
    y = "Genes"
  )

ggsave(filename = "gene_heat.png", path = "../outputs", width = 10, height = 8)

# Hana Rossie
