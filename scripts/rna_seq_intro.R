# RNA-Seq Introduction Notebooks 
# Authors: Acosta, Joe
# Date: 07/11/2026

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


# ---------------------------------------

# Learning Bar plots, density plots, box plots
# scatterplots, and heatmap in R for Bioinformatics
# using ggplot

library(ggplot2)

# base R allows for the generation of these visualizations
# though it is much harder to modify

# Basic format for ggplot function
# ggplot(data, asthetics, )
# data is the data that is going to be used
# asthetics is what I want as my x and y axis

# ggplot(dat.long, aes(x = variable, y = variable))
# like with tidyverse, where we use the pipe opperator to add functionality, 
# with ggplot we use the + opperator to give it functionality

# right now all we are saying is that we want to create a plot that's going to
# use the data dataframe and the aesthetics x and y variables
# The function has no idea what type of plot we want to create so we have to
# give it additional information on what kind of plot we want to give it

# basic format or syntax
ggplot(dat.long, aes(x = , y = )) +
  geom_col()

# can add more information with the plus operator

# Barplot
# looking at the BRCA1 gene
dat.long %>%
  filter(gene == "BRCA1") %>%
  ggplot(., aes(x = samples, y = FPKM)) +
  geom_col()

# don't know which is the tumor or normal; use fill to see which are which
dat.long %>%
  filter(gene == "BRCA1") %>%
  ggplot(., aes(x = samples, y = FPKM, fill = tissue)) +
  geom_col()

# compare the expression of samples for BRCA1

# Density Plot:
# How does the distribution of... compare across the tumor and non tumor
dat.long %>%
  filter(gene == "BRCA1") %>%
  ggplot(., aes(x = FPKM, fill = tissue)) + # fill at tissue since want to
  # compare how distribution of expression compares between tumor and normal
  geom_density() # the normal tissue is obfuscating the brest tumor tissue.

dat.long %>%
  filter(gene == "BRCA1") %>%
  ggplot(., aes(x = FPKM, fill = tissue)) +
  geom_density(alpha = .6)

# Boxplot

# compare between samples having different metastasis status

dat.long %>%
  filter(gene == "BRCA1") %>%
  ggplot(., aes(x = metastasis, y = FPKM)) +
  geom_boxplot()

# convert into a violin plot
dat.long %>%
  filter(gene == "BRCA1") %>%
  ggplot(., aes(x = metastasis, y = FPKM)) +
  geom_violin()


# Scatter Plot
# Compare the expression between two genes

dat.long %>%
  filter(gene == "BRCA1" | gene == "BRCA2") %>%
  spread(key = gene, value = FPKM) %>%
  ggplot(., aes(x = BRCA1, y = BRCA2)) +
  geom_point()

# how to fit a line of best fit to this plot
# can't say if it is significant without a correlation test

dat.long %>%
  filter(gene == "BRCA1" | gene == "BRCA2") %>%
  spread(key = gene, value = FPKM) %>%
  ggplot(., aes(x = BRCA1, y = BRCA2)) +
  geom_point() +
  geom_smooth(method = 'lm', se = FALSE) # lm is straight line
  # se is no confidence intervals


# seperate out the two different tissue types
dat.long %>%
  filter(gene == "BRCA1" | gene == "BRCA2") %>%
  spread(key = gene, value = FPKM) %>%
  ggplot(., aes(x = BRCA1, y = BRCA2, colour = tissue)) + # why do we use the
  # color and not the fill??
  geom_point() +
  geom_smooth(method = 'lm', se = FALSE)

# Heatmaps

# easy to visualize expressions of multiple genes across the samples and compare
# the expressions of all these genese simultaniously

genes.of.interest <- c('BRCA1', 'BRCA2', 'TP53', 'ALK', 'MYCN')

# compare these five genes across all the samples

# start fetching the data

dat.long %>%
  filter(gene %in% genes.of.interest) %>%
  ggplot(., aes(x = samples, y = gene, fill = FPKM)) +
  geom_tile()

# darker the color is lower the expression
# lighter the color is the lighter the expression
# contrast @ with changing the colors of the scale

dat.long %>%
  filter(gene %in% genes.of.interest) %>%
  ggplot(., aes(x = samples, y = gene, fill = FPKM)) +
  geom_tile() +
  scale_fill_gradient(low = "white", high = 'red')

# heatmap allows us to make multiple comparisons and see trends across the data
# we can visualize trends across the data easily

# Saving the plots in multiple formats, publications, or records
# save it to a var or

p <- dat.long %>%
  filter(gene %in% genes.of.interest) %>%
  ggplot(., aes(x = samples, y = gene, fill = FPKM)) +
  geom_tile() +
  scale_fill_gradient(low = "white", high = 'red')

# or funct ggsave()

ggsave(p, filename = "heatmap_save1.pdf", width = 10, height =  8)
# can save as png too, can provide the width and height as in, cm, or other units

# another way to save plots
pdf("heatmap_save2.pdf", width = 10, height = 8)
dat.long %>%
  filter(gene %in% genes.of.interest) %>%
  ggplot(., aes(x = samples, y = gene, fill = FPKM)) +
  geom_tile() +
  scale_fill_gradient(low = "white", high = 'red')
dev.off()

# ---------------------------------------

# Understanding the difference between RPKM/FPKM and TPM

# when to use them and when they are not appropriate to be used

# Replicates: give us confidence that the observed effect is a result of biological
# phenomenon and not as a side effect of a technical issue

# Two types of replicates

# Technical Replicates: are several measurements taken from the same organism/sample

# Biological Replicates: are several measurements each taken from independent
# organisms/samples

# Counts: reads that map to a gene that gives us the expression of the gene in
# that experiment

# Counts Matrix has rows as genes and columns as sample replicates and values are
# the counts

# We must normalize the counts because it doesn't account for biasis

# Gene length
# Sequencing Depth

# RPKM: normalizes for gene length and sequencing depth
# higher the RPKM, higher the expression
# used to quantify transcripts from single-ended reads
# Cannot be used for Differential gene expression;
# RPKM doesn't account biasis in difference in biological conditions

# FPKM: analogous to RPKM
# main difference is that it is used for paired data rather than single ended
# Higher FPKM is higher the expression
# cant be used for DESEq2 either

# TPM: normalizes for gene length and sequencing
# is better suited to compare expression bethween samples
# cant be used for DESEq2 either








