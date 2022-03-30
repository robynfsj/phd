#' Perform Ioannina ordination

library(vegan)  # decostand

source("scripts/ioannina_read_data.R") # also runs libraries readxl, tidyverse




# Clean data --------------------------------------------------------------

# Clean count data.
cleaned_counts <- imported$counts %>%
  filter(
    rowSums(select(., !matches("depth"))) > 0  # remove samples with no diatoms
    & depth != 175.62  # remove outlier with poor preservation
  ) %>%
  column_to_rownames("depth")


# Get names of abundant taxa.
ab <- 0.04  # set required abundance (keeps taxa present at >= this value)
abundant_taxa <- cleaned_counts %>%
  decostand(method = "total", na.rm = TRUE) %>%
  select_if(~any(. >= ab)) %>%
  select(-contains("spp")) %>%
  colnames()


# Isolate data for only abundant taxa.
ab_taxa_counts <- cleaned_counts %>%
  select(all_of(abundant_taxa))




# Transform data ----------------------------------------------------------

ord_data <- decostand(ab_taxa_counts, 
                      method = "total", 
                      na.rm = "TRUE") %>%
  "*"(100) %>%  # convert rel. ab. out of 1.0 to percentage
  sqrt()  # sqrt the percentage

# Note: could have just used method = "hellinger" (sqrt of relative proportions 
# out of 1.0) for same effective result (but eigenvalues would be to a different 
# power of 10). However, have chosen to calculate percentage relative abundance 
# first and then square root to keep result consistent with that obtained when I 
# first ran analyses with the Tilia program.




# Perform ordination ------------------------------------------------------

# Run DCA.
dca <- decorana(ord_data)
dca


# Run PCA.
pca <- rda(ord_data)
pca


# Find significant axes
screeplot(pca, bstick = TRUE, type = "l", main = NULL)
# 3 sig. axes


# Proportion of variance explained by axes.
sum((as.vector(pca$CA$eig) / sum(pca$CA$eig))[1]) * 100  # 25.89 %
sum((as.vector(pca$CA$eig) / sum(pca$CA$eig))[2]) * 100  # 15.19 %
sum((as.vector(pca$CA$eig) / sum(pca$CA$eig))[3]) * 100  # 11.75 %
sum((as.vector(pca$CA$eig) / sum(pca$CA$eig))[4]) * 100  # 7.45 %


# Cumulative variance explained by axes
sum((as.vector(pca$CA$eig) / sum(pca$CA$eig))[1]) * 100  # 25.89 %
sum((as.vector(pca$CA$eig) / sum(pca$CA$eig))[1:2]) * 100  # 41.08 %
sum((as.vector(pca$CA$eig) / sum(pca$CA$eig))[1:3]) * 100  # 52.83 %
sum((as.vector(pca$CA$eig) / sum(pca$CA$eig))[1:4]) * 100  # 60.27 %




# Biplots -----------------------------------------------------------------

# Quick biplot of PC1 and PC2.
biplot(pca, type = "points", choices=c(1,2))



