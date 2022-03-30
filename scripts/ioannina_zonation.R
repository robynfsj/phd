#' Perform Ioannina zonation (via constrained cluster analysis)

library(vegan)  # decostand, designdist
library(rioja)  # chclust
library(ggdendro)  # dendro_data

source("scripts/ioannina_read_data.R") # also runs libraries readxl, tidyverse
source("scripts/borders-for-ggplot2.R")




# Clean data --------------------------------------------------------------

# Clean raw count data.
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

# hellinger = sqrt of relative proportion
coniss_data <-decostand(ab_taxa_counts, 
                        method = "hellinger", 
                        na.rm = "TRUE")




# Cluster analysis --------------------------------------------------------

# Compute distance matrix.
dist_matrix <- designdist(coniss_data,
                          method = "A+B-2*J",
                          terms = "quadratic")


# Run cluster analysis.
coniss_results <- chclust(dist_matrix, method = "coniss")


# Print broken stick model to determine number of zones.
bstick(coniss_results, ng = 30)  
# 12 significant zones




# Plot dendrogram ---------------------------------------------------------

# Extract data for plotting dendrogram.
ddata<- dendro_data(coniss_results, type = "rectangle")


# Modify x values so leaves are plotted by depth in core not sequential order.
new_x <- approxfun(ddata$labels$x, 
                   as.numeric(as.character(ddata$labels$label))) 
ddata$segments$x <- new_x(ddata$segments$x)
ddata$segments$xend <- new_x(ddata$segments$xend)


# Plot dendrogram.
ggplot(segment(ddata)) +
  geom_vline(xintercept = c(156.02, 163.22, 177.62, 190.42, 199.22, 215.22, 
                            233.62, 253.62, 258.42, 275.22, 279.22),
             colour = "lightgrey") +
  geom_segment(aes(x = x, y = y, xend = xend, yend = yend)) +
  coord_flip() +
  scale_y_continuous() +
  scale_x_reverse(breaks = seq(135, 285, 5)) +
  labs(x = "Depth (m)",
       y = "Total sum of squares") +
  theme_bw(10) +
  theme(aspect.ratio = 3,
        legend.position = "none",
        panel.grid = element_blank(),
        panel.border = theme_border(c("left", "bottom")),
        panel.background = element_rect(fill = "transparent"),
        plot.background = element_rect(fill = "transparent")) +
  geom_hline(yintercept = 14,
             linetype = 2)



