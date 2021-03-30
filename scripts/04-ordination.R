# 04 Ordination
# —————————————
# I want to perform an ordination on my diatom results to pick out the main 
# variation in the diatom assemblage throughout my record and to see if I can 
# find out what might be driving this variation.


# 1. Load requirements ----------------------------------------------------

source("scripts/02-manipulate.R")




# 2. Create ordination list for new dfs -----------------------------------

ordination <- list()




# 3. Prepare data for ordination ------------------------------------------

# The dfs have sample_no as the first column and still contain samples within 
# which no diatoms were counted. These must be removed for ordination but will
# remain in the original dfs for any potential future analyses. Although not 
# ideal, sample number is moved to row names so that this information is 
# preserved after ordination.

# Create function that removes samples with no diatoms (and outlier with diatoms 
# but low count) and moves sample number to row names so data is ready for 
# ordination.
ord_data <- function(df) {
  df %>%
    filter(sample_no != "189" & rowSums(.[, -1]) > 0) %>%
    column_to_rownames("sample_no")
}




# 4. DCA ------------------------------------------------------------------

# Ordination is run on percentage relative abundances (i.e. transformed count
# data) and only on abundant taxa present at ≥4 % in at least one sample 
# (dataframe: rel_ab_4).

ordination$dca <- decorana(ord_data(taxa$rel_ab_4))
ordination$dca
# DCA1 = 1.9914

# Gradient length of DCA1 is <2.5 so linear ordination methods are appropriate
# (Legendre & Birks, 2012).




# 5. PCA ------------------------------------------------------------------

ordination$pca <- rda(ord_data(taxa$rel_ab_4))

# Proportion of variance explained by PC1 and PC2.
sum((as.vector(ordination$pca$CA$eig) / sum(ordination$pca$CA$eig))[1:2]) * 100
# 80 %




# 6. Extract PC1 scores ---------------------------------------------------

pc1 = data.frame(scores(ordination$pca, 
                        choices = 1, 
                        display = "sites", 
                        scaling = 0)) %>%
  mutate(sample_no = as.numeric(row.names(.)))

ordination$pca_scrs <- left_join(imported$depths, pc1)




# 7. Plot PCA diagram -----------------------------------------------------

# library(ggplot2)
# 
# 
# # 7.1 PC1 ----
# # ————————————
# plot_pc1 <- ggplot(ordination$pca_scrs, aes(x = depth, y = PC1)) +
#   geom_line() +
#   coord_flip() +
#   scale_x_reverse(expand = c (0, 0), breaks = seq(130, 290, 10)) +
#   scale_y_continuous(expand = c(0, 0)) +
#   labs(x = "Depth (m)") +
#   theme_classic()
# 
# plot_pc1




# # 7.2 PCA (built from scratch using base R) ----
# # ——————————————————————————————————————————————
# # NOTE: This is specific to the pca results. It needs to be manually edited
# # depending on those results. You can plot any results quickly but the pca plot 
# # will be crowded and difficult to read. This is the only way to make it look
# # nice using R.
# 
# # Quick plots to see the data
# par(mfrow = c(2, 2))
# plot(ordination$pca, scaling = 0, main = "raw scores")
# plot(ordination$pca, scaling = 1, main = "sites scores scaled")
# plot(ordination$pca, scaling = 2, main = "species scores scaled")
# plot(ordination$pca, scaling = 3, main = "both scores scaled")
# 
# biplot(ordination$pca, display = "species", scaling = 0, main = "raw scores")
# biplot(ordination$pca, display = "species", scaling = 1, main = "sites scores scaled")
# biplot(ordination$pca, display = "species", scaling = 2, main = "species scores scaled")
# biplot(ordination$pca, display = "species", scaling = 3, main = "both scores scaled")
# 
# 
# 
# 
# 
# 
# scrs <- scores(ordination$pca, scaling = 0)
# 
# scrs_species <- as.data.frame(scrs[["species"]])
# scrs_sites <- as.data.frame(scrs[["sites"]])
# 
# 
# # The quick plot shows that most species plot very close together. Select only
# # the species that drive the most variation in the samples to plot (those with
# # the highest/lowest PC1 and PC2 scores).
# 
# scrs_species <- scrs_species %>%
#   mutate(
#     taxon = row.names(.),
#     PC1_squared = PC1 ^ 2,
#     PC2_squared = PC2 ^ 2
#   ) %>%
#   format(scientific = FALSE) 
# 
# pc1_taxa <- scrs_species %>%
#   slice_max(PC1_squared, n = 8)
# 
# pc2_taxa <- scrs_species %>%
#   slice_max(PC2_squared, n = 8)
# 
# scrs_species <- full_join(pc1_taxa, pc2_taxa)
# 
# labels <- c(expression(italic("C. oc")), 
#             expression(paste(italic("C"),". sp 1")),
#             expression(italic("S. ve")),
#             expression(italic("S. pi")),
#             expression(italic("E. mi")),
#             expression(italic("F. lu")))
# 
# # manually offset species scores for text so that they don't plot overe arrow
# species.scrs.off <- species.scrs * 1.07
# 
# # plot blank axes
# xlim <- range(scrs[["species"]][,1], scrs[["sites"]][,1])
# ylim <- range(scrs[["species"]][,2], scrs[["sites"]][,2])
# 
# plot.new()
# plot.window(xlim = xlim, ylim = ylim, asp = 1)
# abline(h = 0, lty = "dotted")
# abline(v = 0, lty = "dotted")
# axis(side = 1)
# axis(side = 2)
# title(xlab = "PC 1", ylab = "PC 2")
# box()
# 
# # create MIS vector
# mis <- c(rep("MIS 7", 33),
#          rep("MIS 8", 30),
#          rep("MIS 9", 18))
# 
# # plot samples by mis
# points(sites.scrs,
#        cex=0.7)
# 
# 
# # plot arrows of selected species       
# arrows(0, 0, species.scrs$PC1, species.scrs$PC2,
#        length=0.1)
# 
# # add selected species
# text(species.scrs$PC2~species.scrs$PC1,
#      labels = labels,
#      cex = 0.8,
#      col = "black",
#      pos = c(4, 4, 1, 4, 2, 2)) 
# 
# # manually create legend
# legend("topright", 
#        legend = c("MIS 7", "MIS 8", "MIS 9"),
#        col=c("yellowgreen","lightskyblue3","coral"),
#        pch=c(15, 17, 19),
#        bty=1)



# 8. Clean up -------------------------------------------------------------

rm(pc1)




