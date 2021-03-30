# 05 Stratigraphic plot
# —————————————————————

# 1. Load requirements ----------------------------------------------------

library(ggplot2)
library(tidypaleo)
library(cowplot)

source("scripts/02-manipulate.R")
source("scripts/03-zonation.R")
source("scripts/04-ordination.R")
source("scripts/borders-for-ggplot2.R")




# 2. Set theme ------------------------------------------------------------

theme_set(theme_bw(8))
theme_update(legend.position = "none",
             panel.grid = element_blank(),
             panel.border = theme_border(c("left", "bottom"), size = 0.3),
             panel.background = element_rect(fill = "transparent"),
             plot.background = element_rect(fill = "transparent"))




# 3. Plot taxa ------------------------------------------------------------

plot_taxa <- ggplot(taxa$rel_ab_4_tidy, aes(x = pc_ab, y = depth)) +
  geom_col_segsh(size = 0.3) +
  scale_y_reverse(expand = c(0, 0),
                  breaks = seq(135, 285, 5)) +
  facet_abundanceh(vars(fct_inorder(taxa$rel_ab_4_tidy$taxon)),
                   dont_italicize = c("var.", "(≤20 µm)", "(>20 µm)",
                                      "(unidentifiable size)")) +
  labs(x = "%", y = "Depth (cm)")




# 4. Plot life mode -------------------------------------------------------

# As I have gaps in my record, I can't use geom_area and must use geom_ribbon 
# instead. This is a big more complicated.

# Set min and max values for geom_ribbon.
taxa$life_mode_rel_ab_tidy <- taxa$life_mode_rel_ab_tidy %>%
  mutate(
    pc_ab = replace_na(pc_ab, 0),
    ymin = case_when(
      life_mode == "planktonic" ~ 0,
      life_mode == "fac_plank" ~ lag(pc_ab),
      life_mode == "benthic" ~ lag(pc_ab, 2) + lag(pc_ab)),
    ymax = case_when(
      life_mode == "planktonic" ~ pc_ab,
      life_mode == "fac_plank" ~ lag(pc_ab) + pc_ab,
      life_mode == "benthic" ~ lag(pc_ab, 2) + lag(pc_ab) + pc_ab
    )
  ) %>%
  mutate(
    ymin = ifelse(rowSums(select(., c("ymin", "ymax"))) == 0, NA, ymin),
    ymax = ifelse(rowSums(select(., c("ymin", "ymax"))) == 0, NA, ymax)
  )

plot_life_mode <- ggplot(taxa$life_mode_rel_ab_tidy,
                         aes(x = depth,
                             y = pc_ab,
                             fill = life_mode)) +
  geom_ribbon(aes(ymin = ymin,
                  ymax = ymax)) +
  coord_flip(clip = "off") +
  scale_x_reverse(expand = c(0, 0)) +
  scale_y_continuous(expand = c(0, 0)) +
  scale_fill_grey(start = 0,
                  end = 0.8) +
  labs(x = element_blank(),
       y = "%") +
  theme(axis.ticks.y = element_blank(),
        axis.text.y = element_blank())




# 5. Plot concentration ---------------------------------------------------

plot_conc <- ggplot(concentration$calculated, 
                    aes(x = diat_conc_6, y = depth)) +
  geom_col_segsh(size = 0.3) +
  geom_lineh(size = 0.3) +
  scale_x_continuous(expand = c(0, 0),
                     limits = c(0, 50),
                     breaks = seq(0, 50, 25)) +
  scale_y_reverse(expand = c(0, 0)) +
  labs(x = expression("valves (10"^6*" g"^"-1"*")"), # - converts to − in pdf
       y = element_blank()) +
  theme(axis.ticks.y = element_blank(),
        axis.text.y = element_blank()) 




# 6. Plot preservation ----------------------------------------------------

# F-index
plot_f_index <- ggplot(preservation$f_index, aes(x = f_index, y = depth)) +
  geom_lineh(size = 0.3) +
  scale_x_continuous(expand = c(0, 0),
                     limits = c(0, 1),
                     breaks = seq(0, 1, 0.5)) +
  scale_y_reverse(expand = c(0, 0)) +
  labs(x = element_blank(),
       y = element_blank()) +
  theme(axis.ticks.y = element_blank(),
        axis.text.y = element_blank())


# P. ocellata preservation

# Set min and max values for geom_ribbon.
preservation$p_ocellata_tidy <- preservation$p_ocellata_tidy %>%
  mutate(
    pc_ab = replace_na(pc_ab, 0),
    ymin = case_when(
      pres_quality == "pristine" ~ 0,
      pres_quality == "dissolving" ~ lag(pc_ab),
      pres_quality == "centre_only" ~ lag(pc_ab, 2) + lag(pc_ab)),
    ymax = case_when(
      pres_quality == "pristine" ~ pc_ab,
      pres_quality == "dissolving" ~ lag(pc_ab) + pc_ab,
      pres_quality == "centre_only" ~ lag(pc_ab, 2) + lag(pc_ab) + pc_ab
    )
  ) %>%
  # Plot samples correctly where there are no pristine valves but where there 
  # are dissolving and centre only valves
  mutate(
    ymin = ifelse(pres_quality == "pristine" & pc_ab == 0 & lead(pc_ab) != 0,
                  replace_na(pc_ab, 0), ymin),
    ymax = ifelse(pres_quality == "pristine" & pc_ab == 0 & lead(pc_ab) != 0,
                  replace_na(pc_ab, 0), ymax)
  ) %>%
  # Create gaps with horizonal edges where samples are missing or no diatoms 
  # present
  mutate(
    ymin = ifelse(pres_quality == "pristine" & (pc_ab + lead(pc_ab) + lead(pc_ab, 2)) == 0,
                  NA, ymin),
    ymax = ifelse(pres_quality == "pristine" & (pc_ab + lead(pc_ab) + lead(pc_ab, 2)) == 0,
                  NA, ymax)
  ) %>%
  mutate(
    ymin = ifelse(pres_quality == "dissolving" & (pc_ab + lag(pc_ab) + lead(pc_ab)) == 0,
                  NA, ymin),
    ymax = ifelse(pres_quality == "dissolving" & (pc_ab + lag(pc_ab) + lead(pc_ab)) == 0,
                  NA, ymax)
  ) %>%
  mutate(
    ymin = ifelse(pres_quality == "centre_only" & (pc_ab + lag(pc_ab) + lead(pc_ab)) == 0,
                  NA, ymin),
    ymax = ifelse(pres_quality == "centre_only" & (pc_ab + lag(pc_ab) + lag(pc_ab, 2)) == 0,
                  NA, ymax)
  )

plot_p_oc_pres <- ggplot(preservation$p_ocellata_tidy,
                         aes(x = depth,
                             y = pc_ab,
                             fill = pres_quality)) +
  geom_ribbon(aes(ymin = ymin,
                  ymax = ymax)) +
  coord_flip() +
  scale_x_reverse(expand = c(0, 0)) +
  scale_y_continuous(expand = c(0, 0),
                     breaks = seq(0, 100, 50)) +
  scale_fill_grey(start = 0,
                  end = 0.8) +
  labs(x = element_blank(),
       y = "%") +
  theme(axis.ticks.y = element_blank(),
        axis.text.y = element_blank())




# 7. Plot PC1 -------------------------------------------------------------

plot_pc1 <- ggplot(ordination$pca_scrs, 
                   aes(x = PC1, y = depth)) +
  geom_lineh(size = 0.3) +
  scale_x_continuous(expand = c(0, 0),
                     limits = c(-0.2, 0.22),
                     breaks = seq(-0.2, 0.2, 0.2)) +
  scale_y_reverse(expand = c(0, 0)) +
  labs(x = element_blank(),
       y = element_blank()) +
  theme(axis.ticks.y = element_blank(),
        axis.text.y = element_blank())




# 8. Plot CONISS ----------------------------------------------------------

plot_coniss <- ggplot(segment(coniss$ddata)) +
  geom_segment(size = 0.3, aes(x = x, y = y, xend = xend, yend = yend)) +
  coord_flip() +
  scale_y_continuous(expand = c(0, 0)) +
  scale_x_reverse(expand = c(0, 0),
                  breaks = NULL,
                  labels = NULL) +
  labs(x = "",
       y = "Total sum of squares") +
  theme(panel.border = theme_border("bottom", size = 0.5))




# 9. Combine --------------------------------------------------------------

# Quicky plot (does not allow annotations)
# library(patchwork)
# plot_taxa + plot_life_mode + plot_conc + plot_f_index + plot_p_oc_pres + 
#   plot_pc1 + plot_coniss + 
#   plot_layout(nrow = 1, widths = c(12, 2, 1, 1, 1, 1, 2)) 

# Align plots (add NULL plots to create space for labels above plot)
plots <- list(plot_taxa, plot_life_mode, plot_conc, plot_f_index, plot_p_oc_pres, plot_pc1, plot_coniss)
grobs <- lapply(plots, as_grob)
plot_heights <- lapply(grobs, function(x) {x$heights})
aligned_heights <- align_margin(plot_heights, "first")
aligned_heights <- align_margin(aligned_heights, "last")
for (i in seq_along(plots)) {
  grobs[[i]]$heights <- aligned_heights[[i]]
}

cow_plot <- plot_grid(plotlist = grobs, nrow = 1, rel_widths = c(12, 2, 1, 1, 1, 1, 2))

ggdraw(cow_plot) +
  # Manually add labels and lines - might be best to do this in Illustrator 
  # instead as labels will move depending on plot window size.
  draw_label("Planktonic", x = 0.6, y = 0.83, hjust = 0, vjust = 0, size = 7, angle = 45) +
  draw_label("Facultative Planktonic", x = 0.65, y = 0.83, hjust = 0, vjust = 0, size = 7, angle = 45) +
  draw_label("Benthic", x = 0.67, y = 0.83, hjust = 0, vjust = 0, size = 7, angle = 45) +
  draw_label("Concentration", x = 0.695, y = 0.83, hjust = 0, vjust = 0, size = 7, angle = 45) +
  draw_label("F-index", x = 0.745, y = 0.83, hjust = 0, vjust = 0, size = 7, angle = 45) +
  draw_label("Pristine", x = 0.795, y = 0.83, hjust = 0, vjust = 0, size = 7, angle = 45) +
  draw_label("Dissolving", x = 0.81, y = 0.83, hjust = 0, vjust = 0, size = 7, angle = 45) +
  draw_label("Centre only", x = 0.825, y = 0.83, hjust = 0, vjust = 0, size = 7, angle = 45) +
  draw_label("PC1", x = 0.845, y = 0.83, hjust = 0, vjust = 0, size = 7, angle = 45) +
  draw_label("CONISS", x = 0.915, y = 0.83, hjust = 0, vjust = 0, size = 7, angle = 45)


