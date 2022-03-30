#' Plot Ioannina stratigraphic diagram
#' Resulting diagram requires external tweaking e.g. with Illustrator.

library(tidypaleo)
library(cowplot)

source("scripts/ioannina_read_data.R") # calls readxl, tidyverse
source("scripts/ioannina_zonation.R") # calls vegan, rioja, ggdendro
source("scripts/ioannina_ordination.R") # calls vegan




# Set theme ---------------------------------------------------------------

theme_set(theme_bw(8))
theme_update(legend.position = "none",
             panel.grid = element_blank(),
             panel.border = theme_border(c("left", "bottom"), size = 0.3),
             panel.background = element_rect(fill = "transparent"),
             plot.background = element_rect(fill = "transparent"))




# Plot taxa ---------------------------------------------------------------

# Clean count data.
cleaned_counts <- imported$counts %>%
  filter(
    rowSums(select(., !matches("depth"))) > 0  # remove samples with no diatoms
    & depth != 175.62  # remove outlier with poor preservation
  ) %>%
  column_to_rownames("depth")


# Transform counts to percentages.
cleaned_pc_ab <- decostand(cleaned_counts,
                method = "total",
                na.rm = "TRUE") %>%
  "*"(100)


# Isolate abundant taxa.
ab <- 4  # set required % abundance (keeps taxa present at >= this value)
ab_taxa_pc <- cleaned_pc_ab %>%
  select_if(~any(. >= ab)) %>%
  select(-contains("spp"))


# Mutate to tide (long) format.
ab_taxa_pc_tidy <- ab_taxa_pc %>%
  rownames_to_column("depth") %>%
  pivot_longer(-1, 
               names_to = "taxon",
               values_to = "pc_ab") %>%
  mutate(depth = as.numeric(depth),
         taxon = as_factor(taxon))


# Create plot.
plot_taxa <- ggplot(ab_taxa_pc_tidy, aes(x = pc_ab, y = depth)) +
  geom_col_segsh(size = 0.3) +
  scale_y_reverse(expand = c(0, 0),
                  breaks = seq(135, 285, 5)) +
  facet_abundanceh(vars(fct_inorder(ab_taxa_pc_tidy$taxon)),
                   dont_italicize = c("var.", "(≤20 µm)", "(>20 µm)",
                                      "(unidentifiable size)")) +
  labs(x = "%", y = "Depth (cm)")




# Plot life mode summary --------------------------------------------------

# Note: I have gaps in my record so can't use geom_area. Must use geom_ribbon 
# instead. This is a bit more complicated.

# Create function that pulls the names of taxa with different life modes.
pull_taxa <- function(taxon_type) {
  if (taxon_type == "planktonic") {
    filter(imported$taxa_life_modes, type == "planktonic") %>%
      select(taxon) %>%
      pull()
  } else if (taxon_type =="facultative") {
    filter(imported$taxa_life_modes, type == "facultative planktonic") %>%
      select(taxon) %>%
      pull()
  } else if (taxon_type == "benthic") {
    filter(imported$taxa_life_modes, type == "benthic") %>%
      select(taxon) %>%
      pull()
  }
}


# Create tidy dataframe containing abundances grouped by life mode.
life_mode <- imported$counts %>%
  # Sum counts by life mode.
  transmute(
    planktonic = rowSums(select(., contains(pull_taxa("planktonic")))),
    fac_plank = rowSums(select(., contains(pull_taxa("facultative")))),
    benthic = rowSums(select(., contains(pull_taxa("benthic"))))
  ) %>%
  # Convert life mode counts to percentage abundances
  decostand(., method = "total", na.rm = "TRUE") %>%
  "*"(100) %>%
  # Put in tidy (long) format.
  rownames_to_column(., "depth") %>%
  na_if(0) %>%
  pivot_longer(-1, 
               names_to = "life_mode",
               values_to = "pc_ab") %>%
  mutate(depth = as.numeric(depth),
         taxon = as_factor(life_mode))
  
  
# Set min and max values for geom_ribbon.
life_mode_min_max <- life_mode %>%
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


# Create plot.
plot_life_mode <- ggplot(life_mode_min_max,
                         aes(x = depth,
                             y = pc_ab,
                             fill = life_mode)) +
  geom_ribbon(aes(ymin = ymin,
                  ymax = ymax)) +
  coord_flip(clip = "off") +
  scale_x_reverse(expand = c(0, 0)) +
  scale_y_continuous(expand = c(0, 0)) +
  scale_fill_grey(start = 0.9,
                  end = 0.2) +
  labs(x = element_blank(),
       y = "%") +
  theme(axis.ticks.y = element_blank(),
        axis.text.y = element_blank())




# Plot concentration ------------------------------------------------------

# Set concentration of microspheres per ml used to create slides.
# NOTE: working concentration is usually x10^6

micros_conc <- 6810000

# Equation for calculating concentration:

# diatom concentration per 0.1g of sediment =
# total micros introduced * valves counted / micros counted

# NOTE: must multiply by 10 to get concentration per gram (the convention),
# rather than per 0.1g. This ridiculously simple thing alluded me for a stupidly
# long time!

# Calculate concentration and populate results in dataframe.
concentration <- imported$mass %>% 
  mutate(
    micros_introduced = micros_conc / (1 / mass),
    valves_counted = rowSums(select(imported$counts, -depth)),
    micros_counted = imported$microspheres$micros,
    diat_conc = micros_introduced * valves_counted / micros_counted * 10,
    diat_conc_6 = diat_conc / 1000000
  )


# Create plot.
plot_conc <- ggplot(concentration, 
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




# Plot F index ------------------------------------------------------------

# Calculate F index and put in dataframe.
f_index <- imported$f_index %>%
  mutate(f_index = whole / rowSums(select(., -depth))) %>%
  mutate(f_index = ifelse(whole == 0, 0, f_index))  # Convert NAs to 0s.


# Create plot.
plot_f_index <- ggplot(f_index, aes(x = f_index, y = depth)) +
  geom_lineh(size = 0.3) +
  scale_x_continuous(expand = c(0, 0),
                     limits = c(0, 1),
                     breaks = seq(0, 1, 0.5)) +
  scale_y_reverse(expand = c(0, 0)) +
  labs(x = element_blank(),
       y = element_blank()) +
  theme(axis.ticks.y = element_blank(),
        axis.text.y = element_blank())




# Plot P. ocellata preservation -------------------------------------------

# Create tidy dataframe of P. ocellata preservation data.
p_oc_pres <- imported$p_oc_pres %>%
  mutate(
    depth = depth,
    decostand(select(., -depth), 
              method = "total",
              na.rm = "TRUE") * 100
  ) %>%
  pivot_longer(-1, 
               names_to = "pres_quality",
               values_to = "pc_ab") %>%
  mutate(pres_quality = as_factor(pres_quality))


# Set min and max values for geom_ribbon.
p_oc_pres_min_max <- p_oc_pres %>%
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
    ymin = ifelse(pres_quality == "pristine" 
                  & (pc_ab + lead(pc_ab) + lead(pc_ab, 2)) == 0,
                  NA, ymin),
    ymax = ifelse(pres_quality == "pristine" 
                  & (pc_ab + lead(pc_ab) + lead(pc_ab, 2)) == 0,
                  NA, ymax)
  ) %>%
  mutate(
    ymin = ifelse(pres_quality == "dissolving" 
                  & (pc_ab + lag(pc_ab) + lead(pc_ab)) == 0,
                  NA, ymin),
    ymax = ifelse(pres_quality == "dissolving" 
                  & (pc_ab + lag(pc_ab) + lead(pc_ab)) == 0,
                  NA, ymax)
  ) %>%
  mutate(
    ymin = ifelse(pres_quality == "centre_only" 
                  & (pc_ab + lag(pc_ab) + lead(pc_ab)) == 0,
                  NA, ymin),
    ymax = ifelse(pres_quality == "centre_only" 
                  & (pc_ab + lag(pc_ab) + lag(pc_ab, 2)) == 0,
                  NA, ymax)
  )


# Create plot.
plot_p_oc_pres <- ggplot(p_oc_pres_min_max,
                         aes(x = depth,
                             y = pc_ab,
                             fill = pres_quality)) +
  geom_ribbon(aes(ymin = ymin,
                  ymax = ymax)) +
  coord_flip() +
  scale_x_reverse(expand = c(0, 0)) +
  scale_y_continuous(expand = c(0, 0),
                     breaks = seq(0, 100, 50)) +
  scale_fill_grey(start = 0.2,
                  end = 0.9) +
  labs(x = element_blank(),
       y = "%") +
  theme(axis.ticks.y = element_blank(),
        axis.text.y = element_blank())




# Plot PC1 ----------------------------------------------------------------

# Extract PC1 axis scores
pc1_samples = data.frame(scores(pca,
                                choices = 1,
                                display = "sites",
                                scaling = 0)) %>%
  mutate(depth = as.numeric(row.names(.)))


# Create plot.
plot_pc1 <- ggplot(pc1_samples, 
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




# Plot CONISS dendrogram --------------------------------------------------

plot_coniss <- ggplot(segment(ddata)) +
  geom_segment(size = 0.3, aes(x = x, y = y, xend = xend, yend = yend)) +
  coord_flip() +
  scale_y_continuous(expand = c(0, 0)) +
  scale_x_reverse(expand = c(0, 0),
                  breaks = NULL,
                  labels = NULL) +
  labs(x = "",
       y = "Total sum of squares") +
  theme(panel.border = theme_border("bottom", size = 0.5))




# Combine plots -----------------------------------------------------------

## Quicky plot (does not allow annotations)
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


cow_plot <- plot_grid(plotlist = grobs, 
                      nrow = 1, 
                      rel_widths = c(12, 2, 1, 1, 1, 1, 2))


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



