# 02 Manipulate data 
# ——————————————————
# This script manipulates the raw data that was loaded in the first script, 
# creating several new dataframes that will be used later on. Numerous 
# dataframes are created and are placed within lists for clarity. Resulting 
# dataframes after running this script:
#
# p_ocellata
#   |
#   ├ all_splits — diatom counts of the taxon P. ocellata split into size and
#   |              ocelli number
#   |
#   ├ size_split — counts of P. ocellata split into size classes
#   |
#   ├ size_split_rel_ab — percentage relative abundance of size_split, size 
#   |                     class abundances relative to P. ocellata
#   |
#   ├ size_split_rel_ab_all — percentage relative abundance of size_split, size
#   |                         class abundances relative to all taxa
#   |
#   ├ ocelli_split — counts of P. ocellata split into ocelli number classes
#   |
#   ├ ocelli_split_rel_ab — percentage relative abundance of ocelli_split, 
#   |                       ocelli class abundances relative to P. ocellata
#   |
#   ├ ocelli_split_rel_ab_all — percentage relative abundance of ocelli_split, 
#   |                           ocelli class abundances relative to all taxa
#   |
#   ├ size_split_20µm — counts of P. ocellata split into two size classes of
#                       ≤20 µm and >20µm.
#
# taxa
#   |
#   ├ counts — counts of all diatom taxa with P. ocellata not split  
#   |
#   ├ rel_ab — relative abundances of all diatom taxa with P. ocellata not split  
#   |
#   ├ rel_ab_4 — rel_ab but only taxa present at ≥4 % in at least one sample
#   |
#   ├ rel_ab_4_tidy — rel_ab_4 in long form for ggplot
#   |
#   ├ counts_split_20µm — counts of all diatom taxa with P. ocellata split into
#   |                     ≤20 µm, >20 µm and unidentifiable size
#   |
#   ├ rel_ab_split_20µm — relative abundances of all diatom taxa with 
#   |                     P. ocellata split into≤20 µm, >20 µm and 
#   |                     unidentifiable size
#   |
#   ├ rel_ab_split_20µm_4 — rel_ab_split_20µm but only taxa present at ≥4 % in 
#   |                       at least one sample
#   |
#   ├ rel_ab_split_20µm_4_tidy — rel_ab_split_20µm_4 in long form for ggplot
#   |
#   ├ life_mode — counts of all diatoms grouped into three types (planktonic,
#   |             facultative planktonic and benthic)
#   |
#   ├ life_mode_rel_ab — relative abundances of all diatoms grouped into three
#   |                    types (planktonic, facultative planktonic and benthic)
#   |
#   ├ life_mode_rel_ab_tidy — life_mode_rel_ab in long form for ggplot
#
# concentration
#   |
#   ├ measurments — data used to calculate diatom concentrations
#   |
#   ├ calculated — calculated diatom concentrations for each slide
#
# preservation
#   |
#   ├ f_index — proportion of whole valves to all valves counted (out of 1)
#   |
#   ├ p_ocellata — percentage relative abundance of P. ocellata valves that are
#   |              pristine, dissolving or have only the centre remaining
#   |
#   ├ p_ocellata_tidy — p_ocellata in long form for ggplot



# 1. Load requirements ----------------------------------------------------

library(vegan)
library(dplyr)
library(tibble)
library(tidyr)
library(forcats)

source("scripts/01-load.R")




# 2. Create lists for dataframes ------------------------------------------

p_ocellata <- list()
taxa <- list()
concentration <- list()
preservation <- list()




# 3. P. ocellata splits ---------------------------------------------------

# Counts of the taxon Pantocsekiella ocellata were split into different groups
# according to their size and number of ocelli. These dataframes group the 
# P. ocellata counts according to their size and separately according to their
# ocelli number. Percentage relative abundances are then calculated.




# 3.1 Create df: all_splits ----
# ——————————————————————————————
p_ocellata$all_splits <- imported$p_oc_breakdown




# 3.2 Create size split dfs ----
# ——————————————————————————————

p_ocellata$size_split <- p_ocellata$all_splits %>%
  transmute(
    sample_no = sample_no,
    s_1_5 = rowSums(select(., ends_with("1–5"))),
    s_6_10 = rowSums(select(., ends_with("6–10"))),
    s_11_15 = rowSums(select(., ends_with("11–15"))),
    s_16_20 = rowSums(select(., ends_with("16–20"))),
    s_21_25 = rowSums(select(., ends_with("21–25"))),
    s_26_30 = rowSums(select(., ends_with("26–30"))),
    s_31_35 = rowSums(select(., ends_with("31–35"))),
    s_36_40 = rowSums(select(., ends_with("36–40"))),
    s_41_45 = rowSums(select(., ends_with("41–45"))),
    s_u = rowSums(select(., ends_with("U")))
  )

p_ocellata$size_split_rel_ab <- p_ocellata$size_split %>%
  mutate(sample_no = sample_no,
         decostand(select(., -sample_no), 
                   method = "total", 
                   na.rm = "TRUE") * 100)

p_ocellata$size_split_rel_ab_all <- p_ocellata$size_split %>%
  transmute(
    sample_no = sample_no,
    select(., -sample_no) / rowSums(select(imported$counts, -sample_no)) * 100
  )




# 3.3 Create ocelli split dfs ----
# ————————————————————————————————

p_ocellata$ocelli_split <- p_ocellata$all_splits %>%
  transmute(
    sample_no = sample_no,
    o_0 = rowSums(select(., starts_with("0"))),
    o_1 = rowSums(select(., starts_with("1"))),
    o_2 = rowSums(select(., starts_with("2"))),
    o_3 = rowSums(select(., starts_with("3"))),
    o_4 = rowSums(select(., starts_with("4"))),
    o_5 = rowSums(select(., starts_with("5"))),
    o_6 = rowSums(select(., starts_with("6"))),
    o_7 = rowSums(select(., starts_with("7"))),
    o_u = rowSums(select(., starts_with("U")))
  )

p_ocellata$ocelli_split_rel_ab <- p_ocellata$ocelli_split %>%
  mutate(
    sample_no = sample_no,
    decostand(select(., -sample_no), 
              method = "total",
              na.rm = "TRUE") * 100
  )

p_ocellata$ocelli_split_rel_ab_all <- p_ocellata$ocelli_split %>%
  transmute(
    sample_no = sample_no,
    select(., -sample_no) / rowSums(select(imported$counts, -sample_no)) * 100
  )




# 3.4 Create df: split_20µm ----
# ——————————————————————————————

p_ocellata$split_20µm <- p_ocellata$size_split %>%
  transmute(
    sample_no = sample_no,
    lt20 = rowSums(select(., c(s_1_5, s_6_10, s_11_15, s_16_20))),
    mt20 = rowSums(select(., c(s_21_25, s_26_30, s_31_35, s_36_40, s_41_45))),
    u = s_u
  )




# 4. Taxa -----------------------------------------------------------------

# The taxa dataframes contain the counts of the individual diatom taxa and 
# counts grouped according to life mode. Two dataframes are created for the 
# counts of individual diatom taxa. One with all P. ocellata grouped together
# and one where they are split into ≤20 µm and > 20 µm.




# 4.1 Create taxa abundance dfs ----
# ——————————————————————————————————
taxa$counts <- imported$counts

taxa$rel_ab <- taxa$counts %>%
  mutate(
    sample_no = sample_no,
    decostand(select(., -sample_no), 
              method = "total",
              na.rm = "TRUE") * 100
  )

# Create function to get names of abundant taxa.
abundant_taxa <- function(df, min_percentage_abundance) {
  df %>%
    select_if(~any(. >= min_percentage_abundance)) %>%
    select(-contains("spp")) %>%
    select(-sample_no) %>%
    colnames()
}

taxa$rel_ab_4 <- taxa$rel_ab %>%
  transmute(
    sample_no = sample_no,
    select(., all_of(abundant_taxa(., 4)))
  )

taxa$rel_ab_4_tidy <- inner_join(imported$depths, taxa$rel_ab_4) %>%
  pivot_longer(-(1:2), 
               names_to = "taxon",
               values_to = "pc_ab") %>%
  mutate(taxon = as_factor(taxon))




# 4.2 Create taxa abundance dfs with P. ocellata 20µm split ----
# ———————————————————————————————————————————————————————————————
# NOTE: Proper taxon names used for easy plotting later.

taxa$counts_split_20µm <- taxa$counts %>%
  mutate(
    "Pantocsekiella ocellata" = NULL,
    "Pantocsekiella ocellata (≤20 µm)" = p_ocellata$split_20µm$lt20,
    "Pantocsekiella ocellata (>20 µm)" = p_ocellata$split_20µm$mt20,
    "Pantocsekiella ocellata (unidentifiable size)" = p_ocellata$split_20µm$u
  ) %>%
  relocate(c("Pantocsekiella ocellata (≤20 µm)", 
             "Pantocsekiella ocellata (>20 µm)",
             "Pantocsekiella ocellata (unidentifiable size)"), 
           .after = "Pantocsekiella minuscula")

taxa$rel_ab_split_20µm <- taxa$counts_split_20µm %>%
  mutate(
    sample_no = sample_no,
    decostand(select(., -sample_no), 
              method = "total",
              na.rm = "TRUE") * 100
  )

taxa$rel_ab_split_20µm_4 <- taxa$rel_ab_split_20µm %>%
  transmute(
    sample_no = sample_no,
    select(., all_of(abundant_taxa(., 4)))
  )

taxa$rel_ab_split_20µm_4_tidy <- inner_join(imported$depths, 
                                            taxa$rel_ab_split_20µm_4) %>%
  pivot_longer(-(1:2), 
               names_to = "taxon",
               values_to = "pc_ab") %>%
  mutate(taxon = as_factor(taxon))




# 4.3 Create life mode dfs ----
# —————————————————————————————

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

taxa$life_mode <- imported$counts %>%
  transmute(
    sample_no = sample_no,
    planktonic = rowSums(select(., contains(pull_taxa("planktonic")))),
    fac_plank = rowSums(select(., contains(pull_taxa("facultative")))),
    benthic = rowSums(select(., contains(pull_taxa("benthic"))))
  )

taxa$life_mode_rel_ab <- taxa$life_mode %>%
  mutate(
    sample_no = sample_no,
    decostand(select(., -sample_no), 
              method = "total",
              na.rm = "TRUE") * 100
  )

taxa$life_mode_rel_ab_tidy <- inner_join(imported$depths,
                                         taxa$life_mode_rel_ab) %>%
  na_if(0) %>%
  pivot_longer(cols = -(1:2), 
               names_to = "life_mode",
               values_to = "pc_ab") %>%
  mutate(life_mode = as_factor(life_mode))




# 5. Concentration --------------------------------------------------------

# The concentration dataframes contain data on the diatom concentration in 
# each sample. Zero values need to be left in as they indicate preservation was 
# so bad that no diatoms were present.


# 5.1 Set microsphere concentration ----
# ——————————————————————————————————————
# Concentration of microspheres used to create slides.

micros_conc = 68100000 




# 5.2 Create df: measurements ----
# ————————————————————————————————

concentration$measurements <- imported$mass %>% 
  mutate(
    micros_introduced = micros_conc / (1 / mass),
    micros_counted = imported$microspheres$micros
  )




# 5.3 Create df: concentration ----
# —————————————————————————————————

concentration$calculated <- transmute(
  concentration$measurements,
  sample_no = sample_no,
  depth = imported$depths$depth,
  diat_conc = micros_introduced * rowSums(select(imported$counts, -sample_no)) / 
    micros_counted,
  diat_conc_6 = diat_conc / 1000000
)

# Remove temp. objects.
rm(micros_conc)




# 6. Preservation ---------------------------------------------------------

# Two ways of measuring preservation are investigated. The F index is the 
# total number of whole valves as a proportion of all valves counted. The
# preservation of only the P. ocellata valves is measured by the relative
# proportions of their valves that are pristine, dissolving or have only 
# the central area remaining.




# 6.1 Create df: f_index ----
# ———————————————————————————
# Missing samples need to be skipped so they remain as NA. Samples containing no 
# diatoms should have an F index of 0 as they contain evidence that samples 
# would have contained diatoms (e.g. tiny fragments of diatoms are present) but 
# that they have not been preserved very well. The code allows for this by 
# setting the F index of samples containing no whole diatom valves to 0.

preservation$f_index <- imported$f_index %>%
  mutate(f_index = whole / rowSums(select(., -sample_no))) %>%
  mutate(f_index = ifelse(whole == 0, 0, f_index)) %>%
  left_join(imported$depths)




# 6.2 Create P. ocellata preservation dfs ----
# ————————————————————————————————————————————
# Replaces the raw count for each P. ocellata category with percentage 
# abundance. Samples containing no diatoms need to be skipped.

preservation$p_ocellata <- imported$p_oc_pres %>%
  mutate(
    sample_no = sample_no,
    decostand(select(., -sample_no), 
              method = "total",
              na.rm = "TRUE") * 100
  )

preservation$p_ocellata_tidy <- inner_join(imported$depths,
                                           preservation$p_ocellata) %>%
  na_if(0) %>%
  pivot_longer(cols = -c(1, 2), 
               names_to = "pres_quality",
               values_to = "pc_ab") %>%
  mutate(pres_quality = as_factor(pres_quality))



