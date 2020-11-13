# 01 Load data
# ————————————
# This script loads all of the data required for the analysis of the new diatom
# record of the I-284 core (spanning MIS 7 to 9). The data was inputted into 
# several worksheets in the excel spreadsheet "diatom-data.xlsx". This script 
# exports the data as csv files for longevity. All paths are relative to the 
# wd which is the project root (phd folder).




# 1. Load packages --------------------------------------------------------

library(readxl)
library(purrr)
library(readr)
library(tools)




# 2. Read in data ---------------------------------------------------------

# Old manual way. No longer needed. Skip to next section. Left here for future
# reference if needed.

# path <- "data/diatom-data.xlsx"
# depths <- read_xlsx(path, "depths")
# mass <- read_xlsx(path, "mass")
# microspheres <- read_xlsx(path, "microspheres")
# counts <- read_xlsx(path, "counts")
# p_oc_breakdown <- read_xlsx(path, "p_oc_breakdown")
# p_oc_pres <- read_xlsx(path, "p_oc_pres")
# f_index <- read_xlsx(path, "f_index")
# lifestyle <- read_xlsx(path, "lifestyle")




# 3. Load data and write to csv -------------------------------------------

# Create function that reads in a worksheet from an excel file then writes
# it to csv.
read_then_csv <- function(sheet, path) {
  path %>%
    read_excel(sheet = sheet) %>%
    write_csv(paste0("data/csv/imported-", sheet, ".csv"))
}


# Execute function and iterate over all worksheets in excel file.
path <- "data/ioannina.xlsx"
imported <- path %>%
  excel_sheets() %>%
  set_names() %>%
  map(read_then_csv, path = path)




# 4. Clean up -------------------------------------------------------------

# Remove objects that won't be required from now on.
rm(path)
rm(read_then_csv)



