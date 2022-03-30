library(readxl)
library(tidyverse)


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


# Remove objects that won't be required from now on.
rm(path)
rm(read_then_csv)
