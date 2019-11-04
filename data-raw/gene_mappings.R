## code to prepare `genemappings` dataset goes here
library(purrr)

mapping_files <- list.files("mapping_data") # full path of original mapping files
dataset_names <- basename(mapping_files) # use the file names as mapping dataset names

gene_mappings <- mapping_files %>% map(function(f) {
  file.path("mapping_data", f) %>% data.table::fread(header = TRUE)
}) %>% set_names(dataset_names)

usethis::use_data(gene_mappings, internal = TRUE)
