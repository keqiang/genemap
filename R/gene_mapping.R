SPECIES <- c("hs", "mm", "rn", "dr", "dz")

get_mapping_table <- memoise::memoise(
  function(from_species, to_species, version = c("es96")) {
    version <- match.arg(version)
    file_type <- ifelse(from_species == to_species, "desc", "ortholog")
    mapping_direction <- ifelse(
      from_species == to_species, # if mapping within same species
      from_species,
      str_c(from_species, "2", to_species)
    )
    ontholog_file_path <- file.path(
      hdrf_data_dir,
      "speciesids",
      str_interp("${version}_${mapping_direction}_${file_type}_ids")
    )

    mapping_table <- fread(
      ontholog_file_path,
      header = TRUE
    )
  }
)

map_genes <- function(from_species = SPECIES,
                      from_id_list,
                      from_id_type = c("ensemblgid", "ncbigid", "symbol"),
                      to_species = SPECIES,
                      to_id_type = c("ensemblgid", "ncbigid", "symbol")) {
  # validate the arguments
  if (is.null(from_id_list) || length(from_id_list) <= 0) {
    return(
      list(
        "mapped_genes" = c()
      )
    )
  }

  from_species <- match.arg(from_species)
  from_id_type <- match.arg(from_id_type)
  to_species <- match.arg(to_species)
  to_id_type <- match.arg(to_id_type)

  mapping_table <- get_mapping_table(from_species, to_species)

  from_id_type_col <- str_c(from_id_type, from_species, sep = "_")
  to_id_type_col <- str_c(to_id_type, to_species, sep = "_")

  if (from_id_type == "ensemblgid") {
    from_id_list <- from_id_list %>%
      str_split("\\.", n = 2) %>%
      map_chr(~ .[[1]])
  }

  # produce the result
  result_cols <- c(from_id_type_col, to_id_type_col)
  from_symbol_col <- str_c("symbol", from_species, sep = "_")
  if (from_id_type != "symbol") {
    result_cols <- c(result_cols, from_symbol_col) # always add symbol
  }

  to_symbol_col <- str_c("symbol", to_species, sep = "_")
  if (to_id_type != "symbol") {
    result_cols <- c(result_cols, to_symbol_col) # always add symbol
  }
  mapped_genes <- mapping_table[get(from_id_type_col) %in% from_id_list, ..result_cols]

  result_cols <- result_cols[-1]
  mapped_list <- from_id_list %>%
    map(function(id) {
      cur_id_mapped_genes <- mapped_genes[get(from_id_type_col) == id, ..result_cols]
      res <- list(
        "mapped_genes" = cur_id_mapped_genes[[to_id_type_col]]
      )
      if (to_id_type != "symbol") {
        res$mapped_gene_symbols <- cur_id_mapped_genes[[to_symbol_col]]
      }
      if (nrow(cur_id_mapped_genes) > 0) {
        res$gene_symbol <- cur_id_mapped_genes[[from_symbol_col]][1]
      }
      res
    }) %>%
    set_names(from_id_list)
}
