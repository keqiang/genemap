#' @import data.table
#' @importFrom magrittr %>%

#' @export
SPECIES <- c(
  "Human" = "hs",
  "Mouse" = "mm",
  "Rat" = "rn",
  "Zebrafish" = "dr",
  "Fruitfly" = "dm"
)

#' Get the mapping dataset for specified species
#'
#' @param from_species From which species
#' @param to_species To which species
#' @param version Ensembl version
#' @return A data.table that has all the gene mappings for specied parameters
get_mapping_table <- function(from_species, to_species, version = c("es96")) {
  version <- match.arg(version)
  file_type <- ifelse(from_species == to_species, "desc", "ortholog")
  mapping_direction <- ifelse(
    from_species == to_species, # if mapping within same species
    from_species,
    stringr::str_c(from_species, "2", to_species)
  )
  mapping_dataset_name <- stringr::str_interp("${version}_${mapping_direction}_${file_type}_ids")
  if (!mapping_dataset_name %in% names(gene_mapping_datasets)) {
    stop("Wrong parameters")
  }
  gene_mapping_datasets[[mapping_dataset_name]]
}

#' Map gene identifiers within or cross species
#'
#' @param from_species From which species
#' @param from_id_list Gene identifiers to be mapped
#' @param from_id_type The type of the identifiers in \code{from_id_list}
#' @param to_species To which species
#' @param to_id_type The type of the identifiers to map to
#' @return A named list of the same size with \code{from_id_list}. Each gene from \code{from_id_list} will be mapped to an element in this result list.
#'     \code{$mapped_genes} will contain the mapped identifiers and  \code{$mapped_gene_symbols} will contain corresponding symbols of these mapped genes.
#'     \code{$gene_symbol} will be the corresponding gene symbol of the input gene identifier. Name of each element is the original gene identifier used to map.
#' @examples
#' # within species mapping
#' map_genes(
#'   from_species = "hs",
#'   from_id_list = c("ENPP4", "GCLC"),
#'   from_id_type = "symbol",
#'   to_species = "hs"
#' )
#' # cross species mapping
#' map_genes(
#'   from_species = "hs",
#'   from_id_list = c("ENPP4", "GCLC"),
#'   from_id_type = "symbol",
#'   to_species = "mm",
#'   to_id_type = "symbol"
#' )
#' @export
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

  from_id_type_col <- stringr::str_c(from_id_type, from_species, sep = "_")
  to_id_type_col <- stringr::str_c(to_id_type, to_species, sep = "_")

  if (from_id_type == "ensemblgid") {
    from_id_list <- from_id_list %>%
      stringr::str_split("\\.", n = 2) %>%
      purrr::map_chr(~ .[[1]])
  }

  # produce the result
  result_cols <- c(from_id_type_col, to_id_type_col)
  from_symbol_col <- stringr::str_c("symbol", from_species, sep = "_")
  if (from_id_type != "symbol") {
    result_cols <- c(result_cols, from_symbol_col) # always add symbol
  }

  to_symbol_col <- stringr::str_c("symbol", to_species, sep = "_")
  if (to_id_type != "symbol") {
    result_cols <- c(result_cols, to_symbol_col) # always add symbol
  }
  mapped_genes <- mapping_table[get(from_id_type_col) %in% from_id_list, ..result_cols]

  result_cols <- result_cols[-1]
  from_id_list %>%
    purrr::map(function(id) {
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
    purrr::set_names(from_id_list)
}
