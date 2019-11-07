# Package `genemap`
Map gene identifies within or across species

## Usage
```
map_genes(
  from_species = SPECIES,
  from_id_list,
  from_id_type = c("ensemblgid", "ncbigid", "symbol"),
  to_species = SPECIES,
  to_id_type = c("ensemblgid", "ncbigid", "symbol")
)

# Type genemap::SPECIES to see the full list of supported species
genemap::SPECIES
  Human     Mouse       Rat Zebrafish  Fruitfly
   "hs"      "mm"      "rn"      "dr"      "dm"
```

## Arguments
`from_species` From which species

`from_id_list` Gene identifiers to be mapped

`from_id_type` The type of the identifiers in from_id_list

`to_species` To which species

`to_id_type` The type of the identifiers to map to

## Examples

### Same species mapping
```{r}
map_genes(
  from_species = "hs",
  from_id_list = c("ENPP4", "GCLC"),
  from_id_type = "symbol",
  to_species = "hs"
)

$ENPP4
$ENPP4$mapped_genes
[1] "ENSG00000001561"

$ENPP4$mapped_gene_symbols
[1] "ENPP4"

$ENPP4$gene_symbol
[1] "ENPP4"


$GCLC
$GCLC$mapped_genes
[1] "ENSG00000001084"

$GCLC$mapped_gene_symbols
[1] "GCLC"

$GCLC$gene_symbol
[1] "GCLC"
```

### Cross-species mapping
```{r}
map_genes(
  from_species = "hs",
  from_id_list = c("ENPP4", "GCLC"),
  from_id_type = "symbol",
  to_species = "mm",
  to_id_type = "symbol"
)

$ENPP4
$ENPP4$mapped_genes
[1] "Enpp4"


$GCLC
$GCLC$mapped_genes
[1] "Gclc"
```
