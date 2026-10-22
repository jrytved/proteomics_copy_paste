# proteomics-copy-paste
Useful function and scripts that are not important enough to structure into a standalone package

## R

### read_gmt.R

```R
read_gmt <- function(path) {
  lines  <- readLines(path, warn = FALSE)
  fields <- strsplit(lines, "\t", fixed = TRUE)
  fields <- fields[lengths(fields) >= 3]
  sets   <- lapply(fields, function(x) { g <- x[-(1:2)]; g[g != ""] })
  setNames(sets, vapply(fields, `[`, "", 1))
}

```

read_uniprot_fasta.R

## Python
swissprot_hgnc_map.py




