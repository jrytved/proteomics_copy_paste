# proteomics-copy-paste
Useful function and scripts that are not important enough to structure into a standalone package

## R


### tabulate_nicely.R
Any data to a reasonably aesthetically pleasing table

tabulate_nicely <- function(table, header.md, note.md){
  
    table  %>%
      gt::gt() %>%
      gt::tab_header(title = gt::md(header.md)) %>%
      gt::opt_table_font(font = "Arial") %>%
      gt::opt_table_lines(extent = "none") %>%
      gt::tab_source_note(source_note = gt::md(note.md))%>%
      gt::tab_options(
        table.font.size                   = gt::px(12),
        heading.align                     = "left",
        heading.title.font.size           = gt::px(13),
        column_labels.font.weight         = "bold",
        column_labels.border.top.style    = "solid",      
        column_labels.border.top.width    = gt::px(2),
        column_labels.border.top.color    = "black",
        column_labels.border.bottom.style = "solid",        
        column_labels.border.bottom.width = gt::px(1),
        column_labels.border.bottom.color = "black",
        table_body.border.bottom.style    = "solid",      
        table_body.border.bottom.width    = gt::px(2),
        table_body.border.bottom.color    = "black",
        data_row.padding                  = gt::px(3),
        footnotes.font.size               = gt::px(10),
        source_notes.font.size            = gt::px(10)
      ) %>%
     gt::cols_align(align = "right")
  
}

### read_gmt.R

Read a .gmt file into a list of named sets.

```R
read_gmt <- function(path) {
  lines  <- readLines(path, warn = FALSE)
  fields <- strsplit(lines, "\t", fixed = TRUE)
  fields <- fields[lengths(fields) >= 3]
  sets   <- lapply(fields, function(x) { g <- x[-(1:2)]; g[g != ""] })
  setNames(sets, vapply(fields, `[`, "", 1))
}

```

### read_uniprot_fasta.R

Read a fasta retrieved from UniProt (SwissProt reviewd) into a dataframe so sequences can easily be retrieved and mapped to identified peptides.
Adheres to official UniProt fasta structure.

```R
#' Read a UniProt FASTA file into a data.frame
#'
#' Reads the whole file as a single string and splits it into records with
#' Uniprot headers ("sp|A0A087X1C5|CP2D7_HUMAN") are split into
#' database, accession, Protein.Names and description. Ids that don't follow
#' the db|accession|name pattern get NA in those three columns.
#'
#' @param path Path to a FASTA file (plain or .gz).
#' @return data.frame with columns id, database, accession, Protein.Names,
#'         description, sequence (upper case), length.
read_uniprot_fasta <- function(path) {
  if (!file.exists(path)) stop("File not found: ", path)

  # Read file into one string
  if (grepl("\\.gz$", path, ignore.case = TRUE)) {
    con <- gzfile(path, "r")
    on.exit(close(con), add = TRUE)
    txt <- paste(readLines(con, warn = FALSE), collapse = "\n")
  } else {
    txt <- readChar(path, file.info(path)$size, useBytes = TRUE)
  }
  if (grepl("\r", txt, fixed = TRUE)) txt <- gsub("\r", "", txt, fixed = TRUE)

  # Drop anything before the first '>'
  first <- regexpr(">", txt, fixed = TRUE)
  if (first < 1L) stop("Not a valid FASTA: no '>' header found")
  txt <- substr(txt, first + 1L, nchar(txt, type = "bytes"))

  # Split into records
  recs <- strsplit(txt, "\n>", fixed = TRUE)[[1L]]
  rm(txt)

  # Split nto header and sequence
  nl      <- regexpr("\n", recs, fixed = TRUE)
  has_seq <- nl > 0L
  headers <- recs
  headers[has_seq] <- substr(recs[has_seq], 1L, nl[has_seq] - 1L)

  sequences <- character(length(recs))
  sequences[has_seq] <- substr(recs[has_seq], nl[has_seq] + 1L,
                               nchar(recs[has_seq], type = "bytes"))
  rm(recs)
  
  sequences <- toupper(gsub("[ \t\n]", "", sequences, perl = TRUE))
  headers   <- sub("[ \t]+$", "", headers, perl = TRUE)

  # Split header into id, database etc.
  ws   <- regexpr("[ \t]", headers, perl = TRUE)
  has  <- ws > 0L
  id   <- headers
  desc <- character(length(headers))
  id[has]   <- substr(headers[has], 1L, ws[has] - 1L)
  desc[has] <- substr(headers[has], ws[has] + 1L, nchar(headers[has]))

  pat <- "^([^|]*)\\|([^|]*)\\|(.*)$"
  ok  <- grepl(pat, id, perl = TRUE)
  database <- accession <- prot_name <- rep(NA_character_, length(id))
  database[ok]  <- sub(pat, "\\1", id[ok], perl = TRUE)
  accession[ok] <- sub(pat, "\\2", id[ok], perl = TRUE)
  prot_name[ok] <- sub(pat, "\\3", id[ok], perl = TRUE)

  data.frame(id            = id,
             database      = database,
             accession     = accession,
             Protein.Names = prot_name,
             description   = desc,
             sequence      = sequences,
             length        = nchar(sequences),
             stringsAsFactors = FALSE)
}
```


## Python
### swissprot_hgnc_map.py
Stand-alone script to build a map between SwissProt/UniProt accessions and HGNC identifiers from official and updated resources.

### kinase_library_predict.py
Score upstream kinases for a table of identified phosphosites.

```python
def score_dia_nn_sites(sites, out_folder: str, seq_col = "Sequence"):

    pobj = kl.PhosphoProteomics(sites, seq_col=seq_col)

    path = Path(out_folder)
    ser_thr_path = path / "ser_thr_predictions.tsv"
    tyr_path = path / "tyr_predictions.tsv"
    
    ser_thr_predictions = pobj.predict(kin_type='ser_thr')
    tyrosine_predictions = pobj.predict(kin_type='tyrosine')

    ser_thr_predictions.to_csv(ser_thr_path, sep="\t", index=False)
    tyrosine_predictions.to_csv(tyr_path, sep="\t", index=False)
```

Then pivot it longer when importing into R...

```R
ser_thr_kinase_preds_long <- ser_thr_kinase_preds %>%
  pivot_longer(
    cols = matches("_(score|score_rank|percentile|percentile_rank)$"),
    names_to = c("kinase", ".value"),
    names_pattern = "(.*)_(score|score_rank|percentile|percentile_rank)$",
    values_drop_na = TRUE
  )
```




