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