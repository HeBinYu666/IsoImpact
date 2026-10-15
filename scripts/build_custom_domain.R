#!/usr/bin/env Rscript

library(AnnotationHub)
library(ensembldb)
library(GenomicRanges)
library(dplyr)

usage <- function() {
  cat(
    "Usage:\n",
    "  Rscript scripts/build_custom_domain.R --pfam pfam_results.txt --out domain.csv --ensembl 110\n\n",
    "Options:\n",
    "  --pfam      PfamScan result table\n",
    "  --out       Output domain-coordinate file used by IsoImpact with -d/--domain\n",
    "  --ensembl   Ensembl release version number for the AnnotationHub EnsDb query. Default: 110\n",
    "  --species   Species name for the AnnotationHub EnsDb query. Default: Homo sapiens\n",
    "  --gtf       Optional matching GTF path for record keeping; the mapping uses AnnotationHub EnsDb\n",
    "  --help      Show this message\n",
    sep = ""
  )
}

args <- commandArgs(trailingOnly = TRUE)

if ("--help" %in% args || length(args) == 0) {
  usage()
  quit(status = 0)
}

get_arg <- function(flag, default = NULL) {
  idx <- match(flag, args)
  if (is.na(idx)) return(default)
  if (idx == length(args)) stop(paste("Missing value for", flag), call. = FALSE)
  args[[idx + 1]]
}

# 1. Read input and output paths.
pfam_txt <- get_arg("--pfam")
output_csv <- get_arg("--out")
ensembl_release <- get_arg("--ensembl", "110")
species <- get_arg("--species", "Homo sapiens")
gtf_file <- get_arg("--gtf")

if (is.null(pfam_txt) || is.null(output_csv)) {
  usage()
  stop("--pfam and --out are required.", call. = FALSE)
}

if (!file.exists(pfam_txt)) {
  stop(paste("PfamScan result file not found:", pfam_txt), call. = FALSE)
}

if (!is.null(gtf_file)) {
  message("[IsoImpact] Matching GTF file recorded for downstream analysis: ", gtf_file)
}

# 2. Retrieve the EnsDb annotation database for the selected Ensembl release.
message("[IsoImpact] Querying AnnotationHub for ", species, " Ensembl release ", ensembl_release, "...")
ah <- AnnotationHub()
query_result <- query(ah, c("EnsDb", species, ensembl_release))

if (length(query_result) == 0) {
  stop(
    paste(
      "No matching EnsDb record was found in AnnotationHub for",
      species,
      "Ensembl release",
      ensembl_release,
      sep = " "
    ),
    call. = FALSE
  )
}

edb <- query_result[[1]]

# 3. Read PfamScan results.
message("[IsoImpact] Reading PfamScan results...")
pfam_data <- read.table(pfam_txt, stringsAsFactors = FALSE, comment.char = "#", fill = TRUE)

if (nrow(pfam_data) == 0) {
  stop("The PfamScan result file is empty or could not be read.", call. = FALSE)
}

if (ncol(pfam_data) < 7) {
  stop("The PfamScan result file has fewer than seven columns.", call. = FALSE)
}

message("[IsoImpact] Pfam domain records detected: ", nrow(pfam_data))

results <- list()
valid_count <- 0

# 4. Map protein-domain intervals to genomic coordinates.
for (i in 1:nrow(pfam_data)) {
  row <- pfam_data[i, ]

  prot_id_raw <- as.character(row$V1)
  protein_id <- strsplit(prot_id_raw, "\\.")[[1]][1]

  pfam_start <- as.numeric(row$V2)
  pfam_end <- as.numeric(row$V3)
  domain_name <- as.character(row$V7)

  pfam_id_raw <- as.character(row$V6)
  domain_id <- strsplit(pfam_id_raw, "\\.")[[1]][1]

  tryCatch({
    protein_range <- IRanges(start = pfam_start, end = pfam_end, names = protein_id)
    genomic_coords <- proteinToGenome(protein_range, db = edb)

    if (length(genomic_coords) > 0 && length(genomic_coords[[1]]) > 0) {
      g_start <- min(start(genomic_coords[[1]]))
      g_end <- max(end(genomic_coords[[1]]))

      valid_count <- valid_count + 1
      results[[valid_count]] <- data.frame(
        Protein_ID = protein_id,
        Domain_ID = domain_id,
        Domain_Name = domain_name,
        Genomic_Start = g_start,
        Genomic_End = g_end,
        stringsAsFactors = FALSE
      )
    }
  }, error = function(e) {
    # Skip domains that cannot be mapped by the selected EnsDb annotation.
  })

  if (i %% 1000 == 0) {
    message("[IsoImpact] Processed ", i, " records; mapped domains: ", valid_count)
  }
}

# 5. Write the IsoImpact-compatible domain-coordinate file.
if (valid_count > 0) {
  final_df <- bind_rows(results)
  final_df <- distinct(final_df, Protein_ID, Domain_ID, Genomic_Start, Genomic_End, .keep_all = TRUE)
  write.csv(final_df, output_csv, row.names = FALSE)
  message("[IsoImpact] Domain-coordinate file generated: ", output_csv)
  message("[IsoImpact] Mapped non-redundant domains: ", nrow(final_df))
} else {
  warning("No domains could be mapped. Please check whether the protein IDs and annotation release match.")
}
