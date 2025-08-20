#' Preprocess IntAct Database PPI Data
#'
#' Download and preprocess protein-protein interaction data from IntAct database
#' using direct download from the EBI FTP server and standardized processing.
#'
#' @param species_code Character string specifying species (e.g., "9606" for human)
#' @param score_threshold Numeric minimum confidence score threshold (0-1)
#' @param output_file Path to save processed RDS file
#' @param version Character string specifying IntAct version (default: "current")
#' @param verbose Logical whether to print progress messages
#'
#' @return A data frame with columns: geneA, geneB, score, species_code
#'
#' @export
#'
#' @examples
#' # Get human PPI data
#' human_ppi <- preprocessing_intact("9606", score_threshold = 0.4)
#'
#' # Get mouse PPI data
#' mouse_ppi <- preprocessing_intact("10090", score_threshold = 0.4)
preprocessing_intact <- function(species_code = "9606",
                                 output_file = NULL,
                                 version = "2025-03-28",
                                 verbose = TRUE) {
  # Check for required packages
  if (!requireNamespace("dplyr", quietly = TRUE)) {
    stop("dplyr package is required but not installed. Please install with: install.packages('dplyr')")
  }
  if (!requireNamespace("data.table", quietly = TRUE)) {
    stop("data.table package is required but not installed. Please install with: install.packages('data.table')")
  }
  if (!requireNamespace("stringr", quietly = TRUE)) {
    stop("stringr package is required but not installed. Please install with: install.packages('stringr')")
  }

  if (verbose) {
    message(sprintf("Downloading IntAct data (Version: %s)for species %s", version, species_code))
  }

  # IntAct download URL
  intact_url <- sprintf("https://ftp.ebi.ac.uk/pub/databases/intact/%s/psimitab/intact.txt", version)


  # Create temporary file for download
  temp_file <- tempfile(fileext = ".txt")

  if (verbose) {
    message("Downloading from: ", intact_url)
    message("To the file: ", temp_file)
  }

  # Download IntAct data
  utils::download.file(intact_url, temp_file, quiet = !verbose, mode = "w")

  if (verbose) {
    message("Download complete, processing data...")
  }

  # Read IntAct data (large file, read in chunks)
  intact_data <- data.table::fread(temp_file, sep = "\t", check.names = TRUE, quote = "")

  if (verbose) {
    message(sprintf("Loaded %d interactions from IntAct", nrow(intact_data)))
  }

  # Process the data to required format
  # IntAct PSI-MITAB format processing
  processed_data <- intact_data %>%
    dplyr::select(
      Alt..ID.s..interactor.A,
      Alt..ID.s..interactor.B,
      Alias.es..interactor.A,
      Alias.es..interactor.B,
      Taxid.interactor.A,
      Taxid.interactor.B,
      Confidence.value.s.
    ) %>%
    # Extract species information
    dplyr::mutate(
      taxidA = stringr::str_extract(Taxid.interactor.A, "taxid:([0-9]+)") %>%
        stringr::str_remove("taxid:"),
      taxidB = stringr::str_extract(Taxid.interactor.B, "taxid:([0-9]+)") %>%
        stringr::str_remove("taxid:")
    ) %>%
    # Filter by species
    dplyr::filter(
      !is.na(taxidA), !is.na(taxidB),
      taxidA == species_code,
      taxidB == species_code,
      taxidA == taxidB
    ) %>%
    # Extract gene names using multiple patterns
    dplyr::mutate(
      # Extract from gene name annotations
      geneA_name = stringr::str_extract(
        Alias.es..interactor.A,
        "uniprotkb:[A-Za-z0-9-]+\\(gene name\\)"
      ) %>%
        stringr::str_remove("uniprotkb:") %>%
        stringr::str_remove("\\(gene name\\)") %>%
        toupper(),
      geneB_name = stringr::str_extract(
        Alias.es..interactor.B,
        "uniprotkb:[A-Za-z0-9-]+\\(gene name\\)"
      ) %>%
        stringr::str_remove("uniprotkb:") %>%
        stringr::str_remove("\\(gene name\\)") %>%
        toupper(),

      # Extract from display names as fallback
      geneA_display = stringr::str_extract(
        Alias.es..interactor.A,
        "display_short\\):([^|]+)"
      ) %>%
        stringr::str_remove("display_short\\):"),
      geneB_display = stringr::str_extract(
        Alias.es..interactor.B,
        "display_short\\):([^|]+)"
      ) %>%
        stringr::str_remove("display_short\\):"),

      # Use gene name if available, otherwise display name
      gene1 = dplyr::coalesce(geneA_name, geneA_display),
      gene2 = dplyr::coalesce(geneB_name, geneB_display),

      # Extract interaction score
      score = stringr::str_extract(Confidence.value.s., "intact-miscore:([0-9\\.]+)") %>%
        stringr::str_remove("intact-miscore:") %>%
        as.numeric(),
      species_code = taxidA
    ) %>%
    # Clean and filter
    dplyr::select(gene1, gene2, score, species_code) %>%
    dplyr::filter(
      !is.na(gene1), !is.na(gene2), gene1 != "", gene2 != "", gene1 != gene2,
      !is.na(score),
    ) %>%
    dplyr::distinct() %>%
    dplyr::arrange(dplyr::desc(score)) %>%
    dplyr::select(gene1, gene2, species_code, score)

  # Clean up temporary file
  unlink(temp_file)

  if (verbose) {
    message(sprintf("Processed %d unique interactions", nrow(processed_data)))
    message(sprintf("Unique proteins: %d", length(unique(c(processed_data$gene1, processed_data$gene2)))))
  }

  # Save if output file specified
  if (!is.null(output_file)) {
    saveRDS(processed_data, file = output_file)
    if (verbose) {
      message(sprintf("Data saved to: %s", output_file))
    }
  }
  return(processed_data)
}
