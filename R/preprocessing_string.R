#' Preprocess STRING Database PPI Data
#'
#' Download and preprocess protein-protein interaction data from STRING database
#' using the StringDB package for direct API access.
#'
#' @param species_code Character string specifying species (e.g., "9606" for human)
#' @param score_threshold Numeric minimum confidence score threshold (0-1000)
#' @param output_file Path to save processed RDS file
#' @param version Character string specifying STRING version (default: "12.0")
#' @param verbose Logical whether to print progress messages
#'
#' @return A data frame with columns: geneA, geneB, score, species_code
#'
#' @export
#'
#' @examples
#' # Get human PPI data
#' human_ppi <- preprocessing_string("9606", score_threshold = 400)
#'
#' # Get mouse PPI data
#' mouse_ppi <- preprocessing_string("10090", score_threshold = 400)
preprocessing_string <- function(species_code = 9606,
                                 score_threshold = 0,
                                 output_file = NULL,
                                 version = "12",
                                 verbose = TRUE) {
  # Check for required packages
  if (!requireNamespace("STRINGdb", quietly = TRUE)) {
    stop("STRINGdb package is required but not installed. Please install with: install.packages('STRINGdb')")
  }
  if (!requireNamespace("dplyr", quietly = TRUE)) {
    stop("dplyr package is required but not installed. Please install with: install.packages('dplyr')")
  }

  if (verbose) {
    message(sprintf("Downloading STRING v%s data for species %s", version, species_code))
  }


  # Initialize StringDB object
  string_db <- STRINGdb::STRINGdb$new(
    species = species_code,
    version = version,
    score_threshold = score_threshold
  )

  if (verbose) {
    message("STRING database initialized successfully")
  }

  # Get all interactions with scores above threshold
  stringdb_proteins <- string_db$get_proteins()
  # Get all interactions for the species above threshold
  interactions <- string_db$get_interactions(
    stringdb_proteins$protein_external_id
  )

  if (verbose) {
    message(sprintf("Retrieved %d interactions from STRING", nrow(interactions)))
  }
  interactions$gene1 <- stringdb_proteins$preferred_name[match(interactions$from, stringdb_proteins$protein_external_id)]
  interactions$gene2 <- stringdb_proteins$preferred_name[match(interactions$to, stringdb_proteins$protein_external_id)]



  # Process the data to required format
  processed_data <- interactions %>%
    dplyr::mutate(
      # Ensure species_code is character
      species_code = as.character(species_code),
      # Ensure score is numeric
      score = as.numeric(combined_score)
    ) %>%
    dplyr::distinct() %>%
    dplyr::filter(
      !is.na(gene1),
      !is.na(gene2),
      gene1 != gene2,
      gene1 != "",
      gene2 != ""
    ) %>%
    dplyr::select(gene1, gene2, species_code, score)

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
