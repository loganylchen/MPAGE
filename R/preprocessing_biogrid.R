#' Preprocess BioGRID Database PPI Data
#'
#' Download and preprocess protein-protein interaction data from BioGRID database
#' using direct download from the BioGRID FTP server and standardized processing.
#'
#' @param species_code Character string specifying species (e.g., "9606" for human)
#' @param score_threshold Numeric minimum confidence score threshold (not used for BioGRID, included for consistency)
#' @param output_file Path to save processed RDS file
#' @param version Character string specifying BioGRID version (default: "4.4.248")
#' @param verbose Logical whether to print progress messages
#'
#' @return A data frame with columns: geneA, geneB, score, species_code
#'
#' @export
#'
#' @examples
#' # Get human PPI data
#' human_ppi <- preprocessing_biogrid("9606")
#'
#' # Get mouse PPI data
#' mouse_ppi <- preprocessing_biogrid("10090")
preprocessing_biogrid <- function(species_code = "9606",
                                  score_threshold = 0,
                                  output_file = NULL,
                                  version = "4.4.248",
                                  verbose = TRUE) {
  # Check for required packages
  if (!requireNamespace("dplyr", quietly = TRUE)) {
    stop("dplyr package is required but not installed. Please install with: install.packages('dplyr')")
  }
  if (!requireNamespace("readr", quietly = TRUE)) {
    stop("readr package is required but not installed. Please install with: install.packages('readr')")
  }

  if (verbose) {
    message(sprintf("Downloading BioGRID v%s data for species %s", version, species_code))
  }

  # BioGRID download URL
  biogrid_url <- sprintf("https://downloads.thebiogrid.org/Download/BioGRID/Release-Archive/BIOGRID-%s/BIOGRID-ALL-%s.tab3.zip", version, version)

  tryCatch(
    {
      # Create temporary file for download
      temp_zip <- tempfile(fileext = ".zip")
      extract_dir <- tempfile()

      if (verbose) {
        message("Downloading from: ", biogrid_url)
      }

      # Download BioGRID data
      utils::download.file(biogrid_url, temp_zip, quiet = !verbose, mode = "wb")

      if (verbose) {
        message("Download complete, extracting data...")
      }

      # Create extraction directory
      dir.create(extract_dir)

      # Extract zip file
      utils::unzip(temp_zip, exdir = extract_dir)

      # Find the extracted TSV file
      tsv_files <- list.files(extract_dir, pattern = "\\.tab3$", full.names = TRUE)
      if (length(tsv_files) == 0) {
        tsv_files <- list.files(extract_dir, pattern = "\\.tab3\\.txt$", full.names = TRUE)
      }
      if (length(tsv_files) == 0) {
        stop("No BioGRID TSV file found in downloaded archive")
      }

      biogrid_file <- tsv_files[1]

      if (verbose) {
        message("Processing BioGRID data...")
      }

      # Read BioGRID data
      biogrid_data <- data.table::fread(biogrid_file, sep = "\t", check.names = TRUE, quote = "")

      if (verbose) {
        message(sprintf("Loaded %d interactions from BioGRID", nrow(biogrid_data)))
      }

      # Process the data to required format
      # BioGRID format processing
      processed_data <- biogrid_data %>%
        dplyr::select(
          gene1 = Official.Symbol.Interactor.A,
          gene2 = Official.Symbol.Interactor.B,
          species_codeA = Organism.ID.Interactor.A,
          species_codeB = Organism.ID.Interactor.B,
          Experimental.System.Type
        ) %>%
        dplyr::mutate(species_code = species_code, est = factor(Experimental.System.Type, levels = c("physical", "genetic"))) %>%
        dplyr::filter(species_codeA == species_codeB, species_codeA == species_code) %>%
        arrange(est) %>%
        distinct(gene1, gene2, .keep_all = TRUE) %>%
        dplyr::select(gene1, gene2, species_code, est)

      # Clean up temporary files
      unlink(temp_zip)
      unlink(extract_dir, recursive = TRUE)

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
    },
    error = function(e) {
      stop("Error processing BioGRID data: ", e$message)
    }
  )
}
