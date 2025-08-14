# IntAct Database Preprocessing Functions
#'
#' Functions for preprocessing and accessing IntAct PPI database

#' Preprocess IntAct Database
#'
#' Process raw IntAct TSV file into clean interaction data
#'
#' @param input_file Path to IntAct TSV file
#' @param output_file Path to save processed RDS file (optional)
#' @param species_filter Vector of NCBI taxonomy IDs to filter (default: c("9606", "10090"))
#' @param min_score Minimum interaction score threshold (default: 0.0)
#' @param max_missing_proportion Maximum proportion of missing gene names to allow (default: 0.1)
#'
#' @return Cleaned interaction data frame with columns: geneA, geneB, score, taxid
#'
#' @export
#'
#' @examples
#' \dontrun{
#' # Process human and mouse data
#' intact_data <- preprocess_intact_data(
#'   input_file = "intact.txt",
#'   species_filter = c("9606", "10090"),
#'   min_score = 0.3
#' )
#' }
preprocess_intact_data <- function(input_file,
                                   output_file = NULL,
                                   species_filter = c("9606", "10090"),
                                   min_score = 0.0,
                                   max_missing_proportion = 0.1) {
  if (!file.exists(input_file)) {
    stop("Input file not found: ", input_file)
  }

  message("Loading IntAct data from: ", input_file)

  # Read IntAct data with proper handling
  tryCatch(
    {
      df <- data.table::fread(input_file, sep = "\t", check.names = TRUE, quote = "")
    },
    error = function(e) {
      stop("Error reading IntAct file: ", e$message)
    }
  )

  message("Processing ", nrow(df), " interactions...")

  # Validate required columns exist
  required_cols <- c(
    "Alt..ID.s..interactor.A", "Alias.es..interactor.A",
    "Taxid.interactor.A", "Taxid.interactor.B",
    "Confidence.value.s."
  )

  missing_cols <- setdiff(required_cols, colnames(df))
  if (length(missing_cols) > 0) {
    stop("Missing required columns: ", paste(missing_cols, collapse = ", "))
  }

  # Process the data
  processed_data <- df %>%
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
      taxidA %in% species_filter,
      taxidB %in% species_filter,
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
      geneA = dplyr::coalesce(geneA_name, geneA_display),
      geneB = dplyr::coalesce(geneB_name, geneB_display),

      # Extract interaction score
      score = stringr::str_extract(Confidence.value.s., "intact-miscore:([0-9\\.]+)") %>%
        stringr::str_remove("intact-miscore:") %>%
        as.numeric(),
      taxid = taxidA
    ) %>%
    # Clean and filter
    dplyr::select(geneA, geneB, score, taxid) %>%
    dplyr::filter(
      !is.na(geneA), !is.na(geneB), geneA != "", geneB != "",
      !is.na(score), score >= min_score
    ) %>%
    dplyr::distinct() %>%
    dplyr::arrange(dplyr::desc(score))

  # Calculate missing proportion
  total_interactions <- nrow(processed_data)
  missing_genes <- sum(is.na(processed_data$geneA) | is.na(processed_data$geneB))
  missing_proportion <- missing_genes / total_interactions

  if (missing_proportion > max_missing_proportion) {
    warning(sprintf(
      "High proportion of missing gene names: %.2f%%",
      missing_proportion * 100
    ))
  }

  # Create summary statistics
  cat("=== IntAct Processing Summary ===\n")
  cat(sprintf("Total interactions: %d\n", nrow(processed_data)))
  cat(sprintf(
    "Unique proteins: %d\n",
    length(unique(c(processed_data$geneA, processed_data$geneB)))
  ))
  cat(sprintf(
    "Species: %s\n",
    paste(unique(processed_data$taxid), collapse = ", ")
  ))
  cat(sprintf(
    "Score range: [%.3f, %.3f]\n",
    min(processed_data$score, na.rm = TRUE),
    max(processed_data$score, na.rm = TRUE)
  ))
  cat(sprintf("Missing gene names: %.2f%%\n", missing_proportion * 100))

  # Species breakdown
  species_summary <- processed_data %>%
    dplyr::group_by(taxid) %>%
    dplyr::summarise(
      interactions = n(),
      unique_proteins = length(unique(c(geneA, geneB))),
      .groups = "drop"
    ) %>%
    dplyr::arrange(dplyr::desc(interactions))

  cat("\n=== Species Breakdown ===\n")
  print(species_summary)

  # Save if output file provided
  if (!is.null(output_file)) {
    saveRDS(processed_data, output_file)
    message("Saved processed data to: ", output_file)
  }

  return(processed_data)
}

#' Load IntAct Network
#'
#' Load preprocessed IntAct data and create network
#'
#' @param processed_file Path to processed RDS file
#' @param min_score Minimum interaction score to include
#' @param species Filter by species taxonomy ID
#' @param proteins Optional vector of proteins to extract subnetwork
#'
#' @return igraph object representing the PPI network
#'
#' @export
#'
#' @examples
#' \dontrun{
#' # Load full human network
#' intact_network <- load_intact_network(
#'   processed_file = "intact_processed.rds",
#'   species = "9606",
#'   min_score = 0.3
#' )
#' }
load_intact_network <- function(processed_file,
                                min_score = 0.0,
                                species = NULL,
                                proteins = NULL) {
  if (!file.exists(processed_file)) {
    stop("Processed file not found: ", processed_file)
  }

  # Load processed data
  intact_data <- readRDS(processed_file)

  # Filter by score and species
  filtered_data <- intact_data %>%
    dplyr::filter(score >= min_score)

  if (!is.null(species)) {
    filtered_data <- filtered_data %>%
      dplyr::filter(taxid == as.character(species))
  }

  if (nrow(filtered_data) == 0) {
    warning("No interactions found with specified criteria")
    return(igraph::make_empty_graph())
  }

  # Create igraph object
  net <- igraph::graph_from_data_frame(
    filtered_data[, c("geneA", "geneB")],
    directed = FALSE
  )

  # Add edge attributes
  igraph::E(net)$weight <- filtered_data$score
  igraph::E(net)$source <- "IntAct"

  # Extract subnetwork if proteins provided
  if (!is.null(proteins)) {
    available_proteins <- V(net)$name
    valid_proteins <- proteins[proteins %in% available_proteins]

    if (length(valid_proteins) > 0) {
      net <- igraph::induced_subgraph(net, valid_proteins)
      message("Extracted subnetwork with ", length(valid_proteins), " proteins")
    } else {
      warning("None of the specified proteins found in the network")
      return(igraph::make_empty_graph())
    }
  }

  cat("=== IntAct Network Summary ===\n")
  cat(sprintf("Nodes: %d\n", length(igraph::V(net))))
  cat(sprintf("Edges: %d\n", length(igraph::E(net))))
  cat(sprintf("Average degree: %.2f\n", mean(igraph::degree(net))))
  cat(sprintf(
    "Score range: [%.3f, %.3f]\n",
    min(igraph::E(net)$weight, na.rm = TRUE),
    max(igraph::E(net)$weight, na.rm = TRUE)
  ))

  return(net)
}

#' Download and Process IntAct Data
#'
#' Download latest IntAct data and preprocess it
#'
#' @param output_dir Directory to save processed files
#' @param species_filter Species to include
#' @param min_score Minimum interaction score
#' @param overwrite Whether to overwrite existing files
#'
#' @return Path to processed file
#'
#' @export
#'
#' @examples
#' \dontrun{
#' # Download and process latest IntAct data
#' processed_file <- download_intact_data()
#' }
download_intact_data <- function(output_dir = tempdir(),
                                 species_filter = c("9606", "10090"),
                                 min_score = 0.0,
                                 overwrite = FALSE) {
  # URL for IntAct data (may need to be updated)
  intact_url <- "https://ftp.ebi.ac.uk/pub/databases/intact/current/psimitab/intact.zip"

  # Create output directory
  if (!dir.exists(output_dir)) {
    dir.create(output_dir, recursive = TRUE)
  }

  output_file <- file.path(output_dir, "intact_processed.rds")

  if (file.exists(output_file) & !overwrite) {
    message("Processed file already exists: ", output_file)
    return(output_file)
  }

  message("Downloading IntAct data...")

  # Download and extract (placeholder - actual implementation would download)
  temp_file <- file.path(tempdir(), "intact.txt")

  if (!file.exists(temp_file)) {
    warning("IntAct download not implemented - using local file")
    temp_file <- system.file("extdata", "intact.txt", package = "MPAGE")
    if (!file.exists(temp_file)) {
      stop("No IntAct data file found")
    }
  }

  # Process the data
  processed_data <- preprocess_intact_data(
    input_file = temp_file,
    output_file = output_file,
    species_filter = species_filter,
    min_score = min_score
  )

  return(output_file)
}
