# BioGRID Database Preprocessing Functions

#' Preprocess BioGRID Database
#'
#' Process raw BioGRID TSV file into clean interaction data
#'
#' @param input_file Path to BioGRID TSV file (typically BIOGRID-ALL-*.tab3.txt)
#' @param output_file Path to save processed RDS file (optional)
#' @param species_filter Vector of NCBI taxonomy IDs to filter (default: c("9606", "10090"))
#' @param min_score Minimum interaction confidence score (BioGRID uses 0-1000 scale)
#' @param interaction_types Character vector of interaction types to include
#'   (default: c("physical", "genetic", "direct interaction"))
#' @param experimental_systems Character vector of experimental systems to include
#'   (default: includes major detection methods)
#'
#' @return Cleaned interaction data frame with columns: geneA, geneB, score, taxid,
#'   interaction_type, experimental_system, pubmed_id
#'
#' @export
#'
#' @examples
#' \dontrun{
#' # Process human and mouse BioGRID data
#' biogrid_data <- preprocess_biogrid_data(
#'   input_file = "BIOGRID-ALL-4.4.248.tab3.txt",
#'   species_filter = c("9606", "10090"),
#'   min_score = 0
#' )
#' }
preprocess_biogrid_data <- function(input_file,
                                    output_file = NULL,
                                    species_filter = c("9606", "10090"),
                                    min_score = 0,
                                    interaction_types = c("physical", "genetic", "direct interaction", "association"),
                                    experimental_systems = c(
                                      "Two-hybrid", "Affinity Capture-MS", "Co-immunoprecipitation",
                                      "Pull-down", "Yeast Two-hybrid", "Tandem Affinity Purification"
                                    )) {
  if (!file.exists(input_file)) {
    stop("Input file not found: ", input_file)
  }

  message("Loading BioGRID data from: ", input_file)

  # Read BioGRID data
  tryCatch(
    {
      df <- data.table::fread(input_file, sep = "\t", check.names = TRUE, quote = "")
    },
    error = function(e) {
      stop("Error reading BioGRID file: ", e$message)
    }
  )

  message("Processing ", nrow(df), " interactions...")

  # Validate required columns exist
  required_cols <- c(
    "BioGRID Interaction ID", "Entrez Gene Interactor A", "Entrez Gene Interactor B",
    "BioGRID ID Interactor A", "BioGRID ID Interactor B", "Systematic Name Interactor A",
    "Systematic Name Interactor B", "Official Symbol Interactor A", "Official Symbol Interactor B",
    "Organism Interactor A", "Organism Interactor B", "Interaction Type", "Experimental System",
    "Author", "Pubmed ID", "Throughput", "Score"
  )

  missing_cols <- setdiff(required_cols, colnames(df))
  if (length(missing_cols) > 0) {
    warning("Missing expected columns: ", paste(missing_cols, collapse = ", "))
  }

  # Process the data
  processed_data <- df %>%
    dplyr::filter(!is.na(`Official Symbol Interactor A`), !is.na(`Official Symbol Interactor B`)) %>%
    dplyr::mutate(
      geneA = toupper(`Official Symbol Interactor A`),
      geneB = toupper(`Official Symbol Interactor B`),
      taxidA = as.character(`Organism Interactor A`),
      taxidB = as.character(`Organism Interactor B`),
      interaction_type = `Interaction Type`,
      experimental_system = `Experimental System`,
      score = as.numeric(`Score`),
      pubmed_id = as.character(`Pubmed ID`)
    ) %>%
    # Filter by species
    dplyr::filter(
      !is.na(taxidA), !is.na(taxidB),
      taxidA %in% as.character(species_filter),
      taxidB %in% as.character(species_filter),
      taxidA == taxidB
    ) %>%
    # Filter by interaction type and experimental system
    dplyr::filter(
      interaction_type %in% interaction_types | is.na(interaction_type),
      experimental_system %in% experimental_systems | is.na(experimental_system)
    ) %>%
    # Filter by score
    dplyr::filter(score >= min_score) %>%
    # Add taxid column
    dplyr::mutate(taxid = taxidA) %>%
    # Select and clean
    dplyr::select(geneA, geneB, score, taxid, interaction_type, experimental_system, pubmed_id) %>%
    dplyr::filter(!is.na(geneA), !is.na(geneB), geneA != "", geneB != "", geneA != geneB) %>%
    dplyr::distinct() %>%
    dplyr::arrange(dplyr::desc(score))

  # Create summary statistics
  cat("=== BioGRID Processing Summary ===\n")
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
    "Interaction types: %s\n",
    paste(unique(processed_data$interaction_type), collapse = ", ")
  ))
  cat(sprintf(
    "Experimental systems: %s\n",
    paste(unique(processed_data$experimental_system), collapse = ", ")
  ))
  cat(sprintf(
    "Score range: [%d, %d]\n",
    min(processed_data$score, na.rm = TRUE),
    max(processed_data$score, na.rm = TRUE)
  ))

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

  # Interaction type summary
  type_summary <- processed_data %>%
    dplyr::group_by(interaction_type) %>%
    dplyr::summarise(count = n(), .groups = "drop") %>%
    dplyr::arrange(dplyr::desc(count))

  cat("\n=== Interaction Types ===\n")
  print(type_summary)

  # Save if output file provided
  if (!is.null(output_file)) {
    saveRDS(processed_data, output_file)
    message("Saved processed data to: ", output_file)
  }

  return(processed_data)
}

#' Load BioGRID Network
#'
#' Load preprocessed BioGRID data and create network
#'
#' @param processed_file Path to processed RDS file
#' @param min_score Minimum interaction confidence score
#' @param species Filter by species taxonomy ID
#' @param proteins Optional vector of proteins to extract subnetwork
#' @param interaction_types Filter by interaction types
#' @param experimental_systems Filter by experimental systems
#'
#' @return igraph object representing the PPI network
#'
#' @export
#'
#' @examples
#' \dontrun{
#' # Load full human BioGRID network
#' biogrid_network <- load_biogrid_network(
#'   processed_file = "biogrid_processed.rds",
#'   species = "9606",
#'   min_score = 0
#' )
#' }
load_biogrid_network <- function(processed_file,
                                 min_score = 0,
                                 species = NULL,
                                 proteins = NULL,
                                 interaction_types = NULL,
                                 experimental_systems = NULL) {
  if (!file.exists(processed_file)) {
    stop("Processed file not found: ", processed_file)
  }

  # Load processed data
  biogrid_data <- readRDS(processed_file)

  # Apply filters
  filtered_data <- biogrid_data %>%
    dplyr::filter(score >= min_score)

  if (!is.null(species)) {
    filtered_data <- filtered_data %>%
      dplyr::filter(taxid == as.character(species))
  }

  if (!is.null(interaction_types)) {
    filtered_data <- filtered_data %>%
      dplyr::filter(interaction_type %in% interaction_types)
  }

  if (!is.null(experimental_systems)) {
    filtered_data <- filtered_data %>%
      dplyr::filter(experimental_system %in% experimental_systems)
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
  igraph::E(net)$source <- "BioGRID"
  igraph::E(net)$interaction_type <- filtered_data$interaction_type
  igraph::E(net)$experimental_system <- filtered_data$experimental_system
  igraph::E(net)$pubmed_id <- filtered_data$pubmed_id

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

  cat("=== BioGRID Network Summary ===\n")
  cat(sprintf("Nodes: %d\n", length(igraph::V(net))))
  cat(sprintf("Edges: %d\n", length(igraph::E(net))))
  cat(sprintf("Average degree: %.2f\n", mean(igraph::degree(net))))
  cat(sprintf(
    "Score range: [%d, %d]\n",
    min(igraph::E(net)$weight, na.rm = TRUE),
    max(igraph::E(net)$weight, na.rm = TRUE)
  ))

  return(net)
}

#' Download and Process BioGRID Data
#'
#' Download latest BioGRID data and preprocess it
#'
#' @param output_dir Directory to save processed files
#' @param version BioGRID version to download (e.g., "4.4.248")
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
#' # Download and process latest BioGRID data
#' processed_file <- download_biogrid_data()
#' }
download_biogrid_data <- function(output_dir = tempdir(),
                                  version = "4.4.248",
                                  species_filter = c("9606", "10090"),
                                  min_score = 0,
                                  overwrite = FALSE) {
  # URL for BioGRID data
  biogrid_url <- sprintf("https://downloads.thebiogrid.org/Download/BioGRID/Release-Archive/BIOGRID-%s/BIOGRID-ALL-%s.tab3.zip", version, version)

  # Create output directory
  if (!dir.exists(output_dir)) {
    dir.create(output_dir, recursive = TRUE)
  }

  output_file <- file.path(output_dir, sprintf("BIOGRID-ALL-%s.tab3.zip", version))

  if (file.exists(output_file) & !overwrite) {
    message("Processed file already exists: ", output_file)
    return(output_file)
  }

  message("Downloading BioGRID data...")

  # Download and extract (placeholder - actual implementation would download)
  temp_file <- output_file
  download.file(biogrid_url, temp_file, method = "auto")
  if (!file.exists(temp_file)) {
    warning("BioGRID download not implemented - using local file")
    temp_file <- system.file("extdata", "BIOGRID-ALL.tab3.txt", package = "MPAGE")
    if (!file.exists(temp_file)) {
      stop("No BioGRID data file found")
    }
  }

  # Process the data
  processed_data <- preprocess_biogrid_data(
    input_file = temp_file,
    output_file = output_file,
    species_filter = species_filter,
    min_score = min_score
  )

  return(output_file)
}

#' Get BioGRID Interaction Statistics
#'
#' Get summary statistics for BioGRID interactions
#'
#' @param processed_file Path to processed BioGRID RDS file
#' @param species Species taxonomy ID to filter by
#'
#' @return Data frame with summary statistics
#'
#' @export
get_biogrid_stats <- function(processed_file, species = NULL) {
  if (!file.exists(processed_file)) {
    stop("Processed file not found: ", processed_file)
  }

  data <- readRDS(processed_file)

  if (!is.null(species)) {
    data <- data[data$taxid == as.character(species), ]
  }

  stats <- list(
    total_interactions = nrow(data),
    unique_proteins = length(unique(c(data$geneA, data$geneB))),
    species = unique(data$taxid),
    interaction_types = unique(data$interaction_type),
    experimental_systems = unique(data$experimental_system),
    score_range = range(data$score, na.rm = TRUE),
    publications = length(unique(data$pubmed_id))
  )

  return(stats)
}
