# InBioMap Database Preprocessing Functions

#' Preprocess InBioMap Database
#'
#' Process raw InBioMap TSV file into clean interaction data
#'
#' @param input_file Path to InBioMap TSV file
#' @param output_file Path to save processed RDS file (optional)
#' @param species_filter Vector of NCBI taxonomy IDs to filter (default: c("9606"))
#' @param min_score Minimum interaction confidence score (InBioMap uses 0-1 scale)
#' @param interaction_types Character vector of interaction types to include
#'   (default: c("physical association", "direct interaction"))
#'
#' @return Cleaned interaction data frame with columns: geneA, geneB, score, taxid,
#'   interaction_type, evidence
#'
#' @export
#'
#' @examples
#' \dontrun{
#' # Process human InBioMap data
#' inbiomap_data <- preprocess_inbiomap_data(
#'   input_file = "inbiomap.txt",
#'   species_filter = c("9606"),
#'   min_score = 0.5
#' )
#' }
preprocess_inbiomap_data <- function(input_file,
                                     output_file = NULL,
                                     species_filter = c("9606"),
                                     min_score = 0.0,
                                     interaction_types = c("physical association", "direct interaction", "association")) {
  if (!file.exists(input_file)) {
    stop("Input file not found: ", input_file)
  }

  message("Loading InBioMap data from: ", input_file)

  # Read InBioMap data
  tryCatch(
    {
      df <- data.table::fread(input_file, sep = "\t", header = FALSE, quote = "")
    },
    error = function(e) {
      stop("Error reading InBioMap file: ", e$message)
    }
  )

  message("Processing ", nrow(df), " interactions...")


  processed_data <- df %>%
    dplyr::mutate(
      # Extract from gene name annotations
      geneA_name = stringr::str_extract(
        V5,
        "[A-Za-z0-9-]+\\(gene name\\)"
      ) %>%
        stringr::str_remove("\\(gene name\\)"),
      geneB_name = stringr::str_extract(
        V6,
        "[A-Za-z0-9-]+\\(gene name\\)"
      ) %>% stringr::str_remove("\\(gene name\\)"),
      score = stringr::str_extract(V9, "[0-9\\.]+|") %>% stringr::str_remove("|") %>% as.numeric()
    ) %>%
    dplyr::select(geneA_name, geneB_name, score) %>%
    dplyr::filter(!is.na(score), score >= min_score) %>%
    # Clean and filter
    dplyr::distinct() %>%
    dplyr::arrange(dplyr::desc(score))


  # Create summary statistics
  cat("=== InBioMap Processing Summary ===\n")
  cat(sprintf("Total interactions: %d\n", nrow(processed_data)))
  cat(sprintf(
    "Unique proteins: %d\n",
    length(unique(c(processed_data$geneA, processed_data$geneB)))
  ))
  cat(sprintf(
    "Score range: [%.3f, %.3f]\n",
    min(processed_data$score, na.rm = TRUE),
    max(processed_data$score, na.rm = TRUE)
  ))






  # Save if output file provided
  if (!is.null(output_file)) {
    saveRDS(processed_data, output_file)
    message("Saved processed data to: ", output_file)
  }

  return(processed_data)
}

#' Load InBioMap Network
#'
#' Load preprocessed InBioMap data and create network
#'
#' @param processed_file Path to processed RDS file
#' @param min_score Minimum interaction confidence score
#' @param species Filter by species taxonomy ID
#' @param proteins Optional vector of proteins to extract subnetwork
#' @param interaction_types Filter by interaction types
#'
#' @return igraph object representing the PPI network
#'
#' @export
#'
#' @examples
#' \dontrun{
#' # Load full human InBioMap network
#' inbiomap_network <- load_inbiomap_network(
#'   processed_file = "inbiomap_processed.rds",
#'   species = "9606",
#'   min_score = 0.5
#' )
#' }
load_inbiomap_network <- function(processed_file,
                                  min_score = 0.0,
                                  species = NULL,
                                  proteins = NULL,
                                  interaction_types = NULL) {
  if (!file.exists(processed_file)) {
    stop("Processed file not found: ", processed_file)
  }

  # Load processed data
  inbiomap_data <- readRDS(processed_file)

  # Apply filters
  filtered_data <- inbiomap_data %>%
    dplyr::filter(score >= min_score)

  if (!is.null(species)) {
    filtered_data <- filtered_data %>%
      dplyr::filter(taxid == as.character(species))
  }

  if (!is.null(interaction_types)) {
    filtered_data <- filtered_data %>%
      dplyr::filter(interaction_type %in% interaction_types)
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
  igraph::E(net)$source <- "InBioMap"
  igraph::E(net)$interaction_type <- filtered_data$interaction_type
  igraph::E(net)$evidence <- filtered_data$evidence

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

  cat("=== InBioMap Network Summary ===\n")
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

#' Process Local InBioMap Data
#'
#' Process local InBioMap PSI-MITAB data file
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
#' # Process local InBioMap data
#' processed_file <- process_local_inbiomap_data()
#' }
process_local_inbiomap_data <- function(output_dir = "inst/extdata/",
                                        species_filter = c("9606"),
                                        min_score = 0.0,
                                        overwrite = FALSE) {
  # Use local InBioMap data file
  local_file <- file.path("data-raw", "core.psimitab")

  if (!file.exists(local_file)) {
    stop("InBioMap data file not found at: ", local_file)
  }

  # Create output directory
  if (!dir.exists(output_dir)) {
    dir.create(output_dir, recursive = TRUE)
  }

  output_file <- file.path(output_dir, "inbiomap_processed.rds")

  if (file.exists(output_file) & !overwrite) {
    message("Processed file already exists: ", output_file)
    return(output_file)
  }

  message("Processing local InBioMap data...")

  # Process the data using local file
  processed_data <- preprocess_inbiomap_data(
    input_file = local_file,
    output_file = output_file,
    species_filter = species_filter,
    min_score = min_score
  )

  return(output_file)
}

#' Get InBioMap Interaction Statistics
#'
#' Get summary statistics for InBioMap interactions
#'
#' @param processed_file Path to processed InBioMap RDS file
#' @param species Species taxonomy ID to filter by
#'
#' @return List with summary statistics
#'
#' @export
get_inbiomap_stats <- function(processed_file, species = NULL) {
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
    evidence_types = unique(data$evidence),
    score_range = range(data$score, na.rm = TRUE)
  )

  return(stats)
}

#' Validate InBioMap File Format
#'
#' Check if a file has the expected InBioMap format
#'
#' @param input_file Path to the file to validate
#'
#' @return TRUE if valid, FALSE otherwise
#'
#' @export
validate_inbiomap_file <- function(input_file) {
  if (!file.exists(input_file)) {
    return(FALSE)
  }

  tryCatch(
    {
      # Read first few lines
      preview <- data.table::fread(input_file, sep = "\t", nrows = 5, check.names = TRUE)

      # Check for essential columns
      required_patterns <- c("protein", "gene", "score", "taxid", "interaction", "evidence")

      col_names <- tolower(colnames(preview))
      matches <- sapply(required_patterns, function(pattern) {
        any(grepl(pattern, col_names))
      })

      return(sum(matches) >= 3) # At least 3 key columns found
    },
    error = function(e) {
      return(FALSE)
    }
  )
}

#' Convert InBioMap UniProt IDs to Gene Symbols
#'
#' Helper function to convert UniProt IDs to gene symbols using mapping file
#'
#' @param id_mapping_file Path to UniProt to gene symbol mapping file
#' @param gene_ids Vector of gene IDs to map
#'
#' @return Named vector with gene symbols
#'
#' @export
map_inbiomap_ids_to_symbols <- function(id_mapping_file, gene_ids) {
  if (!file.exists(id_mapping_file)) {
    warning("ID mapping file not found: ", id_mapping_file)
    return(gene_ids)
  }

  tryCatch(
    {
      mapping <- data.table::fread(id_mapping_file, sep = "\t", header = TRUE)

      # Standardize column names
      colnames(mapping) <- tolower(colnames(mapping))

      # Find appropriate columns
      uniprot_col <- which(grepl("uniprot|id", colnames(mapping)))[1]
      gene_col <- which(grepl("gene|symbol", colnames(mapping)))[1]

      if (is.na(uniprot_col) || is.na(gene_col)) {
        warning("Could not find appropriate mapping columns")
        return(gene_ids)
      }

      # Create mapping
      mapping_df <- mapping[, c(uniprot_col, gene_col)]
      colnames(mapping_df) <- c("uniprot_id", "gene_symbol")

      # Remove duplicates, keeping first occurrence
      mapping_df <- mapping_df[!duplicated(mapping_df$uniprot_id), ]

      # Apply mapping
      mapped_symbols <- mapping_df$gene_symbol[match(toupper(gene_ids), toupper(mapping_df$uniprot_id))]

      # Return original IDs where mapping not found
      mapped_symbols[is.na(mapped_symbols)] <- gene_ids[is.na(mapped_symbols)]

      return(mapped_symbols)
    },
    error = function(e) {
      warning("Error in ID mapping: ", e$message)
      return(gene_ids)
    }
  )
}

#' Create Sample InBioMap Data
#'
#' Generate sample InBioMap data for testing purposes
#'
#' @param output_file Path to save sample data
#' @param n_interactions Number of interactions to generate
#' @param species Taxonomy ID for sample data
#'
#' @return Path to generated file
#'
#' @export
create_sample_inbiomap <- function(output_file, n_interactions = 1000, species = "9606") {
  # Sample human proteins
  proteins <- c(
    "TP53", "BRCA1", "EGFR", "MYC", "AKT1", "MAPK1", "JUN", "FOS",
    "RELA", "STAT3", "PIK3CA", "MTOR", "PTEN", "CDK2", "CDKN1A"
  )

  # Generate sample interactions
  sample_data <- data.frame(
    geneA = sample(proteins, n_interactions, replace = TRUE),
    geneB = sample(proteins, n_interactions, replace = TRUE),
    score = runif(n_interactions, 0.1, 1.0),
    taxid = species,
    interaction_type = sample(c("physical association", "direct interaction", "association"),
      n_interactions,
      replace = TRUE, prob = c(0.5, 0.3, 0.2)
    ),
    evidence = sample(c("experimental", "curated", "predicted"),
      n_interactions,
      replace = TRUE, prob = c(0.6, 0.3, 0.1)
    )
  )

  # Remove self-interactions and duplicates
  sample_data <- sample_data[sample_data$geneA != sample_data$geneB, ]
  sample_data <- sample_data[!duplicated(sample_data[, c("geneA", "geneB")]), ]

  # Save sample data
  saveRDS(sample_data, output_file)
  message("Created sample InBioMap data with ", nrow(sample_data), " interactions")

  return(output_file)
}

#' InBioMap Data Quality Report
#'
#' Generate comprehensive quality report for InBioMap data
#'
#' @param processed_file Path to processed InBioMap RDS file
#' @param species Species taxonomy ID to filter by
#'
#' @return List with detailed quality metrics
#'
#' @export
generate_inbiomap_report <- function(processed_file, species = NULL) {
  if (!file.exists(processed_file)) {
    stop("Processed file not found: ", processed_file)
  }

  data <- readRDS(processed_file)

  if (!is.null(species)) {
    data <- data[data$taxid == as.character(species), ]
  }

  # Basic statistics
  basic_stats <- list(
    total_interactions = nrow(data),
    unique_proteins = length(unique(c(data$geneA, data$geneB))),
    species = unique(data$taxid),
    score_range = range(data$score, na.rm = TRUE),
    interaction_types = table(data$interaction_type),
    evidence_types = table(data$evidence)
  )

  # Network properties
  if (nrow(data) > 0) {
    net <- igraph::graph_from_data_frame(data[, c("geneA", "geneB")], directed = FALSE)

    network_stats <- list(
      nodes = length(igraph::V(net)),
      edges = length(igraph::E(net)),
      density = igraph::edge_density(net),
      avg_degree = mean(igraph::degree(net)),
      diameter = igraph::diameter(net, directed = FALSE),
      connected_components = igraph::components(net)$no
    )
  } else {
    network_stats <- list(nodes = 0, edges = 0, density = 0, avg_degree = 0, diameter = 0, connected_components = 0)
  }

  # Quality metrics
  quality_metrics <- list(
    missing_values = sum(is.na(data)),
    duplicate_interactions = nrow(data) - nrow(unique(data[, c("geneA", "geneB")])),
    self_interactions = sum(data$geneA == data$geneB, na.rm = TRUE),
    avg_score = mean(data$score, na.rm = TRUE),
    median_score = median(data$score, na.rm = TRUE)
  )

  report <- list(
    basic_stats = basic_stats,
    network_stats = network_stats,
    quality_metrics = quality_metrics,
    timestamp = Sys.time()
  )

  return(report)
}

#' Export InBioMap Data for Visualization
#'
#' Export processed InBioMap data in formats suitable for visualization
#'
#' @param processed_file Path to processed InBioMap RDS file
#' @param output_file Path for output file (CSV format)
#' @param species Species taxonomy ID to filter by
#' @param min_score Minimum score threshold
#'
#' @return Path to exported file
#'
#' @export
export_inbiomap_for_visualization <- function(processed_file, output_file, species = NULL, min_score = 0.0) {
  if (!file.exists(processed_file)) {
    stop("Processed file not found: ", processed_file)
  }

  data <- readRDS(processed_file)

  # Apply filters
  if (!is.null(species)) {
    data <- data[data$taxid == as.character(species), ]
  }

  data <- data[data$score >= min_score, ]

  # Prepare for export
  export_data <- data %>%
    dplyr::mutate(
      score_normalized = (score - min(score)) / (max(score) - min(score)),
      interaction_id = paste(geneA, geneB, sep = "_")
    ) %>%
    dplyr::select(
      interaction_id, geneA, geneB, score, score_normalized,
      interaction_type, evidence, taxid
    )

  # Write to CSV
  write.csv(export_data, output_file, row.names = FALSE)
  message("Exported ", nrow(export_data), " interactions to ", output_file)

  return(output_file)
}

#' Merge InBioMap with Other PPI Databases
#'
#' Merge InBioMap data with other PPI databases for comprehensive analysis
#'
#' @param inbiomap_file Path to processed InBioMap RDS file
#' @param other_files Named list of other PPI database files (e.g., STRING, BioGRID)
#' @param output_file Path for merged output file
#' @param merge_method How to handle duplicate interactions ("max_score", "mean_score", "keep_all")
#'
#' @return Path to merged file
#'
#' @export
merge_ppi_databases <- function(inbiomap_file, other_files, output_file, merge_method = "max_score") {
  # Load InBioMap data
  inbiomap_data <- readRDS(inbiomap_file)
  inbiomap_data$source <- "InBioMap"

  # Load other databases
  all_data <- list(inbiomap_data)

  for (source_name in names(other_files)) {
    if (file.exists(other_files[[source_name]])) {
      data <- readRDS(other_files[[source_name]])
      data$source <- source_name
      all_data[[source_name]] <- data
    }
  }

  # Combine all data
  combined_data <- dplyr::bind_rows(all_data)

  # Handle duplicates based on merge method
  if (merge_method == "max_score") {
    combined_data <- combined_data %>%
      dplyr::group_by(geneA, geneB) %>%
      dplyr::slice_max(score, n = 1) %>%
      dplyr::ungroup()
  } else if (merge_method == "mean_score") {
    combined_data <- combined_data %>%
      dplyr::group_by(geneA, geneB) %>%
      dplyr::summarise(
        score = mean(score, na.rm = TRUE),
        taxid = first(taxid),
        interaction_type = first(interaction_type),
        evidence = first(evidence),
        source = paste(unique(source), collapse = ";"),
        .groups = "drop"
      )
  } else if (merge_method == "keep_all") {
    # Keep all duplicates, add interaction ID
    combined_data$interaction_id <- paste(combined_data$geneA, combined_data$geneB, combined_data$source, sep = "_")
  }

  # Save merged data
  saveRDS(combined_data, output_file)
  message("Merged ", nrow(combined_data), " unique interactions from ", length(all_data), " sources")

  return(output_file)
}

# Helper function for safe symbol usage
`%||%` <- function(x, y) if (is.null(x)) y else x

# Use dplyr::sym for non-standard evaluation
sym <- dplyr::sym

#' @importFrom rlang !!
#' @importFrom rlang sym
NULL
