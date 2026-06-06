#' Build PPI Network (Revised)
#'
#' Construct protein-protein interaction network from multiple data sources using unified interface
#'
#' @param data_sources Character vector of data sources to use ("STRING", "BIOGRID", "INTACT")
#' @param filters List of filtering criteria for each data source. Named list with elements
#'   'string' (numeric score threshold, default: 700), 'intact' (numeric score threshold, default: 0.7),
#'   'biogrid' (character vector of evidence types, default: c("physical"))
#' @param versions List of database versions to use. Named list with elements 'string', 'intact',
#'   'biogrid' specifying version strings for each database
#' @param processed_dir Character string specifying directory path for processed data files (default: "./")
#' @param species Character string specifying species taxonomy ID (default: "9606")
#'
#' @return A named list of igraph objects, one per data source (e.g., \code{string}, \code{biogrid}, \code{intact}).
#'   Use \code{\link{merge_ppi_networks}} to combine them into a single network.
#'
#' @export
#'
#' @examples
#' # Build network using RNA modification proteins
#' ppi_network <- build_ppi_network()
#'
#' # Build network with specific proteins
#' rna_mod_proteins <- get_rna_mod_proteins()
#' ppi_network <- build_ppi_network()
build_ppi_network <- function(data_sources = c("STRING", "BIOGRID", "INTACT"),
                              filters = list(
                                string = 700,
                                intact = 0.7,
                                biogrid = c("physical")
                              ),
                              versions = list(
                                string = "12",
                                intact = "2025-03-28",
                                biogrid = "4.4.248"
                              ),
                              processed_dir = "./",
                              species = "9606") {
  # Validate inputs
  data_sources <- toupper(data_sources)
  valid_sources <- c("STRING", "BIOGRID", "INTACT")
  invalid_sources <- setdiff(data_sources, valid_sources)
  if (length(invalid_sources) > 0) {
    stop(
      "Invalid data sources: ", paste(invalid_sources, collapse = ", "),
      ". Valid sources: ", paste(valid_sources, collapse = ", ")
    )
  }

  # Ensure processed directory exists
  if (!dir.exists(processed_dir)) {
    dir.create(processed_dir, recursive = TRUE)
  }


  # Initialize networks list
  networks <- list()
  processed_files <- list()

  # Process each data source
  for (source in data_sources) {
    message(sprintf("Processing %s database...", source))

    # Check for processed file
    processed_file_name <- sprintf("%s_%s_v%s.rds", tolower(source), as.character(species), gsub("-", "", versions[[tolower(source)]]))
    processed_file <- system.file("extdata", processed_file_name, package = "mpage")



    if (!file.exists(processed_file)) {
      stop(sprintf("Processing %s data not exists: %s", source, processed_file_name))
    }
    processed_df <- readRDS(processed_file)
    message(sprintf("%s Network Summary:", source))
    message(sprintf("Processed %d unique interactions", nrow(processed_df)))
    message(sprintf("Unique proteins: %d", length(unique(c(processed_df$gene1, processed_df$gene2)))))
    if (is.numeric(filters[[tolower(source)]])) {
      processed_df <- processed_df %>% dplyr::filter(score > filters[[tolower(source)]])
    } else {
      processed_df <- processed_df %>% dplyr::filter(est %in% filters[[tolower(source)]])
    }
    network <- igraph::graph_from_data_frame(
      processed_df[, c("gene1", "gene2")],
      directed = FALSE
    )

    # Add edge attributes
    igraph::E(network)$source <- source

    message(sprintf("Nodes: %d", length(igraph::V(network))))
    message(sprintf("Edges: %d", length(igraph::E(network))))
    message(sprintf("Average degree: %.2f", mean(igraph::degree(network))))
    networks[[tolower(source)]] <- network
  }


  # Check if any networks were loaded
  if (length(networks) == 0) {
    warning("No valid networks found from any data sources")
    return(igraph::make_empty_graph())
  }





  return(networks)
}

#' Merge PPI Networks
#'
#' Merge protein-protein interaction networks from multiple data sources
#'
#' @param networks List of igraph objects to merge (typically from build_ppi_network())
#' @param merge_method Method to combine networks ("union" or "intersection"). "union" includes all interactions from any network, "intersection" includes only interactions present in ALL networks.
#' @param add_source_labels Whether to add source labels to edges indicating which databases contributed each interaction
#'
#' @return igraph object representing the merged network
#'
#' @export
#'
#' @examples
#' \dontrun{
#' # Build individual networks
#' networks <- build_ppi_network(data_sources = c("STRING", "BIOGRID", "INTACT"))
#'
#' # Merge all networks using union (all interactions)
#' merged_network <- merge_ppi_networks(networks, merge_method = "union")
#'
#' # Merge using intersection (only interactions present in all networks)
#' intersected_network <- merge_ppi_networks(networks, merge_method = "intersection")
#'
#' # Visualize merged network
#' plot(merged_network, vertex.size = 5, edge.width = 1)
#' }
merge_ppi_networks <- function(networks, merge_method = "union", add_source_labels = TRUE) {
  if (length(networks) == 0) {
    warning("No networks provided for merging")
    return(igraph::make_empty_graph())
  }

  if (length(networks) == 1) {
    message("Only one network provided, returning it directly")
    return(networks[[1]])
  }

  message(sprintf("Merging %d networks using method: %s", length(networks), merge_method))

  # Convert all networks to data frames for merging
  network_data <- list()

  for (source_name in names(networks)) {
    net <- networks[[source_name]]
    if (!is.null(net) && igraph::vcount(net) > 0) {
      edges <- igraph::as_data_frame(net, what = "edges")

      # Ensure consistent column names and formatting
      if (nrow(edges) > 0) {
        edges$source <- source_name
        edges$from <- toupper(as.character(edges$from))
        edges$to <- toupper(as.character(edges$to))

        # Handle missing weights
        if (!"weight" %in% colnames(edges) || all(is.na(edges$weight))) {
          edges$weight <- 1
        }

        # Ensure consistent ordering of gene pairs to avoid duplicates
        edges <- edges %>%
          dplyr::mutate(
            geneA = pmin(from, to),
            geneB = pmax(from, to)
          ) %>%
          dplyr::select(geneA, geneB, weight, source)

        network_data[[source_name]] <- edges
      }
    }
  }

  if (length(network_data) == 0) {
    warning("No valid network data found for merging")
    return(igraph::make_empty_graph())
  }

  # Combine all edge data from different sources
  combined_edges <- dplyr::bind_rows(network_data)

  if (nrow(combined_edges) == 0) {
    warning("No edges found in combined network data")
    return(igraph::make_empty_graph())
  }

  # Handle merge based on method
  if (merge_method == "union") {
    # Union: include all interactions from any network
    merged_edges <- combined_edges %>%
      dplyr::group_by(geneA, geneB) %>%
      dplyr::summarise(
        weight = mean(weight, na.rm = TRUE),
        sources = paste(unique(source), collapse = ";"),
        .groups = "drop"
      )
  } else if (merge_method == "intersection") {
    # Intersection: only include interactions present in ALL networks
    source_names <- names(network_data)

    # Count how many sources each interaction appears in
    edge_counts <- combined_edges %>%
      dplyr::group_by(geneA, geneB) %>%
      dplyr::summarise(
        source_count = length(unique(source)),
        weight = mean(weight, na.rm = TRUE),
        sources = paste(unique(source), collapse = ";"),
        .groups = "drop"
      )

    # Keep only interactions present in all networks
    merged_edges <- edge_counts %>%
      dplyr::filter(source_count == length(source_names))
  } else {
    stop("Invalid merge_method. Use 'union' or 'intersection'")
  }

  # Create merged network
  merged_network <- igraph::graph_from_data_frame(
    merged_edges[, c("geneA", "geneB")],
    directed = FALSE
  )

  # Add comprehensive edge attributes
  igraph::E(merged_network)$weight <- merged_edges$weight
  if (add_source_labels) {
    igraph::E(merged_network)$sources <- merged_edges$sources
  }

  # Add vertex attributes
  igraph::V(merged_network)$degree <- igraph::degree(merged_network)

  # Summary statistics
  message("=== Merged Network Summary ===")
  message(sprintf("Data sources: %s", paste(names(network_data), collapse = ", ")))
  message(sprintf("Merge method: %s", merge_method))
  message(sprintf("Total nodes: %d", igraph::vcount(merged_network)))
  message(sprintf("Total edges: %d", igraph::ecount(merged_network)))
  message(sprintf("Network density: %.4f", igraph::edge_density(merged_network)))
  message(sprintf("Average degree: %.2f", mean(igraph::degree(merged_network))))

  # Detailed source contribution analysis
  if (add_source_labels && "sources" %in% colnames(merged_edges)) {
    source_list <- strsplit(merged_edges$sources, ";")
    source_summary <- table(unlist(source_list))
    message("=== Source Contribution ===")
    print(source_summary)

    # Multi-source interactions
    multi_source <- lengths(source_list) > 1
    message(sprintf(
      "Interactions from multiple sources: %d (%.1f%%)",
      sum(multi_source), 100 * sum(multi_source) / nrow(merged_edges)
    ))
  }

  return(merged_network)
}
