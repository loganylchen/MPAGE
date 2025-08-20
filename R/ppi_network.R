#' Build PPI Network (Revised)
#'
#' Construct protein-protein interaction network from multiple data sources using unified interface
#'
#' @param data_sources Character vector of data sources to use ("STRING", "BIOGRID", "INTACT")
#' @param filters
#' @param versions
#' @param processed_dir
#' @param species Character string specifying species taxonomy ID (default: "9606")
#'
#' @return An igraph object representing the PPI network
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
    processed_file <- system.file("extdata", processed_file_name, package = "MPAGE")



    if (!file.exists(processed_file)) {
      stop(sprintf("Processing %s data not exists: %s", source, processed_file_name))
    }
    processed_df <- readRDS(processed_file)
    cat(sprintf("%s Network Summary:\n", source))
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
    igraph::E(network)$source <- "STRING"

    cat("Nodes:", length(igraph::V(network)), "\n")
    cat("Edges:", length(igraph::E(network)), "\n")
    cat("Average degree:", mean(igraph::degree(network)), "\n")
    networks[[tolower(source)]] <- network
  }


  # Check if any networks were loaded
  if (length(networks) == 0) {
    warning("No valid networks found from any data sources")
    return(igraph::make_empty_graph())
  }





  return(networks)
}
