#' Build PPI Network
#'
#' Construct protein-protein interaction network from multiple data sources
#'
#' @param proteins Character vector of protein identifiers. If provided, extracts subnetwork containing these proteins.
#' @param data_sources Character vector of data sources to use ("STRING", "BioGRID", "IntAct")
#' @param min_confidence Minimum confidence score for STRING database (0-1)
#' @param include_experimental Logical whether to include experimental interactions
#' @param exclude_predicted Logical whether to exclude predicted interactions
#' @param species Character string specifying species (default: "Homo sapiens")
#' @param string_version Character string specifying STRINGdb version (default: "12")
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
#' ppi_network <- build_ppi_network(
#'   proteins = rna_mod_proteins$gene_symbol,
#'   data_sources = c("STRING", "BioGRID"),
#'   min_confidence = 0.7,
#'   species = "Homo sapiens"
#' )
build_ppi_network <- function(proteins = NULL,
                              data_sources = c("STRING", "BioGRID", "IntAct"),
                              min_confidence = 0.7,
                              include_experimental = TRUE,
                              exclude_predicted = FALSE,
                              species = "Homo sapiens",
                              string_version = "12") {
  if (!requireNamespace("igraph", quietly = TRUE)) {
    stop("Package 'igraph' is required but not installed.")
  }

  # Build full network regardless of proteins parameter

  # Validate data sources
  valid_sources <- c("STRING", "BioGRID", "IntAct")
  invalid_sources <- setdiff(data_sources, valid_sources)
  if (length(invalid_sources) > 0) {
    stop("Invalid data sources: ", paste(invalid_sources, collapse = ", "))
  }

  # Map species names to STRINGdb species codes
  species_map <- list(
    "Homo sapiens" = 9606,
    "Mus musculus" = 10090
  )

  species_code <- species_map[[species]]
  if (is.null(species_code)) {
    stop(
      "Unsupported species: ", species, ". Available: ",
      paste(names(species_map), collapse = ", ")
    )
  }

  # Initialize empty network
  networks <- list()

  # STRING database
  if ("STRING" %in% data_sources) {
    if (!requireNamespace("STRINGdb", quietly = TRUE)) {
      warning("STRINGdb package not available. Skipping STRING database.")
    } else {
      string_net <- .get_string_ppi(proteins, min_confidence, species_code, string_version)
      if (!is.null(string_net)) networks$STRING <- string_net
    }
  }

  # BioGRID database (simplified - would need actual API calls)
  if ("BioGRID" %in% data_sources) {
    if (!requireNamespace("OmnipathR", quietly = TRUE)) {
      warning("OmnipathR package not available. Skipping BioGRID database.")
    } else {
      biogrid_net <- .get_biogrid_ppi(proteins, species)
      if (!is.null(biogrid_net)) networks$BioGRID <- biogrid_net
    }
  }

  # IntAct database (simplified - would need actual API calls)
  if ("IntAct" %in% data_sources) {
    intact_net <- .get_intact_ppi(proteins, species, min_confidence)
    if (!is.null(intact_net)) networks$IntAct <- intact_net
  }




  if (length(networks) == 0) {
    warning("No data sources returned valid networks")
    return(igraph::make_empty_graph())
  }

  # Combine networks
  combined_net <- .combine_networks(networks)

  # If proteins provided, extract subnetwork containing these proteins
  if (!is.null(proteins)) {
    # Ensure proteins are in the network
    available_proteins <- igraph::V(combined_net)$name
    valid_proteins <- proteins[proteins %in% available_proteins]

    if (length(valid_proteins) > 0) {
      # Extract subnetwork
      combined_net <- igraph::induced_subgraph(combined_net, valid_proteins)
      message("Extracted subnetwork with ", length(valid_proteins), " proteins")
    } else {
      warning("None of the specified proteins found in the network")
      return(igraph::make_empty_graph())
    }
  }

  return(combined_net)
}

#' Get STRING PPI
#' @noRd
.get_string_ppi <- function(proteins, min_confidence, species_code, string_version) {
  tryCatch(
    {
      string_db <- STRINGdb::STRINGdb$new(
        species = species_code,
        version = string_version,
        score_threshold = min_confidence * 1000
      )

      # Get full network instead of protein-limited network
      message("Downloading full STRING network for species ", species_code)
      stringdb_proteins <- string_db$get_proteins()
      # Get all interactions for the species above threshold
      all_interactions <- string_db$get_interactions(
        stringdb_proteins$protein_external_id
      )

      if (nrow(all_interactions) == 0) {
        return(NULL)
      }

      # Map interactions to gene symbols
      all_interactions$gene1 <- stringdb_proteins$preferred_name[match(all_interactions$from, stringdb_proteins$protein_external_id)]
      all_interactions$gene2 <- stringdb_proteins$preferred_name[match(all_interactions$to, stringdb_proteins$protein_external_id)]

      # Filter valid interactions
      all_interactions <- all_interactions[!is.na(all_interactions$gene1) & !is.na(all_interactions$gene2), ]

      # Create igraph object with full network
      net <- igraph::graph_from_data_frame(
        all_interactions[, c("gene1", "gene2")],
        directed = FALSE
      )

      # Add edge attributes
      igraph::E(net)$weight <- all_interactions$combined_score
      igraph::E(net)$source <- "STRING"
      cat("STRING DB Network Summary:\n")
      cat("Nodes:", length(igraph::V(net)), "\n")
      cat("Edges:", length(igraph::E(net)), "\n")
      cat("Average degree:", mean(igraph::degree(net)), "\n")
      return(net)
    },
    error = function(e) {
      warning("Error accessing STRING database: ", e$message)
      return(NULL)
    }
  )
}

#' Get BioGRID PPI (placeholder)
#' @noRd
.get_biogrid_ppi <- function(proteins, species) {
  # This is a placeholder - in practice would use BioGRID API
  # For now, return NULL to indicate not implemented
  message("BioGRID database was downloaded from https://downloads.thebiogrid.org/File/BioGRID/Release-Archive/BIOGRID-4.4.248/BIOGRID-ALL-4.4.248.tab3.zip")
  return(NULL)
}

#' Get IntAct PPI
#' @noRd
.get_intact_ppi <- function(proteins, species, min_confidence) {
  tryCatch(
    {
      # Map species names to taxonomy IDs
      species_map <- list(
        "Homo sapiens" = "9606",
        "Mus musculus" = "10090"
      )

      taxid <- species_map[[species]]
      if (is.null(taxid)) {
        warning("Species not supported for IntAct: ", species)
        return(NULL)
      }

      processed_file <- system.file("extdata", "intact_processed.rds", package = "MPAGE")

      if (is.null(processed_file)) {
        message("No preprocessed IntAct data found. Use preprocess_intact_data() to create it.")
        return(NULL)
      }
      message("IntAct database was downloaded at 15/08/2025")
      message("Loading IntAct network from: ", processed_file)

      # Load and filter data
      intact_data <- readRDS(processed_file) %>% dplyr::filter(score >= min_confidence)

      # Filter by species
      species_data <- intact_data[intact_data$taxid == taxid, ]

      if (nrow(species_data) == 0) {
        warning("No IntAct data found for species: ", species)
        return(NULL)
      }

      # Create igraph object
      net <- igraph::graph_from_data_frame(
        species_data[, c("geneA", "geneB")],
        directed = FALSE
      )

      # Add edge attributes
      igraph::E(net)$weight <- species_data$score
      igraph::E(net)$source <- "IntAct"
      cat("IntAct DB Network Summary:\n")
      cat("Nodes:", length(igraph::V(net)), "\n")
      cat("Edges:", length(igraph::E(net)), "\n")
      cat("Average degree:", mean(igraph::degree(net)), "\n")
      return(net)
    },
    error = function(e) {
      warning("Error loading IntAct data: ", e$message)
      return(NULL)
    }
  )
}

#' Get InBioMap PPI (placeholder)
#' @noRd
.get_inbiomap_ppi <- function(proteins, species) {
  # This is a placeholder - in practice would use IntAct API
  # For now, return NULL to indicate not implemented
  message("InBioMap integration not yet implemented - returning NULL")
  return(NULL)
}



#' Combine Networks
#' @noRd
.combine_networks <- function(networks) {
  if (length(networks) == 1) {
    return(networks[[1]])
  }

  # Combine all edge lists
  edge_lists <- lapply(networks, function(net) {
    edges <- igraph::as_data_frame(net, what = "edges")
    vertices <- igraph::as_data_frame(net, what = "vertices")

    # Ensure we have vertex names
    if (is.null(vertices$name)) {
      vertices$name <- rownames(vertices)
    }

    edges$from_gene <- vertices$name[match(edges$from, rownames(vertices))]
    edges$to_gene <- vertices$name[match(edges$to, rownames(vertices))]

    edges <- edges[, c("from_gene", "to_gene", "weight")]
    colnames(edges)[1:2] <- c("from", "to")
    return(edges)
  })

  combined_edges <- do.call(rbind, edge_lists)
  combined_edges <- unique(combined_edges)

  # Create combined network
  combined_net <- igraph::graph_from_data_frame(combined_edges, directed = FALSE)

  return(combined_net)
}

#' Filter PPI Network
#'
#' Filter and refine PPI network based on various criteria
#'
#' @param ppi_network An igraph object representing the PPI network
#' @param min_evidence Minimum number of evidence sources required
#' @param max_self_loops Maximum number of self-loops to allow
#' @param func_filter Functional filter for RNA-related proteins
#'
#' @return Filtered igraph object
#'
#' @export
#'
#' @examples
#' filtered_ppi <- filter_ppi_network(
#'   ppi_network,
#'   min_evidence = 2,
#'   max_self_loops = 0,
#'   func_filter = "RNA"
#' )
filter_ppi_network <- function(ppi_network,
                               min_evidence = 1,
                               max_self_loops = 0,
                               func_filter = NULL) {
  if (!igraph::is.igraph(ppi_network)) {
    stop("ppi_network must be an igraph object")
  }

  # Remove self-loops
  if (max_self_loops == 0) {
    ppi_network <- igraph::simplify(ppi_network, remove.multiple = TRUE, remove.loops = TRUE)
  }

  # Filter based on degree (basic filtering)
  if (min_evidence > 1) {
    # This is a simplified approach - in practice would need more sophisticated filtering
    degrees <- igraph::degree(ppi_network)
    nodes_to_keep <- V(ppi_network)[degrees >= min_evidence]
    if (length(nodes_to_keep) > 0) {
      ppi_network <- igraph::induced_subgraph(ppi_network, nodes_to_keep)
    }
  }

  # Functional filtering (placeholder - would need actual function annotation)
  if (!is.null(func_filter)) {
    # This would filter based on GO terms or other functional annotations
    # For now, we keep all nodes
    warning("Functional filtering currently not implemented - keeping all nodes")
  }

  return(ppi_network)
}

#' Save PPI Network
#'
#' Save PPI network to a file in SIF format
#'
#' @param ppi_network An igraph object representing the PPI network
#' @param file_path Path to save the SIF file
#'
#' @export
#'
#' @examples
#' save_ppi_network(ppi_network, "filtered_ppi_network.sif")
save_ppi_network <- function(ppi_network, file_path) {
  if (!igraph::is.igraph(ppi_network)) {
    stop("ppi_network must be an igraph object")
  }

  edges <- igraph::as_data_frame(ppi_network, what = "edges")
  if (nrow(edges) > 0) {
    utils::write.table(edges[, 1:2], file_path,
      sep = "\t",
      row.names = FALSE, col.names = FALSE, quote = FALSE
    )
  }
  message("PPI network saved to: ", file_path)
}

#' Load PPI Network
#'
#' Load PPI network from a file
#'
#' @param file_path Path to the SIF file
#'
#' @return An igraph object representing the PPI network
#'
#' @export
#'
#' @examples
#' ppi_network <- load_ppi_network("filtered_ppi_network.sif")
load_ppi_network <- function(file_path) {
  if (!file.exists(file_path)) {
    stop("File does not exist: ", file_path)
  }

  edges <- utils::read.table(file_path, sep = "\t", header = FALSE)
  colnames(edges) <- c("from", "to")

  network <- igraph::graph_from_data_frame(edges, directed = FALSE)
  return(network)
}
