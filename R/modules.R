set.seed(422)
#' Iterative Module Identification with Target Proteins
#'
#' Identify network modules using louvain and fast_greedy algorithms iteratively,
#' focusing on modules that contain target proteins
#'
#' @param ppi_network An igraph object representing the PPI network
#' @param min_module_size Minimum number of nodes in a module
#' @param max_module_size Maximum number of nodes in a module
#' @param target_proteins Character vector of target proteins to focus on
#'
#' @return A list containing all identified modules that contain target proteins
#'
#' @export
#'
#' @examples
#' modules <- identify_modules_iterative(
#'   ppi_network = ppi_network,
#'   min_module_size = 10,
#'   max_module_size = 1000,
#'   target_proteins = c("DNMT1", "DNMT3A", "TET1")
#' )
identify_modules_iterative <- function(ppi_network,
                                       min_module_size = 10,
                                       max_module_size = 1000,
                                       target_proteins = NULL) {
  if (!igraph::is.igraph(ppi_network)) {
    stop("ppi_network must be an igraph object")
  }

  if (is.null(target_proteins) || length(target_proteins) == 0) {
    stop("target_proteins must be provided")
  }

  # Get all vertex names in the network
  all_vertices <- igraph::V(ppi_network)$name

  # Check if target proteins are in the network
  target_in_network <- target_proteins[target_proteins %in% all_vertices]
  if (length(target_in_network) == 0) {
    warning("None of the target proteins are found in the network")
    return(list())
  }

  # Initialize result list
  final_modules <- list()

  # Helper function to check if a module contains target proteins
  contains_targets <- function(module_genes) {
    any(module_genes %in% target_in_network)
  }

  # Helper function for recursive clustering
  recursive_clustering <- function(subgraph, targets_in_subgraph, depth = 1) {
    modules <- list()

    # Use fast_greedy on the subgraph
    tryCatch(
      {
        set.seed(422)
        fc <- igraph::cluster_fast_greedy(subgraph)
        module_list <- igraph::groups(fc)

        # Process each module
        for (i in seq_along(module_list)) {
          module_genes <- module_list[[i]]
          module_size <- length(module_genes)

          # Check if module contains target proteins
          if (contains_targets(module_genes)) {
            # Check size criteria
            if (module_size >= min_module_size && module_size <= max_module_size) {
              # Valid module - add to results
              modules[[length(modules) + 1]] <- list(
                module_id = paste0("RC", depth, "_", length(modules) + 1),
                genes = module_genes,
                size = module_size,
                algorithm = "FASTGREEDY",
                depth = depth,
                target_proteins = module_genes[module_genes %in% targets_in_subgraph]
              )
            } else if (module_size > max_module_size) {
              # Module too large - recurse further
              subgraph_module <- igraph::induced_subgraph(subgraph, module_genes)

              # Only recurse if we have targets in this subgraph
              targets_here <- targets_in_subgraph[targets_in_subgraph %in% module_genes]
              if (length(targets_here) > 0) {
                sub_modules <- recursive_clustering(subgraph_module, targets_here, depth + 1)
                modules <- c(modules, sub_modules)
              }
            }
            # Skip modules smaller than min_module_size
          }
        }
      },
      error = function(e) {
        warning("Fast Greedy clustering failed at depth ", depth, ": ", e$message)
      }
    )

    return(modules)
  }

  # Step 1: Initial clustering with Louvain
  tryCatch(
    {
      set.seed(422)
      louvain_cluster <- igraph::cluster_louvain(ppi_network)
      louvain_modules <- igraph::groups(louvain_cluster)

      # Process each louvain module
      candidate_modules <- list()

      for (i in seq_along(louvain_modules)) {
        module_genes <- louvain_modules[[i]]
        module_size <- length(module_genes)

        # Check if module contains target proteins
        if (contains_targets(module_genes)) {
          if (module_size >= min_module_size && module_size <= max_module_size) {
            # Valid module - keep as is
            candidate_modules[[length(candidate_modules) + 1]] <- list(
              module_id = paste0("L", i),
              genes = module_genes,
              size = module_size,
              algorithm = "LOUVAIN",
              depth = 0,
              target_proteins = module_genes[module_genes %in% target_in_network]
            )
          } else if (module_size > max_module_size) {
            # Too large - apply fast_greedy recursively
            subgraph <- igraph::induced_subgraph(ppi_network, module_genes)
            targets_here <- target_in_network[target_in_network %in% module_genes]

            if (length(targets_here) > 0) {
              sub_modules <- recursive_clustering(subgraph, targets_here, 1)
              candidate_modules <- c(candidate_modules, sub_modules)
            }
          }
          # Skip modules smaller than min_module_size
        }
      }
    },
    error = function(e) {
      warning("Louvain clustering failed: ", e$message)
      candidate_modules <- list()
    }
  )

  # Return all modules with target proteins (no deduplication)
  final_modules <- candidate_modules

  # Add summary information
  if (length(final_modules) > 0) {
    message("Identified ", length(final_modules), " modules containing target proteins")
    message(
      "Total target proteins covered: ",
      length(unique(unlist(lapply(final_modules, function(x) x$target_proteins))))
    )
  }

  return(final_modules)
}
