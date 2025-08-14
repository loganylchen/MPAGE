#' Calculate Intra-pair Relationships
#'
#' Calculate relationships within modules for gene pairs
#'
#' @param expression_data Matrix with genes as rows and samples as columns
#' @param modules List of modules
#' @param metrics Character vector of metrics to calculate ("ratio", "difference", "correlation")
#'
#' @return List containing intra-module gene pair relationships
#'
#' @export
#'
#' @examples
#' intra_pair_scores <- calculate_intra_pair_relationships(
#'   expression_data = expression_data,
#'   modules = classified_modules,
#'   metrics = c("ratio", "difference", "correlation")
#' )
calculate_intra_pair_relationships <- function(expression_data, modules, 
                                         metrics = c("ratio", "difference", "correlation")) {
  
  if (!is.matrix(expression_data) && !is.data.frame(expression_data)) {
    stop("expression_data must be a matrix or data frame")
  }
  
  if (!is.list(modules)) {
    stop("modules must be a list")
  }
  
  # Convert to matrix
  expr_mat <- as.matrix(expression_data)
  
  # Initialize results
  intra_results <- list()
  
  for (module in modules) {
    if (!"genes" %in% names(module)) next
    
    module_genes <- intersect(module$genes, rownames(expr_mat))
    
    if (length(module_genes) < 2) {
      warning(paste("Module", module$module_id, "has fewer than 2 genes in expression data"))
      next
    }
    
    # Generate all gene pairs within the module
    gene_pairs <- t(combn(module_genes, 2))
    
    if (nrow(gene_pairs) == 0) next
    
    # Calculate relationships for each pair
    pair_scores <- apply(gene_pairs, 1, function(pair) {
      gene1 <- pair[1]
      gene2 <- pair[2]
      
      scores <- list()
      
      if ("ratio" %in% metrics) {
        scores$ratio <- expr_mat[gene1, ] / expr_mat[gene2, ]
      }
      
      if ("difference" %in% metrics) {
        scores$difference <- expr_mat[gene1, ] - expr_mat[gene2, ]
      }
      
      if ("correlation" %in% metrics) {
        scores$correlation <- cor(expr_mat[gene1, ], expr_mat[gene2, ], use = "pairwise.complete.obs")
      }
      
      return(scores)
    })
    
    # Organize results
    module_results <- list(
      module_id = module$module_id,
      genes = module_genes,
      gene_pairs = gene_pairs,
      relationships = pair_scores,
      metrics = metrics
    )
    
    intra_results[[module$module_id]] <- module_results
  }
  
  return(intra_results)
}

#' Calculate Inter-pair Relationships
#'
#' Calculate relationships between modules for gene pairs
#'
#' @param expression_data Matrix with genes as rows and samples as columns
#' @param modules List of modules
#' @param module_pairs List of module pairs to analyze
#' @param metrics Character vector of metrics to calculate ("ratio", "difference")
#'
#' @return List containing inter-module gene pair relationships
#'
#' @export
#'
#' @examples
#' inter_pair_scores <- calculate_inter_pair_relationships(
#'   expression_data = expression_data,
#'   modules = classified_modules,
#'   module_pairs = list(c("m6A_writer", "m6A_reader")),
#'   metrics = c("ratio", "difference")
#' )
calculate_inter_pair_relationships <- function(expression_data, modules, 
                                         module_pairs = NULL,
                                         metrics = c("ratio", "difference")) {
  
  if (!is.matrix(expression_data) && !is.data.frame(expression_data)) {
    stop("expression_data must be a matrix or data frame")
  }
  
  if (!is.list(modules)) {
    stop("modules must be a list")
  }
  
  # Convert to matrix
  expr_mat <- as.matrix(expression_data)
  
  # Create module mapping
  module_map <- setNames(lapply(modules, function(m) m$genes), 
                         sapply(modules, function(m) m$module_id))
  
  # If no module pairs specified, create all possible pairs
  if (is.null(module_pairs)) {
    module_ids <- names(module_map)
    module_pairs <- combn(module_ids, 2, simplify = FALSE)
  }
  
  # Initialize results
  inter_results <- list()
  
  for (pair in module_pairs) {
    if (length(pair) != 2) {
      warning("Each module pair must contain exactly 2 module IDs")
      next
    }
    
    module1_id <- pair[1]
    module2_id <- pair[2]
    
    if (!module1_id %in% names(module_map) || !module2_id %in% names(module_map)) {
      warning(paste("Module pair", paste(pair, collapse = "-"), "not found"))
      next
    }
    
    # Get genes for each module
    genes1 <- intersect(module_map[[module1_id]], rownames(expr_mat))
    genes2 <- intersect(module_map[[module2_id]], rownames(expr_mat))
    
    if (length(genes1) == 0 || length(genes2) == 0) {
      warning(paste("No overlapping genes found for module pair", paste(pair, collapse = "-")))
      next
    }
    
    # Generate all gene pairs between modules
    gene_pairs <- expand.grid(genes1, genes2)
    colnames(gene_pairs) <- c("gene1", "gene2")
    
    if (nrow(gene_pairs) == 0) next
    
    # Calculate relationships for each pair
    pair_scores <- apply(gene_pairs, 1, function(pair) {
      gene1 <- pair["gene1"]
      gene2 <- pair["gene2"]
      
      scores <- list()
      
      if ("ratio" %in% metrics) {
        scores$ratio <- expr_mat[gene1, ] / expr_mat[gene2, ]
      }
      
      if ("difference" %in% metrics) {
        scores$difference <- expr_mat[gene1, ] - expr_mat[gene2, ]
      }
      
      return(scores)
    })
    
    # Organize results
    pair_results <- list(
      module_pair = paste(module1_id, module2_id, sep = "_"),
      module1_id = module1_id,
      module2_id = module2_id,
      genes1 = genes1,
      genes2 = genes2,
      gene_pairs = gene_pairs,
      relationships = pair_scores,
      metrics = metrics
    )
    
    inter_results[[paste(module1_id, module2_id, sep = "_")]] <- pair_results
  }
  
  return(inter_results)
}

#' Calculate Synergy Scores
#'
#' Calculate synergy scores between module pairs
#'
#' @param intra_pair_scores Output from calculate_intra_pair_relationships()
#' @param inter_pair_scores Output from calculate_inter_pair_relationships()
#'
#' @return List containing synergy scores for module pairs
#'
#' @export
#'
#' @examples
#' synergy_scores <- calculate_synergy_scores(
#'   intra_pair_scores,
#'   inter_pair_scores
#' )
calculate_synergy_scores <- function(intra_pair_scores, inter_pair_scores) {
  
  synergy_results <- list()
  
  # Calculate synergy for each module pair
  for (pair_name in names(inter_pair_scores)) {
    pair_data <- inter_pair_scores[[pair_name]]
    
    # Get corresponding intra-module scores
    module1_id <- pair_data$module1_id
    module2_id <- pair_data$module2_id
    
    if (module1_id %in% names(intra_pair_scores) && 
        module2_id %in% names(intra_pair_scores)) {
      
      intra1 <- intra_pair_scores[[module1_id]]
      intra2 <- intra_pair_scores[[module2_id]]
      
      # Calculate synergy metrics
      synergy_metrics <- list()
      
      # Example: Compare inter-module vs intra-module relationships
      for (metric in pair_data$metrics) {
        if (metric %in% c("ratio", "difference")) {
          # Aggregate inter-module relationships
          inter_values <- sapply(pair_data$relationships, function(x) mean(x[[metric]], na.rm = TRUE))
          
          # Aggregate intra-module relationships
          intra1_values <- sapply(intra1$relationships, function(x) mean(x[[metric]], na.rm = TRUE))
          intra2_values <- sapply(intra2$relationships, function(x) mean(x[[metric]], na.rm = TRUE))
          
          # Calculate synergy score (example: relative difference)
          synergy_score <- mean(inter_values, na.rm = TRUE) - 
                          mean(c(intra1_values, intra2_values), na.rm = TRUE)
          
          synergy_metrics[[paste0(metric, "_synergy")]] <- synergy_score
        }
      }
      
      synergy_results[[pair_name]] <- list(
        module_pair = pair_name,
        module1_id = module1_id,
        module2_id = module2_id,
        synergy_scores = synergy_metrics,
        num_inter_pairs = nrow(pair_data$gene_pairs),
        num_intra1_pairs = nrow(intra1$gene_pairs),
        num_intra2_pairs = nrow(intra2$gene_pairs)
      )
    }
  }
  
  return(synergy_results)
}

#' Plot Pair Relationships
#'
#' Visualize gene pair relationships and synergy scores
#'
#' @param synergy_scores Output from calculate_synergy_scores()
#' @param output_dir Directory to save plots
#'
#' @export
#'
#' @examples
#' plot_pair_relationships(synergy_scores, output_dir = "pair_relationship_plots/")
plot_pair_relationships <- function(synergy_scores, output_dir = "pair_relationship_plots") {
  if (!dir.exists(output_dir)) {
    dir.create(output_dir, recursive = TRUE)
  }
  
  # Create synergy heatmap
  if (length(synergy_scores) > 0) {
    synergy_df <- do.call(rbind, lapply(synergy_scores, function(x) {
      data.frame(
        module_pair = x$module_pair,
        module1 = x$module1_id,
        module2 = x$module2_id,
        synergy_score = x$synergy_scores[[1]] %||% 0
      )
    }))
    
    if (nrow(synergy_df) > 0) {
      # Extract module names for visualization
      synergy_df$module1 <- sub("_.*", "", synergy_df$module_pair)
      synergy_df$module2 <- sub(".*_", "", synergy_df$module_pair)
      
      # Create heatmap
      library(reshape2)
      heatmap_data <- reshape2::acast(synergy_df, module1 ~ module2, value.var = "synergy_score")
      
      # Save heatmap
      png(file.path(output_dir, "synergy_heatmap.png"), width = 800, height = 600)
      if (!is.null(heatmap_data) && nrow(heatmap_data) > 0) {
        pheatmap::pheatmap(
          heatmap_data,
          main = "Module Pair Synergy Scores",
          color = colorRampPalette(c("blue", "white", "red"))(50),
          cluster_rows = TRUE,
          cluster_cols = TRUE
        )
      }
      dev.off()
    }
  }
  
  # Create barplot of top synergistic pairs
  synergy_df <- do.call(rbind, lapply(synergy_scores, function(x) {
    data.frame(
      module_pair = x$module_pair,
      synergy_score = x$synergy_scores[[1]] %||% 0
    )
  }))
  
  if (nrow(synergy_df) > 0) {
    top_pairs <- head(synergy_df[order(abs(synergy_df$synergy_score), decreasing = TRUE), ], 10)
    
    p <- ggplot2::ggplot(top_pairs, ggplot2::aes(x = reorder(module_pair, synergy_score), y = synergy_score)) +
      ggplot2::geom_bar(stat = "identity", fill = "steelblue") +
      ggplot2::labs(title = "Top Synergistic Module Pairs",
                    x = "Module Pair", y = "Synergy Score") +
      ggplot2::theme_minimal() +
      ggplot2::theme(axis.text.x = ggplot2::element_text(angle = 45, hjust = 1))
    
    ggplot2::ggsave(file.path(output_dir, "top_synergistic_pairs.png"), 
                   p, width = 10, height = 6)
  }
  
  message("Pair relationship plots saved to: ", output_dir)
}

#' Helper function for NULL-coalescing
#' @noRd
`%||%` <- function(x, y) {
  if (is.null(x)) y else x
}