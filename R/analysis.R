#' Classify Modules
#'
#' Classify modules by modification type or functional role
#'
#' @param modules List of annotated modules
#' @param classification Type of classification ("modification_type" or "functional_role")
#' @param rna_proteins Data frame containing RNA modification proteins
#'
#' @return Classified modules with additional metadata
#'
#' @export
#'
#' @examples
#' classified_modules <- classify_modules(
#'   annotated_modules,
#'   classification = "modification_type",
#'   rna_proteins = rna_mod_proteins
#' )
classify_modules <- function(modules, classification = "modification_type", rna_proteins = NULL) {
  if (!is.list(modules)) {
    stop("modules must be a list")
  }
  
  classified_modules <- modules
  
  for (i in seq_along(modules)) {
    module <- modules[[i]]
    
    if (classification == "modification_type" && !is.null(rna_proteins)) {
      # Classify by dominant modification type
      rna_in_module <- module$genes[module$genes %in% rna_proteins$gene_symbol]
      
      if (length(rna_in_module) > 0) {
        # Get modification types for RNA proteins in module
        mod_types <- rna_proteins$modification_type[rna_proteins$gene_symbol %in% rna_in_module]
        if (length(mod_types) > 0) {
          dominant_type <- names(sort(table(mod_types), decreasing = TRUE))[1]
          module$classification <- dominant_type
        }
      }
    }
    
    if (classification == "functional_role" && !is.null(rna_proteins)) {
      # Classify by dominant functional role
      rna_in_module <- module$genes[module$genes %in% rna_proteins$gene_symbol]
      
      if (length(rna_in_module) > 0) {
        # Get functional roles for RNA proteins in module
        roles <- rna_proteins$functional_role[rna_proteins$gene_symbol %in% rna_in_module]
        if (length(roles) > 0) {
          dominant_role <- names(sort(table(roles), decreasing = TRUE))[1]
          module$classification <- dominant_role
        }
      }
    }
    
    classified_modules[[i]] <- module
  }
  
  return(classified_modules)
}

#' Analyze Module Features
#'
#' Analyze topological and functional features of modules
#'
#' @param modules List of classified modules
#' @param ppi_network Original PPI network (igraph object)
#'
#' @return Data frame with module features
#'
#' @export
#'
#' @examples
#' module_analysis <- analyze_module_features(
#'   classified_modules,
#'   ppi_network = ppi_network
#' )
analyze_module_features <- function(modules, ppi_network) {
  if (!is.list(modules)) {
    stop("modules must be a list")
  }
  
  features <- lapply(seq_along(modules), function(i) {
    module <- modules[[i]]
    
    # Calculate features
    subgraph <- igraph::induced_subgraph(ppi_network, module$genes)
    
    list(
      module_id = module$module_id,
      size = length(module$genes),
      density = igraph::edge_density(subgraph),
      avg_degree = mean(igraph::degree(subgraph)),
      max_degree = max(igraph::degree(subgraph)),
      diameter = ifelse(length(module$genes) > 1, 
                        igraph::diameter(subgraph), 0),
      clustering_coefficient = ifelse(length(module$genes) > 2,
                                      mean(igraph::transitivity(subgraph, type = "local"), na.rm = TRUE), 0),
      betweenness_centrality = ifelse(length(module$genes) > 1,
                                      mean(igraph::betweenness(subgraph), na.rm = TRUE), 0),
      rna_count = ifelse("rna_count" %in% names(module), module$rna_count, 0),
      rna_ratio = ifelse("rna_ratio" %in% names(module), module$rna_ratio, 0),
      classification = ifelse("classification" %in% names(module), module$classification, NA)
    )
  })
  
  features_df <- do.call(rbind, lapply(features, function(x) {
    data.frame(x, stringsAsFactors = FALSE)
  }))
  
  return(features_df)
}

#' Plot Module Features
#'
#' Visualize module properties and characteristics
#'
#' @param module_analysis Data frame from analyze_module_features()
#' @param output_dir Directory to save plots
#'
#' @export
#'
#' @examples
#' plot_module_features(module_analysis, output_dir = "module_plots/")
plot_module_features <- function(module_analysis, output_dir = "module_plots") {
  if (!dir.exists(output_dir)) {
    dir.create(output_dir, recursive = TRUE)
  }
  
  # Size distribution
  p1 <- ggplot2::ggplot(module_analysis, ggplot2::aes(x = size)) +
    ggplot2::geom_histogram(bins = 30, fill = "steelblue", alpha = 0.7) +
    ggplot2::labs(title = "Module Size Distribution", x = "Module Size", y = "Count") +
    ggplot2::theme_minimal()
  
  # Density vs Size
  p2 <- ggplot2::ggplot(module_analysis, ggplot2::aes(x = size, y = density)) +
    ggplot2::geom_point(alpha = 0.7) +
    ggplot2::labs(title = "Density vs Module Size", x = "Module Size", y = "Density") +
    ggplot2::theme_minimal()
  
  # RNA content distribution
  if (all(c("rna_ratio", "classification") %in% colnames(module_analysis))) {
    p3 <- ggplot2::ggplot(module_analysis, ggplot2::aes(x = classification, y = rna_ratio)) +
      ggplot2::geom_boxplot(fill = "lightgreen", alpha = 0.7) +
      ggplot2::labs(title = "RNA Content by Classification", 
                    x = "Classification", y = "RNA Protein Ratio") +
      ggplot2::theme_minimal() +
      ggplot2::theme(axis.text.x = ggplot2::element_text(angle = 45, hjust = 1))
    
    ggplot2::ggsave(file.path(output_dir, "rna_content_by_classification.png"), p3, width = 8, height = 6)
  }
  
  # Save plots
  ggplot2::ggsave(file.path(output_dir, "module_size_distribution.png"), p1, width = 8, height = 6)
  ggplot2::ggsave(file.path(output_dir, "density_vs_size.png"), p2, width = 8, height = 6)
  
  message("Module feature plots saved to: ", output_dir)
}

#' Score Modules
#'
#' Calculate module activity scores for samples using various methods
#'
#' @param expression_data Matrix with genes as rows and samples as columns
#' @param modules List of modules
#' @param methods Character vector of scoring methods ("ssGSEA", "GSVA", "Zscore")
#' @param ssgsea_params Parameters for ssGSEA
#'
#' @return Matrix of module scores (modules x samples)
#'
#' @export
#'
#' @examples
#' expression_data <- matrix(rnorm(1000*20), nrow = 1000)
#' rownames(expression_data) <- paste0("GENE", 1:1000)
#' colnames(expression_data) <- paste0("SAMPLE", 1:20)
#'
#' module_scores <- score_modules(
#'   expression_data = expression_data,
#'   modules = classified_modules,
#'   methods = c("ssGSEA", "GSVA")
#' )
score_modules <- function(expression_data, modules, 
                      methods = c("ssGSEA", "GSVA", "Zscore"),
                      ssgsea_params = list(kcdf = "Gaussian")) {
  
  if (!is.matrix(expression_data) && !is.data.frame(expression_data)) {
    stop("expression_data must be a matrix or data frame")
  }
  
  if (!is.list(modules)) {
    stop("modules must be a list")
  }
  
  # Convert to matrix
  expr_mat <- as.matrix(expression_data)
  
  # Create gene sets from modules
  gene_sets <- lapply(modules, function(m) m$genes)
  names(gene_sets) <- sapply(modules, function(m) m$module_id)
  
  # Remove empty gene sets
  gene_sets <- gene_sets[sapply(gene_sets, length) > 0]
  
  if (length(gene_sets) == 0) {
    warning("No valid gene sets found")
    return(NULL)
  }
  
  # Initialize result list
  results <- list()
  
  # ssGSEA
  if ("ssGSEA" %in% methods) {
    if (!requireNamespace("GSVA", quietly = TRUE)) {
      warning("GSVA package not available for ssGSEA. Skipping.")
    } else {
      tryCatch({
        gset_idx <- GSEABase::GeneSetCollection(
          lapply(names(gene_sets), function(name) {
            GSEABase::GeneSet(geneIds = gene_sets[[name]], geneIdType = GSEABase::SymbolIdentifier())
          })
        )
        
        ssgsea_scores <- GSVA::gsva(
          expr_mat,
          gset_idx,
          method = "ssgsea",
          kcdf = ssgsea_params$kcdf %||% "Gaussian"
        )
        
        results$ssGSEA <- ssgsea_scores
      }, error = function(e) {
        warning("ssGSEA calculation failed: ", e$message)
      })
    }
  }
  
  # GSVA
  if ("GSVA" %in% methods) {
    if (!requireNamespace("GSVA", quietly = TRUE)) {
      warning("GSVA package not available for GSVA. Skipping.")
    } else {
      tryCatch({
        gset_idx <- GSEABase::GeneSetCollection(
          lapply(names(gene_sets), function(name) {
            GSEABase::GeneSet(geneIds = gene_sets[[name]], geneIdType = GSEABase::SymbolIdentifier())
          })
        )
        
        gsva_scores <- GSVA::gsva(
          expr_mat,
          gset_idx,
          method = "gsva"
        )
        
        results$GSVA <- gsva_scores
      }, error = function(e) {
        warning("GSVA calculation failed: ", e$message)
      })
    }
  }
  
  # Z-score
  if ("Zscore" %in% methods) {
    tryCatch({
      z_scores <- .calculate_z_scores(expr_mat, gene_sets)
      results$Zscore <- z_scores
    }, error = function(e) {
      warning("Z-score calculation failed: ", e$message)
    })
  }
  
  return(results)
}

#' Calculate Z-scores
#' @noRd
.calculate_z_scores <- function(expr_mat, gene_sets) {
  z_scores <- matrix(nrow = length(gene_sets), ncol = ncol(expr_mat))
  rownames(z_scores) <- names(gene_sets)
  colnames(z_scores) <- colnames(expr_mat)
  
  for (i in seq_along(gene_sets)) {
    genes_in_set <- intersect(gene_sets[[i]], rownames(expr_mat))
    
    if (length(genes_in_set) > 0) {
      for (j in seq_len(ncol(expr_mat))) {
        sample_vals <- expr_mat[genes_in_set, j]
        bg_vals <- expr_mat[, j]
        
        # Calculate z-score
        z_scores[i, j] <- (mean(sample_vals) - mean(bg_vals)) / sd(bg_vals)
      }
    }
  }
  
  return(z_scores)
}

#' Compare Module Scores
#'
#' Compare module scores across sample groups
#'
#' @param module_scores List of module score matrices
#' @param sample_groups Factor or vector indicating sample groups
#' @param method Scoring method to use
#'
#' @return List containing comparison results
#'
#' @export
#'
#' @examples
#' sample_groups <- factor(rep(c("Control", "Treatment"), each = 10))
#' score_comparison <- compare_module_scores(
#'   module_scores,
#'   sample_groups = sample_groups,
#'   method = "ssGSEA"
#' )
compare_module_scores <- function(module_scores, sample_groups, method = "ssGSEA") {
  if (!method %in% names(module_scores)) {
    stop("Method ", method, " not found in module_scores")
  }
  
  scores <- module_scores[[method]]
  
  if (length(sample_groups) != ncol(scores)) {
    stop("Length of sample_groups must match number of columns in scores")
  }
  
  # Create data frame for analysis
  score_df <- reshape2::melt(t(scores))
  colnames(score_df) <- c("Sample", "Module", "Score")
  score_df$Group <- sample_groups[match(score_df$Sample, names(sample_groups))]
  
  # Statistical comparison
  comparison_results <- list()
  
  for (module in unique(score_df$Module)) {
    module_data <- score_df[score_df$Module == module, ]
    
    # ANOVA or t-test depending on number of groups
    groups <- unique(module_data$Group)
    if (length(groups) == 2) {
      # t-test
      group1 <- module_data$Score[module_data$Group == groups[1]]
      group2 <- module_data$Score[module_data$Group == groups[2]]
      
      test_result <- t.test(group1, group2)
      comparison_results[[module]] <- list(
        test_type = "t-test",
        statistic = test_result$statistic,
        p_value = test_result$p.value,
        mean_group1 = mean(group1),
        mean_group2 = mean(group2)
      )
    } else if (length(groups) > 2) {
      # ANOVA
      formula <- Score ~ Group
      anova_result <- aov(formula, data = module_data)
      summary_result <- summary(anova_result)
      
      comparison_results[[module]] <- list(
        test_type = "ANOVA",
        f_statistic = summary_result[[1]]$`F value`[1],
        p_value = summary_result[[1]]$`Pr(>F)`[1]
      )
    }
  }
  
  return(list(
    scores = scores,
    comparison_data = score_df,
    statistical_tests = comparison_results
  ))
}

#' Plot Module Scores
#'
#' Visualize module scores across sample groups
#'
#' @param score_comparison Output from compare_module_scores()
#' @param output_dir Directory to save plots
#'
#' @export
#'
#' @examples
#' plot_module_scores(score_comparison, output_dir = "score_plots/")
plot_module_scores <- function(score_comparison, output_dir = "score_plots") {
  if (!dir.exists(output_dir)) {
    dir.create(output_dir, recursive = TRUE)
  }
  
  score_df <- score_comparison$comparison_data
  
  # Boxplot for each module
  p <- ggplot2::ggplot(score_df, ggplot2::aes(x = Module, y = Score, fill = Group)) +
    ggplot2::geom_boxplot(alpha = 0.7) +
    ggplot2::labs(title = "Module Scores by Group", 
                  x = "Module", y = "Score") +
    ggplot2::theme_minimal() +
    ggplot2::theme(axis.text.x = ggplot2::element_text(angle = 45, hjust = 1))
  
  ggplot2::ggsave(file.path(output_dir, "module_scores_by_group.png"), 
                 p, width = 10, height = 6)
  
  message("Module score plots saved to: ", output_dir)
}