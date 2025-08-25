#' Intra-gene Pair Analysis
#'
#' Analyze gene expression patterns within a module between case and control conditions
#' using Fisher exact test to detect significant expression pattern changes.
#'
#' @param expression_matrix Gene expression matrix (genes as rows, samples as columns)
#' @param sample_condition Named vector indicating sample conditions ("case" or "control")
#' @param module_genes Vector of gene names in the module to analyze
#' @param expression_threshold Quantile threshold for high/low expression classification (default: 0.5)
#' @param min_samples Minimum samples per condition for valid analysis (default: 3)
#'
#' @return Data frame with gene pair analysis results including Fisher test statistics
#'
#' @export
#'
#' @examples
#' \dontrun{
#' # Create example data
#' expr_mat <- matrix(rnorm(1000), nrow = 100, ncol = 10)
#' colnames(expr_mat) <- paste0("sample_", 1:10)
#' rownames(expr_mat) <- paste0("gene_", 1:100)
#' 
#' # Define conditions
#' conditions <- c(sample_1 = "case", sample_2 = "case", sample_3 = "case",
#'                  sample_4 = "case", sample_5 = "case",
#'                  sample_6 = "control", sample_7 = "control", sample_8 = "control",
#'                  sample_9 = "control", sample_10 = "control")
#' 
#' # Analyze intra-gene pairs
#' module_genes <- c("gene_1", "gene_2", "gene_3", "gene_4", "gene_5")
#' results <- intra_gene_pair(expr_mat, conditions, module_genes)
#' }
intra_gene_pair <- function(expression_matrix, sample_condition, module_genes, 
                            expression_threshold = 0.5, min_samples = 3) {
  
  # Input validation
  if (!is.matrix(expression_matrix) && !is.data.frame(expression_matrix)) {
    stop("expression_matrix must be a matrix or data frame")
  }
  
  if (!is.character(sample_condition) || length(unique(sample_condition)) != 2) {
    stop("sample_condition must have exactly two conditions: 'case' and 'control'")
  }
  
  # Ensure conditions are named
  if (is.null(names(sample_condition))) {
    names(sample_condition) <- colnames(expression_matrix)
  }
  
  # Convert to matrix and subset genes
  expr_mat <- as.matrix(expression_matrix)
  available_genes <- intersect(module_genes, rownames(expr_mat))
  
  if (length(available_genes) < 2) {
    warning("Fewer than 2 genes available in expression matrix")
    return(data.frame())
  }
  
  expr_mat <- expr_mat[available_genes, ]
  
  # Match samples between expression matrix and conditions
  common_samples <- intersect(colnames(expr_mat), names(sample_condition))
  expr_mat <- expr_mat[, common_samples]
  sample_condition <- sample_condition[common_samples]
  
  # Split samples by condition
  case_samples <- names(sample_condition)[sample_condition == "case"]
  control_samples <- names(sample_condition)[sample_condition == "control"]
  
  if (length(case_samples) < min_samples || length(control_samples) < min_samples) {
    warning("Insufficient samples in one or both conditions")
    return(data.frame())
  }
  
  # Generate all gene pairs
  gene_pairs <- t(combn(available_genes, 2))
  if (nrow(gene_pairs) == 0) {
    return(data.frame())
  }
  
  # Function to classify expression pattern for Fisher exact test
  classify_expression_pattern <- function(gene1_expr, gene2_expr, threshold) {
    # Calculate quantiles for high/low classification
    gene1_high <- gene1_expr > quantile(gene1_expr, threshold, na.rm = TRUE)
    gene2_high <- gene2_expr > quantile(gene2_expr, threshold, na.rm = TRUE)
    
    # Create 2x2 table for Fisher test
    # Pattern: gene1 high + gene2 low vs other combinations
    pattern1 <- gene1_high & !gene2_high
    pattern2 <- !gene1_high & gene2_high
    pattern3 <- gene1_high & gene2_high
    pattern4 <- !gene1_high & !gene2_high
    
    return(list(
      pattern1 = pattern1,
      pattern2 = pattern2,
      pattern3 = pattern3,
      pattern4 = pattern4
    ))
  }
  
  # Analyze each gene pair
  results_list <- list()
  
  for (i in 1:nrow(gene_pairs)) {
    gene1 <- gene_pairs[i, 1]
    gene2 <- gene_pairs[i, 2]
    
    # Get expression for both genes
    gene1_expr <- expr_mat[gene1, ]
    gene2_expr <- expr_mat[gene2, ]
    
    # Separate by condition
    gene1_case <- gene1_expr[case_samples]
    gene2_case <- gene2_expr[case_samples]
    gene1_control <- gene1_expr[control_samples]
    gene2_control <- gene2_expr[control_samples]
    
    # Classify expression patterns for each condition
    case_patterns <- classify_expression_pattern(gene1_case, gene2_case, expression_threshold)
    control_patterns <- classify_expression_pattern(gene1_control, gene2_control, expression_threshold)
    
    # Create contingency tables for Fisher exact test
    case_table <- table(case_patterns$pattern1, case_patterns$pattern2)
    control_table <- table(control_patterns$pattern1, control_patterns$pattern2)
    
    # Ensure 2x2 tables
    if (nrow(case_table) != 2 || ncol(case_table) != 2 || 
        nrow(control_table) != 2 || ncol(control_table) != 2) {
      next
    }
    
    # Fisher exact test
    fisher_result <- fisher.test(rbind(case_table, control_table))
    
    # Calculate additional metrics
    case_ratio <- sum(case_patterns$pattern1) / length(case_samples)
    control_ratio <- sum(control_patterns$pattern1) / length(control_samples)
    
    # Calculate correlation within each condition
    case_cor <- cor(gene1_case, gene2_case, method = "pearson", use = "pairwise.complete.obs")
    control_cor <- cor(gene1_control, gene2_control, method = "pearson", use = "pairwise.complete.obs")
    
    results_list[[length(results_list) + 1]] <- data.frame(
      gene1 = gene1,
      gene2 = gene2,
      p_value = fisher_result$p.value,
      odds_ratio = fisher_result$estimate,
      case_pattern_ratio = case_ratio,
      control_pattern_ratio = control_ratio,
      case_correlation = case_cor,
      control_correlation = control_cor,
      correlation_change = case_cor - control_cor,
      stringsAsFactors = FALSE
    )
  }
  
  if (length(results_list) == 0) {
    return(data.frame())
  }
  
  results <- do.call(rbind, results_list)
  results$adj_p_value <- p.adjust(results$p_value, method = "BH")
  results$significant <- results$adj_p_value < 0.05
  
  return(results)
}

#' Inter-gene Pair Analysis
#'
#' Analyze gene expression patterns between genes from different modules
#' using Fisher exact test to detect significant expression pattern changes.
#'
#' @param expression_matrix Gene expression matrix (genes as rows, samples as columns)
#' @param sample_condition Named vector indicating sample conditions ("case" or "control")
#' @param modules List of modules, each containing genes (e.g., list(module1 = c("gene1", "gene2"), module2 = c("gene3", "gene4")))
#' @param expression_threshold Quantile threshold for high/low expression classification (default: 0.5)
#' @param min_samples Minimum samples per condition for valid analysis (default: 3)
#'
#' @return Data frame with inter-module gene pair analysis results
#'
#' @export
#'
#' @examples
#' \dontrun{
#' # Create example data
#' expr_mat <- matrix(rnorm(1000), nrow = 100, ncol = 10)
#' colnames(expr_mat) <- paste0("sample_", 1:10)
#' rownames(expr_mat) <- paste0("gene_", 1:100)
#' 
#' # Define conditions
#' conditions <- c(sample_1 = "case", sample_2 = "case", sample_3 = "case",
#'                  sample_4 = "case", sample_5 = "case",
#'                  sample_6 = "control", sample_7 = "control", sample_8 = "control",
#'                  sample_9 = "control", sample_10 = "control")
#' 
#' # Define modules
#' modules <- list(
#'   module1 = c("gene_1", "gene_2", "gene_3"),
#'   module2 = c("gene_4", "gene_5", "gene_6"),
#'   module3 = c("gene_7", "gene_8", "gene_9")
#' )
#' 
#' # Analyze inter-gene pairs
#' results <- inter_gene_pair(expr_mat, conditions, modules)
#' }
inter_gene_pair <- function(expression_matrix, sample_condition, modules,
                           expression_threshold = 0.5, min_samples = 3) {
  
  # Input validation
  if (!is.matrix(expression_matrix) && !is.data.frame(expression_matrix)) {
    stop("expression_matrix must be a matrix or data frame")
  }
  
  if (!is.list(modules)) {
    stop("modules must be a list")
  }
  
  if (!is.character(sample_condition) || length(unique(sample_condition)) != 2) {
    stop("sample_condition must have exactly two conditions: 'case' and 'control'")
  }
  
  # Ensure conditions are named
  if (is.null(names(sample_condition))) {
    names(sample_condition) <- colnames(expression_matrix)
  }
  
  # Convert to matrix
  expr_mat <- as.matrix(expression_matrix)
  
  # Match samples
  common_samples <- intersect(colnames(expr_mat), names(sample_condition))
  expr_mat <- expr_mat[, common_samples]
  sample_condition <- sample_condition[common_samples]
  
  # Split by condition
  case_samples <- names(sample_condition)[sample_condition == "case"]
  control_samples <- names(sample_condition)[sample_condition == "control"]
  
  if (length(case_samples) < min_samples || length(control_samples) < min_samples) {
    warning("Insufficient samples in one or both conditions")
    return(data.frame())
  }
  
  # Filter modules for available genes
  module_genes <- lapply(modules, function(genes) {
    intersect(genes, rownames(expr_mat))
  })
  
  # Remove empty modules
  module_genes <- module_genes[sapply(module_genes, length) > 0]
  
  if (length(module_genes) < 2) {
    warning("Fewer than 2 modules with available genes")
    return(data.frame())
  }
  
  # Function to classify expression pattern
  classify_expression_pattern <- function(gene1_expr, gene2_expr, threshold) {
    gene1_high <- gene1_expr > quantile(gene1_expr, threshold, na.rm = TRUE)
    gene2_high <- gene2_expr > quantile(gene2_expr, threshold, na.rm = TRUE)
    
    pattern1 <- gene1_high & !gene2_high
    pattern2 <- !gene1_high & gene2_high
    pattern3 <- gene1_high & gene2_high
    pattern4 <- !gene1_high & !gene2_high
    
    return(list(
      pattern1 = pattern1,
      pattern2 = pattern2,
      pattern3 = pattern3,
      pattern4 = pattern4
    ))
  }
  
  # Analyze all inter-module gene pairs
  results_list <- list()
  module_names <- names(module_genes)
  
  for (i in 1:(length(module_names) - 1)) {
    for (j in (i + 1):length(module_names)) {
      module1 <- module_names[i]
      module2 <- module_names[j]
      
      genes1 <- module_genes[[module1]]
      genes2 <- module_genes[[module2]]
      
      # Generate all gene pairs between modules
      gene_pairs <- expand.grid(genes1, genes2, stringsAsFactors = FALSE)
      colnames(gene_pairs) <- c("gene1", "gene2")
      
      for (k in 1:nrow(gene_pairs)) {
        gene1 <- gene_pairs[k, "gene1"]
        gene2 <- gene_pairs[k, "gene2"]
        
        # Get expression data
        gene1_expr <- expr_mat[gene1, ]
        gene2_expr <- expr_mat[gene2, ]
        
        # Separate by condition
        gene1_case <- gene1_expr[case_samples]
        gene2_case <- gene2_expr[case_samples]
        gene1_control <- gene1_expr[control_samples]
        gene2_control <- gene2_expr[control_samples]
        
        # Classify patterns
        case_patterns <- classify_expression_pattern(gene1_case, gene2_case, expression_threshold)
        control_patterns <- classify_expression_pattern(gene1_control, gene2_control, expression_threshold)
        
        # Create contingency tables
        case_table <- table(case_patterns$pattern1, case_patterns$pattern2)
        control_table <- table(control_patterns$pattern1, control_patterns$pattern2)
        
        # Ensure 2x2 tables
        if (nrow(case_table) != 2 || ncol(case_table) != 2 || 
            nrow(control_table) != 2 || ncol(control_table) != 2) {
          next
        }
        
        # Fisher exact test
        fisher_result <- fisher.test(rbind(case_table, control_table))
        
        # Calculate metrics
        case_ratio <- sum(case_patterns$pattern1) / length(case_samples)
        control_ratio <- sum(control_patterns$pattern1) / length(control_samples)
        
        case_cor <- cor(gene1_case, gene2_case, method = "pearson", use = "pairwise.complete.obs")
        control_cor <- cor(gene1_control, gene2_control, method = "pearson", use = "pairwise.complete.obs")
        
        results_list[[length(results_list) + 1]] <- data.frame(
          gene1 = gene1,
          gene2 = gene2,
          module1 = module1,
          module2 = module2,
          p_value = fisher_result$p.value,
          odds_ratio = fisher_result$estimate,
          case_pattern_ratio = case_ratio,
          control_pattern_ratio = control_ratio,
          case_correlation = case_cor,
          control_correlation = control_cor,
          correlation_change = case_cor - control_cor,
          stringsAsFactors = FALSE
        )
      }
    }
  }
  
  if (length(results_list) == 0) {
    return(data.frame())
  }
  
  results <- do.call(rbind, results_list)
  results$adj_p_value <- p.adjust(results$p_value, method = "BH")
  results$significant <- results$adj_p_value < 0.05
  
  return(results)
}

#' Inter-module Pair Analysis
#'
#' Analyze module-level expression patterns between case and control conditions
#' using GSVA (Gene Set Variation Analysis) to calculate module activity scores,
#' followed by Fisher exact test to detect significant pattern changes.
#'
#' @param expression_matrix Gene expression matrix (genes as rows, samples as columns)
#' @param sample_condition Named vector indicating sample conditions ("case" or "control")
#' @param modules List of modules, each containing genes (e.g., list(module1 = c("gene1", "gene2"), module2 = c("gene3", "gene4")))
#' @param expression_threshold Quantile threshold for high/low activity classification (default: 0.5)
#' @param min_samples Minimum samples per condition for valid analysis (default: 3)
#' @param method GSVA method to use (default: "gsva")
#'
#' @return Data frame with inter-module analysis results
#'
#' @export
#'
#' @examples
#' \dontrun{
#' # Create example data
#' expr_mat <- matrix(rnorm(1000), nrow = 100, ncol = 10)
#' colnames(expr_mat) <- paste0("sample_", 1:10)
#' rownames(expr_mat) <- paste0("gene_", 1:100)
#' 
#' # Define conditions
#' conditions <- c(sample_1 = "case", sample_2 = "case", sample_3 = "case",
#'                  sample_4 = "case", sample_5 = "case",
#'                  sample_6 = "control", sample_7 = "control", sample_8 = "control",
#'                  sample_9 = "control", sample_10 = "control")
#' 
#' # Define modules
#' modules <- list(
#'   module1 = c("gene_1", "gene_2", "gene_3"),
#'   module2 = c("gene_4", "gene_5", "gene_6"),
#'   module3 = c("gene_7", "gene_8", "gene_9")
#' )
#' 
#' # Analyze inter-module pairs
#' results <- inter_module_pair(expr_mat, conditions, modules)
#' }
inter_module_pair <- function(expression_matrix, sample_condition, modules,
                            expression_threshold = 0.5, min_samples = 3, method = "gsva") {
  
  # Input validation
  if (!is.matrix(expression_matrix) && !is.data.frame(expression_matrix)) {
    stop("expression_matrix must be a matrix or data frame")
  }
  
  if (!is.list(modules)) {
    stop("modules must be a list")
  }
  
  if (!is.character(sample_condition) || length(unique(sample_condition)) != 2) {
    stop("sample_condition must have exactly two conditions: 'case' and 'control'")
  }
  
  # Check for GSVA package
  if (!requireNamespace("GSVA", quietly = TRUE)) {
    stop("GSVA package is required. Please install with: BiocManager::install('GSVA')")
  }
  
  # Ensure conditions are named
  if (is.null(names(sample_condition))) {
    names(sample_condition) <- colnames(expression_matrix)
  }
  
  # Convert to matrix
  expr_mat <- as.matrix(expression_matrix)
  
  # Match samples
  common_samples <- intersect(colnames(expr_mat), names(sample_condition))
  expr_mat <- expr_mat[, common_samples]
  sample_condition <- sample_condition[common_samples]
  
  # Split by condition
  case_samples <- names(sample_condition)[sample_condition == "case"]
  control_samples <- names(sample_condition)[sample_condition == "control"]
  
  if (length(case_samples) < min_samples || length(control_samples) < min_samples) {
    warning("Insufficient samples in one or both conditions")
    return(data.frame())
  }
  
  # Filter modules for available genes
  module_genes <- lapply(modules, function(genes) {
    intersect(genes, rownames(expr_mat))
  })
  
  # Remove empty modules
  module_genes <- module_genes[sapply(module_genes, length) > 3] # Minimum 3 genes for GSVA
  
  if (length(module_genes) < 2) {
    warning("Fewer than 2 modules with sufficient genes for GSVA")
    return(data.frame())
  }
  
  # Prepare gene sets for GSVA
  gene_sets <- module_genes
  
  tryCatch({
    # Calculate GSVA scores for all modules
    if (method == "gsva") {
      gsva_scores <- GSVA::gsva(expr_mat, gene_sets, method = "gsva", kcdf = "Gaussian")
    } else {
      gsva_scores <- GSVA::gsva(expr_mat, gene_sets, method = method)
    }
    
    # Separate scores by condition
    case_scores <- gsva_scores[, case_samples, drop = FALSE]
    control_scores <- gsva_scores[, control_samples, drop = FALSE]
    
    # Function to classify activity patterns for Fisher test
    classify_activity_pattern <- function(scores, threshold) {
      high_activity <- scores > quantile(scores, threshold, na.rm = TRUE)
      low_activity <- scores <= quantile(scores, threshold, na.rm = TRUE)
      return(list(
        high = high_activity,
        low = low_activity
      ))
    }
    
    # Analyze all module pairs
    results_list <- list()
    module_names <- colnames(gsva_scores)
    
    for (i in 1:(length(module_names) - 1)) {
      for (j in (i + 1):length(module_names)) {
        module1 <- module_names[i]
        module2 <- module_names[j]
        
        # Get activity scores
        module1_case <- case_scores[module1, ]
        module2_case <- case_scores[module2, ]
        module1_control <- control_scores[module1, ]
        module2_control <- control_scores[module2, ]
        
        # Classify activity patterns
        case_patterns1 <- classify_activity_pattern(module1_case, expression_threshold)
        case_patterns2 <- classify_activity_pattern(module2_case, expression_threshold)
        control_patterns1 <- classify_activity_pattern(module1_control, expression_threshold)
        control_patterns2 <- classify_activity_pattern(module2_control, expression_threshold)
        
        # Create 2x2 tables for Fisher test
        # Pattern: module1 high + module2 low vs other combinations
        case_table <- table(
          module1_high = case_patterns1$high,
          module2_low = case_patterns2$low
        )
        control_table <- table(
          module1_high = control_patterns1$high,
          module2_low = control_patterns2$low
        )
        
        # Ensure 2x2 tables
        if (all(dim(case_table) == c(2, 2)) && all(dim(control_table) == c(2, 2))) {
          # Fisher exact test
          fisher_result <- fisher.test(rbind(case_table, control_table))
          
          # Calculate metrics
          case_ratio <- sum(case_patterns1$high & case_patterns2$low) / length(case_samples)
          control_ratio <- sum(control_patterns1$high & control_patterns2$low) / length(control_samples)
          
          case_cor <- cor(module1_case, module2_case, method = "pearson", use = "pairwise.complete.obs")
          control_cor <- cor(module1_control, module2_control, method = "pearson", use = "pairwise.complete.obs")
          
          # Calculate effect size (difference in means)
          module1_effect <- mean(module1_case) - mean(module1_control)
          module2_effect <- mean(module2_case) - mean(module2_control)
          
          results_list[[length(results_list) + 1]] <- data.frame(
            module1 = module1,
            module2 = module2,
            p_value = fisher_result$p.value,
            odds_ratio = fisher_result$estimate,
            case_pattern_ratio = case_ratio,
            control_pattern_ratio = control_ratio,
            case_correlation = case_cor,
            control_correlation = control_cor,
            correlation_change = case_cor - control_cor,
            module1_effect_size = module1_effect,
            module2_effect_size = module2_effect,
            module_interaction = module1_effect * module2_effect,
            stringsAsFactors = FALSE
          )
        }
      }
    }
    
    if (length(results_list) == 0) {
      return(data.frame())
    }
    
    results <- do.call(rbind, results_list)
    results$adj_p_value <- p.adjust(results$p_value, method = "BH")
    results$significant <- results$adj_p_value < 0.05
    
    return(results)
    
  }, error = function(e) {
    warning("Error in GSVA analysis: ", e$message)
    return(data.frame())
  })
}

#' Helper function for NULL-coalescing
#' @noRd
`%||%` <- function(x, y) {
  if (is.null(x)) y else x
}