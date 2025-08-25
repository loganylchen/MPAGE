#' Module Activity Comparison
#'
#' Compare module activity scores between case and control conditions using Fisher exact test.
#' First calculates module activity for each sample using GSVA or specified method,
#' then tests if module1 activity is significantly greater than module2 activity.
#'
#' @param expression_matrix Gene expression matrix (genes as rows, samples as columns)
#' @param sample_condition Named vector indicating sample conditions ("case" or "control")
#' @param modules List of modules, each containing genes (e.g., list(module1 = c("gene1", "gene2"), module2 = c("gene3", "gene4")))
#' @param module_pairs List of module pairs to compare. If NULL, compares all possible pairs.
#' @param min_samples Minimum samples per condition for valid analysis (default: 3)
#' @param activity_method Method to calculate module activity (default: "gsva", options: "gsva", "ssgsea", "zscore", "plage")
#'
#' @return Data frame with module activity comparison results including Fisher test statistics
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
#' conditions <- c(
#'   sample_1 = "case", sample_2 = "case", sample_3 = "case",
#'   sample_4 = "case", sample_5 = "case",
#'   sample_6 = "control", sample_7 = "control", sample_8 = "control",
#'   sample_9 = "control", sample_10 = "control"
#' )
#'
#' # Define modules
#' modules <- list(
#'   module1 = c("gene_1", "gene_2", "gene_3"),
#'   module2 = c("gene_4", "gene_5", "gene_6"),
#'   module3 = c("gene_7", "gene_8", "gene_9")
#' )
#'
#' # Compare module activities
#' results <- module_activity_comparison(expr_mat, conditions, modules)
#'
#' # Compare specific module pairs
#' pairs <- list(c("module1", "module2"), c("module2", "module3"))
#' results <- module_activity_comparison(expr_mat, conditions, modules, module_pairs = pairs)
#' }
module_activity_comparison <- function(expression_matrix, sample_condition, modules,
                                   module_pairs = NULL, min_samples = 3,
                                   activity_method = "gsva") {
  
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
  
  # Remove empty modules and modules with too few genes
  min_genes <- 3
  module_genes <- module_genes[sapply(module_genes, length) >= min_genes]
  
  if (length(module_genes) < 2) {
    warning("Fewer than 2 modules with sufficient genes (minimum ", min_genes, " genes)")
    return(data.frame())
  }
  
  # Calculate module activity scores
  activity_scores <- calculate_module_activity(expr_mat, module_genes, method = activity_method)
  
  if (is.null(activity_scores) || nrow(activity_scores) == 0) {
    warning("Could not calculate module activity scores")
    return(data.frame())
  }
  
  # Create module pairs if not provided
  if (is.null(module_pairs)) {
    module_names <- colnames(activity_scores)
    module_pairs <- combn(module_names, 2, simplify = FALSE)
  }
  
  # Function to classify activity patterns for Fisher test
  classify_activity_pattern <- function(scores1, scores2) {
    module1_greater <- scores1 > scores2
    module2_greater <- scores2 > scores1
    equal_activity <- scores1 == scores2
    
    return(list(
      module1_greater = module1_greater,
      module2_greater = module2_greater,
      equal_activity = equal_activity
    ))
  }
  
  # Analyze each module pair
  results_list <- list()
  
  for (pair in module_pairs) {
    if (length(pair) != 2) {
      warning("Each module pair must contain exactly 2 module names")
      next
    }
    
    module1 <- pair[1]
    module2 <- pair[2]
    
    if (!module1 %in% colnames(activity_scores) || !module2 %in% colnames(activity_scores)) {
      warning(paste("Module pair", paste(pair, collapse = "-"), "not found"))
      next
    }
    
    # Get activity scores
    module1_scores <- activity_scores[, module1]
    module2_scores <- activity_scores[, module2]
    
    # Separate by condition
    module1_case <- module1_scores[case_samples]
    module2_case <- module2_scores[case_samples]
    module1_control <- module1_scores[control_samples]
    module2_control <- module2_scores[control_samples]
    
    # Classify activity patterns
    case_patterns <- classify_activity_pattern(module1_case, module2_case)
    control_patterns <- classify_activity_pattern(module1_control, module2_control)
    
    # Create 2x2 contingency table for Fisher exact test
    # Rows: module1 > module2 (TRUE/FALSE)
    # Columns: condition (case/control)
    contingency_table <- matrix(0, nrow = 2, ncol = 2)
    rownames(contingency_table) <- c("module1_greater", "module1_not_greater")
    colnames(contingency_table) <- c("case", "control")
    
    # Fill the table
    contingency_table[1, 1] <- sum(case_patterns$module1_greater)  # module1 > module2 in case
    contingency_table[2, 1] <- sum(case_patterns$module2_greater)  # module2 > module1 in case
    contingency_table[1, 2] <- sum(control_patterns$module1_greater)  # module1 > module2 in control
    contingency_table[2, 2] <- sum(control_patterns$module2_greater)  # module2 > module1 in control
    
    # Ensure table has no zero rows/columns
    if (any(contingency_table == 0)) {
      contingency_table <- contingency_table + 0.5  # Add small pseudocount
    }
    
    # Fisher exact test
    fisher_result <- fisher.test(contingency_table)
    
    # Calculate additional metrics
    case_ratio <- sum(case_patterns$module1_greater) / length(case_samples)
    control_ratio <- sum(control_patterns$module1_greater) / length(control_samples)
    
    case_cor <- cor(module1_case, module2_case, method = "pearson", use = "pairwise.complete.obs")
    control_cor <- cor(module1_control, module2_control, method = "pearson", use = "pairwise.complete.obs")
    
    # Calculate effect sizes
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
  
  if (length(results_list) == 0) {
    return(data.frame())
  }
  
  results <- do.call(rbind, results_list)
  results$adj_p_value <- p.adjust(results$p_value, method = "BH")
  results$significant <- results$adj_p_value < 0.05
  
  return(results)
}

#' Calculate Module Activity Scores
#'
#' Calculate module activity scores for each sample using various methods.
#'
#' @param expression_matrix Gene expression matrix (genes as rows, samples as columns)
#' @param modules List of modules, each containing genes
#' @param method Method to calculate module activity ("gsva", "ssgsea", "zscore", "plage")
#'
#' @return Matrix with module activity scores (modules as rows, samples as columns)
#'
#' @noRd
calculate_module_activity <- function(expression_matrix, modules, method = "gsva") {
  
  # Check required packages
  if (method == "gsva" && !requireNamespace("GSVA", quietly = TRUE)) {
    stop("GSVA package is required for GSVA method. Please install with: BiocManager::install('GSVA')")
  }
  
  if (method %in% c("ssgsea", "zscore", "plage") && !requireNamespace("GSVA", quietly = TRUE)) {
    stop("GSVA package is required for these methods. Please install with: BiocManager::install('GSVA')")
  }
  
  expr_mat <- as.matrix(expression_matrix)
  
  # Filter modules for available genes
  module_genes <- lapply(modules, function(genes) {
    intersect(genes, rownames(expr_mat))
  })
  
  # Remove empty modules
  module_genes <- module_genes[sapply(module_genes, length) > 0]
  
  if (length(module_genes) == 0) {
    warning("No modules with available genes")
    return(NULL)
  }
  
  tryCatch({
    # Calculate activity scores based on method
    switch(method,
      "gsva" = {
        scores <- GSVA::gsva(expr_mat, module_genes, method = "gsva", kcdf = "Gaussian")
      },
      "ssgsea" = {
        scores <- GSVA::gsva(expr_mat, module_genes, method = "ssgsea")
      },
      "zscore" = {
        scores <- GSVA::gsva(expr_mat, module_genes, method = "zscore")
      },
      "plage" = {
        scores <- GSVA::gsva(expr_mat, module_genes, method = "plage")
      },
      {
        warning("Unknown method, using GSVA as default")
        scores <- GSVA::gsva(expr_mat, module_genes, method = "gsva", kcdf = "Gaussian")
      }
    )
    
    return(scores)
    
  }, error = function(e) {
    warning("Error calculating module activity: ", e$message)
    return(NULL)
  })
}

#' Final Integrated Analysis
#'
#' Automatically run all types of gene pair analyses (intra-module, inter-module, 
#' and module activity comparison) and return comprehensive results.
#'
#' @param expression_matrix Gene expression matrix (genes as rows, samples as columns)
#' @param sample_condition Named vector indicating sample conditions ("case" or "control")
#' @param modules List of modules, each containing genes (e.g., list(module1 = c("gene1", "gene2"), module2 = c("gene3", "gene4")))
#' @param gene_pairs List of specific gene pairs to analyze. If NULL, analyzes all possible pairs.
#' @param module_pairs List of specific module pairs to compare. If NULL, compares all possible pairs.
#' @param min_samples Minimum samples per condition for valid analysis (default: 3)
#' @param activity_method Method to calculate module activity (default: "gsva")
#'
#' @return List containing results from all three analysis types with metadata
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
#' conditions <- c(
#'   sample_1 = "case", sample_2 = "case", sample_3 = "case",
#'   sample_4 = "case", sample_5 = "case",
#'   sample_6 = "control", sample_7 = "control", sample_8 = "control",
#'   sample_9 = "control", sample_10 = "control"
#' )
#'
#' # Define modules
#' modules <- list(
#'   module1 = c("gene_1", "gene_2", "gene_3"),
#'   module2 = c("gene_4", "gene_5", "gene_6"),
#'   module3 = c("gene_7", "gene_8", "gene_9")
#' )
#'
#' # Run integrated analysis
#' final_results <- final_integrated_analysis(expr_mat, conditions, modules)
#' 
#' # View results
#' print(final_results$summary)
#' head(final_results$intra_module_results)
#' head(final_results$inter_module_results)
#' head(final_results$module_activity_results)
#' }
final_integrated_analysis <- function(expression_matrix, sample_condition, modules,
                                   gene_pairs = NULL, module_pairs = NULL, 
                                   min_samples = 3, activity_method = "gsva") {
  
  # Validate inputs
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
  
  # Check minimum samples
  if (length(case_samples) < min_samples || length(control_samples) < min_samples) {
    stop("Insufficient samples in one or both conditions")
  }
  
  # Filter modules for available genes
  module_genes <- lapply(modules, function(genes) {
    intersect(genes, rownames(expr_mat))
  })
  
  # Remove empty modules and modules with too few genes
  min_genes <- 3
  module_genes <- module_genes[sapply(module_genes, length) >= min_genes]
  
  if (length(module_genes) < 1) {
    stop("No modules with sufficient genes (minimum ", min_genes, " genes)")
  }
  
  message("Starting integrated analysis...")
  message("Input genes: ", nrow(expr_mat))
  message("Input samples: ", ncol(expr_mat))
  message("Valid modules: ", length(module_genes))
  message("Case samples: ", length(case_samples))
  message("Control samples: ", length(control_samples))
  
  # Initialize results list
  results <- list()
  
  # 1. Intra-module analysis
  message("\nRunning intra-module analysis...")
  intra_results <- list()
  
  for (module_name in names(module_genes)) {
    module_genes_list <- module_genes[[module_name]]
    
    if (length(module_genes_list) >= 2) {
      # Run intra-gene pair analysis for this module
      intra_df <- intra_gene_pair(
        expression_matrix = expr_mat,
        sample_condition = sample_condition,
        module_genes = module_genes_list
      )
      
      if (nrow(intra_df) > 0) {
        intra_df$comparison_type <- "intra-module"
        intra_df$module1 <- module_name
        intra_df$module2 <- module_name
        intra_results[[module_name]] <- intra_df
      }
    }
  }
  
  if (length(intra_results) > 0) {
    results$intra_module_results <- do.call(rbind, intra_results)
    message("Intra-module analysis completed: ", nrow(results$intra_module_results), " gene pairs")
  } else {
    results$intra_module_results <- data.frame()
    message("No intra-module results generated")
  }
  
  # 2. Inter-module analysis
  message("\nRunning inter-module analysis...")
  inter_results <- list()
  
  # Create module pairs if not provided
  if (is.null(gene_pairs)) {
    # Generate all possible gene pairs between different modules
    module_names <- names(module_genes)
    
    for (i in 1:(length(module_names) - 1)) {
      for (j in (i + 1):length(module_names)) {
        module1_name <- module_names[i]
        module2_name <- module_names[j]
        
        module1_genes <- module_genes[[module1_name]]
        module2_genes <- module_genes[[module2_name]]
        
        # Run inter-gene pair analysis
        inter_df <- inter_gene_pair(
          expression_matrix = expr_mat,
          sample_condition = sample_condition,
          module1_genes = module1_genes,
          module2_genes = module2_genes
        )
        
        if (nrow(inter_df) > 0) {
          inter_df$comparison_type <- "inter-module"
          inter_df$module1 <- module1_name
          inter_df$module2 <- module2_name
          inter_results[[paste(module1_name, module2_name, sep = "_")]] <- inter_df
        }
      }
    }
  } else {
    # Use provided gene pairs
    message("Using provided gene pairs for inter-module analysis")
    # This would need custom handling based on gene_pairs structure
  }
  
  if (length(inter_results) > 0) {
    results$inter_module_results <- do.call(rbind, inter_results)
    message("Inter-module analysis completed: ", nrow(results$inter_module_results), " gene pairs")
  } else {
    results$inter_module_results <- data.frame()
    message("No inter-module results generated")
  }
  
  # 3. Module activity comparison
  message("\nRunning module activity comparison...")
  activity_results <- module_activity_comparison(
    expression_matrix = expr_mat,
    sample_condition = sample_condition,
    modules = module_genes,
    module_pairs = module_pairs,
    min_samples = min_samples,
    activity_method = activity_method
  )
  
  if (nrow(activity_results) > 0) {
    activity_results$comparison_type <- "module-activity"
    results$module_activity_results <- activity_results
    message("Module activity comparison completed: ", nrow(activity_results), " module pairs")
  } else {
    results$module_activity_results <- data.frame()
    message("No module activity results generated")
  }
  
  # Create comprehensive summary
  summary_df <- data.frame(
    analysis_type = c("intra-module", "inter-module", "module-activity"),
    total_pairs = c(
      nrow(results$intra_module_results),
      nrow(results$inter_module_results), 
      nrow(results$module_activity_results)
    ),
    significant_pairs = c(
      sum(results$intra_module_results$significant, na.rm = TRUE),
      sum(results$inter_module_results$significant, na.rm = TRUE),
      sum(results$module_activity_results$significant, na.rm = TRUE)
    ),
    stringsAsFactors = FALSE
  )
  
  # Module involvement summary
  all_modules <- unique(c(
    results$intra_module_results$module1,
    results$inter_module_results$module1,
    results$inter_module_results$module2,
    results$module_activity_results$module1,
    results$module_activity_results$module2
  ))
  all_modules <- all_modules[all_modules != ""]
  
  module_summary <- data.frame(
    module = all_modules,
    involved_in_intra = all_modules %in% unique(results$intra_module_results$module1),
    involved_in_inter = all_modules %in% unique(c(results$inter_module_results$module1, results$inter_module_results$module2)),
    involved_in_activity = all_modules %in% unique(c(results$module_activity_results$module1, results$module_activity_results$module2)),
    stringsAsFactors = FALSE
  )
  
  # Store summary
  results$summary <- list(
    overall_summary = summary_df,
    module_involvement = module_summary,
    input_parameters = list(
      total_genes = nrow(expr_mat),
      total_samples = ncol(expr_mat),
      case_samples = length(case_samples),
      control_samples = length(control_samples),
      valid_modules = length(module_genes),
      analysis_method = activity_method
    )
  )
  
  # Add metadata to results
  results$metadata <- list(
    timestamp = Sys.time(),
    r_version = R.version.string,
    package_versions = list(
      gene_pairs = "1.0.0"  # Assuming this is the package version
    )
  )
  
  message("\n=== Integrated Analysis Complete ===")
  message("Total results generated:")
  message("  Intra-module pairs: ", nrow(results$intra_module_results))
  message("  Inter-module pairs: ", nrow(results$inter_module_results))
  message("  Module activity comparisons: ", nrow(activity_results))
  message("  Modules involved: ", nrow(module_summary))
  
  return(results)
}