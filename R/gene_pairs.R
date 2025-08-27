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
  module_genes <- lapply(modules, function(module) {
    intersect(module$genes, rownames(expr_mat))
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
    module_names <- rownames(activity_scores)
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

  # Check if progress package is available
  has_progress <- requireNamespace("progress", quietly = TRUE)

  # Setup progress bar if available
  if (has_progress) {
    pb <- progress::progress_bar$new(
      format = "  [:bar] :percent (:current/:total) ETA: :eta",
      total = length(module_pairs), clear = FALSE, width = 60
    )
  }

  # Analyze each module pair
  results_list <- list()

  for (i in seq_along(module_pairs)) {
    pair <- module_pairs[[i]]

    if (has_progress) pb$tick()

    if (length(pair) != 2) {
      warning("Each module pair must contain exactly 2 module names")
      next
    }

    module1 <- pair[1]
    module2 <- pair[2]

    if (!module1 %in% rownames(activity_scores) || !module2 %in% rownames(activity_scores)) {
      warning(paste("Module pair", paste(pair, collapse = "-"), "not found"))
      next
    }

    # Get activity scores
    module1_scores <- activity_scores[module1, ]
    module2_scores <- activity_scores[module2, ]

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
    contingency_table[1, 1] <- sum(case_patterns$module1_greater) # module1 > module2 in case
    contingency_table[2, 1] <- sum(case_patterns$module2_greater) # module2 > module1 in case
    contingency_table[1, 2] <- sum(control_patterns$module1_greater) # module1 > module2 in control
    contingency_table[2, 2] <- sum(control_patterns$module2_greater) # module2 > module1 in control

    # Ensure table has no zero rows/columns
    if (any(contingency_table == 0)) {
      contingency_table <- contingency_table + 0.5 # Add small pseudocount
    }

    # Fisher exact test with warning suppression
    fisher_result <- suppressWarnings(fisher.test(contingency_table))

    # Calculate additional metrics
    case_ratio <- sum(case_patterns$module1_greater) / length(case_samples)
    control_ratio <- sum(control_patterns$module1_greater) / length(control_samples)

    case_cor <- cor(module1_case, module2_case, method = "pearson", use = "pairwise.complete.obs")
    control_cor <- cor(module1_control, module2_control, method = "pearson", use = "pairwise.complete.obs")

    # Calculate effect sizes
    module1_effect <- mean(module1_case) - mean(module1_control)
    module2_effect <- mean(module2_case) - mean(module2_control)

    results_list[[length(results_list) + 1]] <- data.frame(
      activity_method = activity_method,
      event1 = module1,
      event2 = module2,
      event_name = paste0(module1, ":", module2),
      p_value = fisher_result$p.value,
      odds_ratio = fisher_result$estimate,
      case_pattern_ratio = case_ratio,
      control_pattern_ratio = control_ratio,
      case_correlation = case_cor,
      control_correlation = control_cor,
      correlation_change = case_cor - control_cor,
      event1_effect_size = module1_effect,
      event2_effect_size = module2_effect,
      event_interaction = module1_effect * module2_effect,
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

#' Intra-Module Gene Pair Analysis
#'
#' Analyze gene pairs within the same module to identify differential expression patterns
#' between case and control conditions using Fisher exact test.
#'
#' @param expression_matrix Gene expression matrix (genes as rows, samples as columns)
#' @param sample_condition Named vector indicating sample conditions ("case" or "control")
#' @param module_genes Vector of gene symbols in the module
#' @param min_samples Minimum samples per condition for valid analysis (default: 3)
#'
#' @return Data frame with intra-module gene pair analysis results
#'
#' @export
#'
#' @examples
#' \dontrun{
#' # Create example data
#' expr_mat <- matrix(rnorm(100), nrow = 10, ncol = 10)
#' colnames(expr_mat) <- paste0("sample_", 1:10)
#' rownames(expr_mat) <- paste0("gene_", 1:10)
#'
#' # Define conditions
#' conditions <- c(
#'   sample_1 = "case", sample_2 = "case", sample_3 = "case",
#'   sample_4 = "case", sample_5 = "case",
#'   sample_6 = "control", sample_7 = "control", sample_8 = "control",
#'   sample_9 = "control", sample_10 = "control"
#' )
#'
#' # Analyze intra-module pairs
#' results <- intra_gene_pair(expr_mat, conditions, c("gene_1", "gene_2", "gene_3"))
#' }
intra_gene_pair <- function(expression_matrix, sample_condition, modules, min_samples = 3) {
  # Input validation
  if (!is.matrix(expression_matrix) && !is.data.frame(expression_matrix)) {
    stop("expression_matrix must be a matrix or data frame")
  }

  if (!is.character(sample_condition) || length(unique(sample_condition)) != 2) {
    stop("sample_condition must have exactly two conditions: 'case' and 'control'")
  }

  # Convert to matrix
  expr_mat <- as.matrix(expression_matrix)

  # Ensure conditions are named
  if (is.null(names(sample_condition))) {
    names(sample_condition) <- colnames(expression_matrix)
  }

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
  module_genes <- lapply(modules, function(module) {
    intersect(module$genes, rownames(expr_mat))
  })

  # Remove empty modules and modules with too few genes
  min_genes <- 3
  module_genes <- module_genes[sapply(module_genes, length) >= min_genes]



  # Function to classify expression pattern for Fisher test
  classify_expression_pattern <- function(expr1, expr2) {
    gene1_greater <- expr1 > expr2
    gene2_greater <- expr2 > expr1
    equal_expression <- expr1 == expr2

    return(list(
      gene1_greater = gene1_greater,
      gene2_greater = gene2_greater,
      equal_expression = equal_expression
    ))
  }

  # Check if progress package is available
  has_progress <- requireNamespace("progress", quietly = TRUE)

  # Calculate total number of gene pairs for progress bar
  total_pairs <- 0
  for (module_name in names(module_genes)) {
    available_genes <- intersect(module_genes[[module_name]], rownames(expr_mat))
    if (length(available_genes) >= 2) {
      total_pairs <- total_pairs + choose(length(available_genes), 2)
    }
  }

  # Setup progress bar if available
  if (has_progress && total_pairs > 0) {
    pb <- progress::progress_bar$new(
      format = "  [:bar] :percent (:current/:total) ETA: :eta",
      total = total_pairs, clear = FALSE, width = 60
    )
  }

  # Analyze each gene pair
  results_list <- list()
  processed_pairs <- 0

  for (module_name in names(module_genes)) {
    # Filter module genes for available genes
    available_genes <- intersect(module_genes[[module_name]], rownames(expr_mat))

    if (length(available_genes) < 2) {
      warning("Fewer than 2 genes available in module")
      next
    }

    # Generate all possible gene pairs within the module
    gene_pairs <- combn(available_genes, 2, simplify = FALSE)

    for (pair in gene_pairs) {
      processed_pairs <- processed_pairs + 1
      if (has_progress && total_pairs > 0) pb$tick()

      gene1 <- pair[1]
      gene2 <- pair[2]

      # Get expression values
      expr1 <- expr_mat[gene1, ]
      expr2 <- expr_mat[gene2, ]

      # Separate by condition
      expr1_case <- expr1[case_samples]
      expr2_case <- expr2[case_samples]
      expr1_control <- expr1[control_samples]
      expr2_control <- expr2[control_samples]

      # Classify expression patterns
      case_patterns <- classify_expression_pattern(expr1_case, expr2_case)
      control_patterns <- classify_expression_pattern(expr1_control, expr2_control)

      # Create 2x2 contingency table for Fisher exact test
      # Rows: gene1 > gene2 (TRUE/FALSE)
      # Columns: condition (case/control)
      contingency_table <- matrix(0, nrow = 2, ncol = 2)
      rownames(contingency_table) <- c("gene1_greater", "gene1_not_greater")
      colnames(contingency_table) <- c("case", "control")

      # Fill the table
      contingency_table[1, 1] <- sum(case_patterns$gene1_greater) # gene1 > gene2 in case
      contingency_table[2, 1] <- sum(case_patterns$gene2_greater) # gene2 > gene1 in case
      contingency_table[1, 2] <- sum(control_patterns$gene1_greater) # gene1 > gene2 in control
      contingency_table[2, 2] <- sum(control_patterns$gene2_greater) # gene2 > gene1 in control

      # Ensure table has no zero rows/columns
      if (any(contingency_table == 0)) {
        contingency_table <- contingency_table + 0.5 # Add small pseudocount
      }

      # Fisher exact test with warning suppression
      fisher_result <- suppressWarnings(fisher.test(contingency_table))

      # Calculate additional metrics
      case_ratio <- sum(case_patterns$gene1_greater) / length(case_samples)
      control_ratio <- sum(control_patterns$gene1_greater) / length(control_samples)

      case_cor <- cor(expr1_case, expr2_case, method = "pearson", use = "pairwise.complete.obs")
      control_cor <- cor(expr1_control, expr2_control, method = "pearson", use = "pairwise.complete.obs")

      # Calculate effect sizes
      gene1_effect <- mean(expr1_case) - mean(expr1_control)
      gene2_effect <- mean(expr2_case) - mean(expr2_control)

      results_list[[length(results_list) + 1]] <- data.frame(
        activity_method = "Intra-Module",
        event1 = gene1,
        event2 = gene2,
        event_name = module_name,
        p_value = fisher_result$p.value,
        odds_ratio = fisher_result$estimate,
        case_pattern_ratio = case_ratio,
        control_pattern_ratio = control_ratio,
        case_correlation = case_cor,
        control_correlation = control_cor,
        correlation_change = case_cor - control_cor,
        event1_effect_size = gene1_effect,
        event2_effect_size = gene2_effect,
        event_interaction = gene1_effect * gene2_effect,
        stringsAsFactors = FALSE
      )
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

#' Inter-Module Gene Pair Analysis
#'
#' Analyze gene pairs between different modules to identify differential expression patterns
#' between case and control conditions using Fisher exact test.
#'
#' @param expression_matrix Gene expression matrix (genes as rows, samples as columns)
#' @param sample_condition Named vector indicating sample conditions ("case" or "control")
#' @param module1_genes Vector of gene symbols in the first module
#' @param module2_genes Vector of gene symbols in the second module
#' @param min_samples Minimum samples per condition for valid analysis (default: 3)
#'
#' @return Data frame with inter-module gene pair analysis results
#'
#' @export
#'
#' @examples
#' \dontrun{
#' # Create example data
#' expr_mat <- matrix(rnorm(100), nrow = 10, ncol = 10)
#' colnames(expr_mat) <- paste0("sample_", 1:10)
#' rownames(expr_mat) <- paste0("gene_", 1:10)
#'
#' # Define conditions
#' conditions <- c(
#'   sample_1 = "case", sample_2 = "case", sample_3 = "case",
#'   sample_4 = "case", sample_5 = "case",
#'   sample_6 = "control", sample_7 = "control", sample_8 = "control",
#'   sample_9 = "control", sample_10 = "control"
#' )
#'
#' # Analyze inter-module pairs
#' results <- inter_gene_pair(
#'   expr_mat,
#'   conditions,
#'   c("gene_1", "gene_2", "gene_3"),
#'   c("gene_4", "gene_5", "gene_6")
#' )
#' }
inter_gene_pair <- function(expression_matrix, sample_condition, modules, min_samples = 3) {
  # Input validation
  if (!is.matrix(expression_matrix) && !is.data.frame(expression_matrix)) {
    stop("expression_matrix must be a matrix or data frame")
  }

  if (!is.character(sample_condition) || length(unique(sample_condition)) != 2) {
    stop("sample_condition must have exactly two conditions: 'case' and 'control'")
  }

  # Convert to matrix
  expr_mat <- as.matrix(expression_matrix)

  # Ensure conditions are named
  if (is.null(names(sample_condition))) {
    names(sample_condition) <- colnames(expression_matrix)
  }

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
  module_genes <- lapply(modules, function(module) {
    intersect(module$genes, rownames(expr_mat))
  })

  module_names <- names(modules)
  module_pairs <- combn(module_names, 2, simplify = FALSE)

  # Check if progress package is available
  has_progress <- requireNamespace("progress", quietly = TRUE)

  # Calculate total number of gene pairs for progress bar
  total_pairs <- 0
  for (pair in module_pairs) {
    if (length(pair) == 2) {
      module1 <- pair[1]
      module2 <- pair[2]
      module1_available <- intersect(module_genes[[module1]], rownames(expr_mat))
      module2_available <- intersect(module_genes[[module2]], rownames(expr_mat))
      if (length(module1_available) >= 1 && length(module2_available) >= 1) {
        total_pairs <- total_pairs + length(module1_available) * length(module2_available)
      }
    }
  }

  # Setup progress bar if available
  if (has_progress && total_pairs > 0) {
    pb <- progress::progress_bar$new(
      format = "  [:bar] :percent (:current/:total) ETA: :eta",
      total = total_pairs, clear = FALSE, width = 60
    )
  }

  # Function to classify expression pattern for Fisher test
  classify_expression_pattern <- function(expr1, expr2) {
    gene1_greater <- expr1 > expr2
    gene2_greater <- expr2 > expr1
    equal_expression <- expr1 == expr2

    return(list(
      gene1_greater = gene1_greater,
      gene2_greater = gene2_greater,
      equal_expression = equal_expression
    ))
  }

  # Analyze each module pair
  results_list <- list()
  processed_pairs <- 0

  for (pair in module_pairs) {
    if (length(pair) != 2) {
      warning("Each module pair must contain exactly 2 module names")
      next
    }

    module1 <- pair[1]
    module2 <- pair[2]

    # Filter module genes for available genes
    module1_available <- intersect(module_genes[[module1]], rownames(expr_mat))
    module2_available <- intersect(module_genes[[module2]], rownames(expr_mat))

    if (length(module1_available) < 1 || length(module2_available) < 1) {
      warning(sprintf("No genes available in one or both modules:%s:%s", module1, module2))
      next
    }

    # Generate all possible gene pairs between modules
    gene_pairs <- expand.grid(module1 = module1_available, module2 = module2_available, stringsAsFactors = FALSE)

    for (i in 1:nrow(gene_pairs)) {
      processed_pairs <- processed_pairs + 1
      if (has_progress && total_pairs > 0) pb$tick()

      gene1 <- gene_pairs$module1[i]
      gene2 <- gene_pairs$module2[i]

      # Get expression values
      expr1 <- expr_mat[gene1, ]
      expr2 <- expr_mat[gene2, ]

      # Separate by condition
      expr1_case <- expr1[case_samples]
      expr2_case <- expr2[case_samples]
      expr1_control <- expr1[control_samples]
      expr2_control <- expr2[control_samples]

      # Classify expression patterns
      case_patterns <- classify_expression_pattern(expr1_case, expr2_case)
      control_patterns <- classify_expression_pattern(expr1_control, expr2_control)

      # Create 2x2 contingency table for Fisher exact test
      # Rows: gene1 > gene2 (TRUE/FALSE)
      # Columns: condition (case/control)
      contingency_table <- matrix(0, nrow = 2, ncol = 2)
      rownames(contingency_table) <- c("gene1_greater", "gene1_not_greater")
      colnames(contingency_table) <- c("case", "control")

      # Fill the table
      contingency_table[1, 1] <- sum(case_patterns$gene1_greater) # gene1 > gene2 in case
      contingency_table[2, 1] <- sum(case_patterns$gene2_greater) # gene2 > gene1 in case
      contingency_table[1, 2] <- sum(control_patterns$gene1_greater) # gene1 > gene2 in control
      contingency_table[2, 2] <- sum(control_patterns$gene2_greater) # gene2 > gene1 in control

      # Ensure table has no zero rows/columns
      if (any(contingency_table == 0)) {
        contingency_table <- contingency_table + 0.5 # Add small pseudocount
      }

      # Fisher exact test with warning suppression
      fisher_result <- suppressWarnings(fisher.test(contingency_table))

      # Calculate additional metrics
      case_ratio <- sum(case_patterns$gene1_greater) / length(case_samples)
      control_ratio <- sum(control_patterns$gene1_greater) / length(control_samples)

      case_cor <- cor(expr1_case, expr2_case, method = "pearson", use = "pairwise.complete.obs")
      control_cor <- cor(expr1_control, expr2_control, method = "pearson", use = "pairwise.complete.obs")

      # Calculate effect sizes
      gene1_effect <- mean(expr1_case) - mean(expr1_control)
      gene2_effect <- mean(expr2_case) - mean(expr2_control)

      results_list[[length(results_list) + 1]] <- data.frame(
        activity_method = "Inter-Module",
        event1 = gene1,
        event2 = gene2,
        event_name = paste0(module1, ":", module2),
        p_value = fisher_result$p.value,
        odds_ratio = fisher_result$estimate,
        case_pattern_ratio = case_ratio,
        control_pattern_ratio = control_ratio,
        case_correlation = case_cor,
        control_correlation = control_cor,
        correlation_change = case_cor - control_cor,
        event1_effect_size = gene1_effect,
        event2_effect_size = gene2_effect,
        event_interaction = gene1_effect * gene2_effect,
        stringsAsFactors = FALSE
      )
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
calculate_module_activity <- function(expression_matrix, module_genes, method = "gsva") {
  # Check required packages
  if (method == "gsva" && !requireNamespace("GSVA", quietly = TRUE)) {
    stop("GSVA package is required for GSVA method. Please install with: BiocManager::install('GSVA')")
  }

  if (method %in% c("ssgsea", "zscore", "plage") && !requireNamespace("GSVA", quietly = TRUE)) {
    stop("GSVA package is required for these methods. Please install with: BiocManager::install('GSVA')")
  }

  expr_mat <- as.matrix(expression_matrix)



  if (length(module_genes) == 0) {
    warning("No modules with available genes")
    return(NULL)
  }

  tryCatch(
    {
      # Calculate activity scores based on method
      switch(method,
        "gsva" = {
          params <- GSVA::gsvaParam(expr_mat, module_genes)
        },
        "ssgsea" = {
          params <- GSVA::ssgseaParam(expr_mat, module_genes)
        },
        "zscore" = {
          params <- GSVA::zscoreParam(expr_mat, module_genes)
        },
        "plage" = {
          params <- GSVA::plageParam(expr_mat, module_genes)
        },
        {
          warning("Unknown method, using GSVA as default")
          params <- GSVA::gsvaParam(expr_mat, module_genes)
        }
      )
      scores <- GSVA::gsva(params)
      return(scores)
    },
    error = function(e) {
      warning("Error calculating module activity: ", e$message)
      return(NULL)
    }
  )
}

#' Final Integrated Analysis
#'
#' Automatically run all types of gene pair analyses (intra-module, inter-module,
#' and module activity comparison) and return comprehensive results.
#'
#' @param expression_matrix Gene expression matrix (genes as rows, samples as columns)
#' @param sample_condition Named vector indicating sample conditions ("case" or "control")
#' @param modules List of modules, each containing genes (e.g., list(module1 = c("gene1", "gene2"), module2 = c("gene3", "gene4")))
#' @param module_pairs List of specific module pairs to compare. If NULL, compares all possible pairs.
#' @param min_samples Minimum samples per condition for valid analysis (default: 3)
#' @param activity_method List of methods to calculate module activity. Default: c("gsva", "ssgsea", "zscore", "plage")
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
                                      module_pairs = NULL,
                                      min_samples = 3, activity_method = c("gsva", "ssgsea", "zscore", "plage")) {
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

  # Extract gene lists from module structure
  if (is.list(modules) && length(modules) > 0) {
    # Check if modules have the structured format with $genes
    if (is.list(modules[[1]]) && "genes" %in% names(modules[[1]])) {
      # Handle structured format: list of lists with $module_id and $genes
      module_genes <- lapply(modules, function(x) intersect(x$genes, rownames(expr_mat)))
      names(module_genes) <- sapply(modules, function(x) {
        if (!is.null(x$module_id)) x$module_id[1] else names(modules)[which(modules == x)[1]]
      })
    } else {
      # Handle simple list format: list(module_name = c(genes...))
      module_genes <- lapply(modules, function(genes) intersect(genes, rownames(expr_mat)))
    }
  } else {
    stop("Invalid module format")
  }

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
  has_progress <- requireNamespace("progress", quietly = TRUE)
  # Initialize results list
  results <- list()

  # 1. Intra-module analysis
  message("\nRunning intra-module analysis...")
  intra_results <- intra_gene_pair(expr_mat, sample_condition, modules)



  if (length(intra_results) > 0) {
    results$intra_module_results <- intra_results
    message("Intra-module analysis completed: ", nrow(results$intra_module_results), " gene pairs")
  } else {
    results$intra_module_results <- data.frame()
    message("No intra-module results generated")
  }

  # 2. Inter-module analysis
  message("\nRunning inter-module analysis...")
  inter_results <- inter_gene_pair(expr_mat, sample_condition, modules)
  if (length(inter_results) > 0) {
    results$inter_module_results <- inter_results
    message("Inter-module analysis completed: ", nrow(results$inter_module_results), " gene pairs")
  } else {
    results$inter_module_results <- data.frame()
    message("No inter-module results generated")
  }


  # 3. Module activity comparison
  message("\nRunning module activity comparison...")

  # Use all specified activity methods
  activity_methods <- activity_method
  all_activity_results <- list()

  if (has_progress && length(activity_methods) > 0) {
    pb_activity <- progress::progress_bar$new(
      format = "  [:bar] :percent (:current/:total) ETA: :eta",
      total = length(activity_methods), clear = FALSE, width = 60
    )
  }

  processed_methods <- 0
  for (method in activity_methods) {
    processed_methods <- processed_methods + 1
    if (has_progress && length(activity_methods) > 0) pb_activity$tick()

    message("  Running module activity comparison with method: ", method)
    activity_results <- module_activity_comparison(
      expression_matrix = expr_mat,
      sample_condition = sample_condition,
      modules = modules,
      module_pairs = module_pairs,
      min_samples = min_samples,
      activity_method = method
    )

    if (nrow(activity_results) > 0) {
      # activity_results$comparison_type <- "module-activity"
      # activity_results$activity_method <- method
      all_activity_results[[method]] <- activity_results
      message("  Method ", method, " completed: ", nrow(activity_results), " module pairs")
    } else {
      message("  Method ", method, " completed: 0 module pairs")
    }
  }

  if (length(all_activity_results) > 0) {
    results$module_activity_results <- do.call(rbind, all_activity_results)
    message("Module activity comparison completed: ", nrow(results$module_activity_results), " total results across all methods")
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



  # Store summary
  results$summary <- list(
    overall_summary = summary_df,
    input_parameters = list(
      total_genes = nrow(expr_mat),
      total_samples = ncol(expr_mat),
      case_samples = length(case_samples),
      control_samples = length(control_samples),
      valid_modules = length(module_genes),
      analysis_methods = activity_methods
    )
  )

  # Add metadata to results
  results$metadata <- list(
    timestamp = Sys.time(),
    r_version = R.version.string,
    package_versions = list(
      gene_pairs = "1.0.0" # Assuming this is the package version
    )
  )

  message("\n=== Integrated Analysis Complete ===")
  message("Total results generated:")
  message("  Intra-module pairs: ", nrow(results$intra_module_results))
  message("  Inter-module pairs: ", nrow(results$inter_module_results))
  message("  Module activity comparisons: ", nrow(activity_results))

  return(results)
}
