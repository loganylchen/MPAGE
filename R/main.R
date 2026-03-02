#' Main Analysis Function
#'
#' Comprehensive main function to perform all types of gene pair analyses (intra-module, inter-module,
#' and module activity comparison) with default modules loaded from package if not provided.
#'
#' @param expression_data Gene expression matrix (genes as rows, samples as columns)
#' @param sample_condition Named vector indicating sample conditions ("case" or "control")
#' @param modules List of modules, each containing genes. If NULL, loads default modules from package.
#' @param module_pairs List of specific module pairs to compare. If NULL, compares all possible pairs.
#' @param min_samples Minimum samples per condition for valid analysis (default: 3)
#' @param activity_methods List of methods to calculate module activity (default: c("gsva", "ssgsea", "zscore", "plage"))
#'
#' @return List containing comprehensive results from all three analysis types with metadata
#'
#' @export
#'
#' @examples
#' \dontrun{
#' # Basic usage with default modules
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
#' # Run main analysis with default modules
#' results <- main_analysis(expr_mat, conditions)
#'
#' # Run with custom modules
#' custom_modules <- list(
#'   module1 = c("gene_1", "gene_2", "gene_3"),
#'   module2 = c("gene_4", "gene_5", "gene_6")
#' )
#' results <- main_analysis(expr_mat, conditions, modules = custom_modules)
#'
#' # Run with specific module pairs
#' results <- main_analysis(expr_mat, conditions, module_pairs = list(c("module1", "module2")))
#' }
run_mpage <- function(expression_data, sample_condition, modules = NULL,
                      module_pairs = NULL, min_samples = 3,
                      activity_methods = c("gsva", "ssgsea", "zscore", "plage"), temp_dir = "./temp") {
  # Input validation
  if (!is.matrix(expression_data) && !is.data.frame(expression_data)) {
    stop("expression_data must be a matrix or data frame")
  }

  if (!is.character(sample_condition) || length(unique(sample_condition)) != 2) {
    stop("sample_condition must have exactly two conditions: 'case' and 'control'")
  }

  # Ensure conditions are named
  if (is.null(names(sample_condition))) {
    names(sample_condition) <- colnames(expression_data)
  }

  # Convert to matrix
  expr_mat <- as.matrix(expression_data)

  # Load default modules if not provided
  if (is.null(modules)) {
    message("Loading default modules from package...")
    modules <- load_default_modules()
    if (is.null(modules) || length(modules) == 0) {
      warning("No default modules available, please provide custom modules")
      return(list(error = "No modules provided or available"))
    }
  } else {
    message("Using provided modules...")
  }

  # Validate modules format
  if (!is.list(modules)) {
    stop("modules must be a list")
  }

  # Check module format and convert if necessary
  if (length(modules) > 0) {
    if (is.list(modules[[1]]) && "genes" %in% names(modules[[1]])) {
      message("Using structured module format...")
    } else {
      message("Converting simple list format to structured format...")
      # Convert simple format to structured format
      structured_modules <- lapply(names(modules), function(name) {
        list(
          module_id = name,
          genes = modules[[name]],
          size = length(modules[[name]]),
          algorithm = "USER_PROVIDED",
          depth = 1,
          target_proteins = character(0)
        )
      })
      names(structured_modules) <- names(modules)
      modules <- structured_modules
    }
  }

  # Display input summary
  message("=== Input Summary ===")
  message("Expression matrix dimensions: ", nrow(expr_mat), " genes x ", ncol(expr_mat), " samples")
  message("Sample conditions: ", table(sample_condition))
  message("Number of modules: ", length(modules))

  # Get module genes for analysis
  module_genes <- lapply(modules, function(x) intersect(x$genes, rownames(expr_mat)))

  # Filter modules with sufficient genes
  min_genes <- 3
  valid_modules <- module_genes[lengths(module_genes) >= min_genes]

  if (length(valid_modules) == 0) {
    stop("No modules with sufficient genes (minimum ", min_genes, " genes)")
  }

  message("Valid modules after filtering: ", length(valid_modules))
  message("Genes per valid module: ", paste(lengths(valid_modules), collapse = ", "))

  # Run comprehensive analysis
  message("\n=== Starting Comprehensive Analysis ===")

  tryCatch(
    {
      results <- final_integrated_analysis(
        expression_matrix = expr_mat,
        sample_condition = sample_condition,
        modules = modules,
        module_pairs = module_pairs,
        min_samples = min_samples,
        activity_method = activity_methods, temp_dir = temp_dir
      )

      # Add main function metadata
      results$main_analysis_metadata <- list(
        function_name = "main_analysis",
        timestamp = Sys.time(),
        r_version = R.version.string,
        package_version = "1.0.0",
        input_parameters = list(
          expression_data_dim = dim(expr_mat),
          sample_conditions = table(sample_condition),
          modules_provided = !is.null(modules),
          custom_module_pairs = !is.null(module_pairs),
          min_samples = min_samples,
          activity_methods = activity_methods
        )
      )

      message("\n=== Analysis Complete ===")
      message("Results summary:")
      message("  Intra-module pairs: ", nrow(results$intra_module_results))
      message("  Inter-module pairs: ", nrow(results$inter_module_results))
      message("  Module activity comparisons: ", nrow(results$module_activity_results))
      message(
        "  Significant results (adj p < 0.05): ",
        sum(results$intra_module_results$significant, na.rm = TRUE) +
          sum(results$inter_module_results$significant, na.rm = TRUE) +
          sum(results$module_activity_results$significant, na.rm = TRUE)
      )

      return(results)
    },
    error = function(e) {
      message("Error during analysis: ", e$message)
      return(list(error = e$message))
    }
  )
}

#' Load Default Modules
#'
#' Load predefined modules from the package
#'
#' @return List of modules in structured format
#' @keywords internal
#' @noRd
load_default_modules <- function() {
  # Check if package has default modules
  default_modules_file <- system.file("extdata", "rmp_modules.rds", package = "MPAGE")

  if (file.exists(default_modules_file)) {
    default_modules <- readRDS(default_modules_file)
    return(default_modules)
  } else {
    stop("No found default modules")
  }
}
