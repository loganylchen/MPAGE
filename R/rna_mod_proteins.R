#' Get RNA Modification Proteins
#'
#' Retrieve curated list of RNA modification-related proteins including writers, readers, and erasers
#' with detailed annotations and references. Supports both built-in database and user-provided files.
#' @param use_built_in Logical whether to use built-in data frame (default: TRUE)
#' @param file_path Character string specifying path to user-provided CSV file (required if use_built_in = FALSE)
#'
#' @return A data frame containing RNA modification proteins with columns:
#'   When use_built_in = TRUE: gene_symbol, modification_type, functional_role
#'   When use_built_in = FALSE: gene_symbol, modification_type, functional_role
#'
#' @details
#' This function can use either the built-in comprehensive database or a user-provided CSV file.
#'
#' **Built-in database columns**:
#'   * gene_symbol, modification_type, functional_role
#'
#' **User file required columns**:
#'   * gene_symbol, modification_type, functional_role
#'
#' @export
#'
#' @examples
#' # Use built-in database
#' m6a_proteins <- get_rna_mod_proteins(
#'   use_built_in = TRUE
#' )
#'
#' # Use user-provided file
#' user_proteins <- get_rna_mod_proteins(
#'   use_built_in = FALSE,
#'   file_path = "my_proteins.csv"
#' )
#'
get_rna_mod_proteins <- function(use_built_in = TRUE,
                                 file_path = NULL) {
  if (!use_built_in) {
    return(.get_user_rna_mod_proteins(
      file_path = file_path
    ))
  }

  # Use built-in database
  return(.get_built_in_rna_mod_proteins())
}

#' Get RNA Modification Proteins from User File
#' @noRd
.get_user_rna_mod_proteins <- function(file_path) {
  if (is.null(file_path)) {
    stop("file_path must be provided when use_built_in = FALSE")
  }

  if (!file.exists(file_path)) {
    stop("File not found: ", file_path)
  }

  # Read user-provided file
  proteins <- utils::read.csv(file_path, stringsAsFactors = FALSE)

  # Check required columns
  required_cols <- c("gene_symbol", "modification_type", "functional_role")
  missing_cols <- setdiff(required_cols, colnames(proteins))

  if (length(missing_cols) > 0) {
    stop("Missing required columns: ", paste(missing_cols, collapse = ", "))
  }




  # Return only the required columns for user-provided data
  result <- proteins[, c("gene_symbol", "modification_type", "functional_role")]
  rownames(result) <- NULL
  return(unique(result))
}

#' Get RNA Modification Proteins from Built-in Database
#' @noRd
.get_built_in_rna_mod_proteins <- function() {
  # Read the built-in database
  protein_data <- system.file("extdata", "rnamod_proteins.rds", package = "mpage")

  if (!file.exists(protein_data)) {
    stop("Built-in RNA modification protein database not found")
  }

  proteins <- readRDS(protein_data)

  # Return all columns from built-in database
  result <- proteins[, c(
    "gene_symbol", "modification_type", "functional_role"
  )]
  rownames(result) <- NULL
  return(unique(result))
}
