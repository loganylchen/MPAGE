#' Get RNA Modification Proteins
#'
#' Retrieve curated list of RNA modification-related proteins including writers, readers, and erasers
#' with detailed annotations and references. Supports both built-in database and user-provided files.
#'
#' @param modification_types Character vector of modification types to include. 
#'   Available: "m6A", "m5C", "Psi", "m1A", "m7G", "A-to-I"
#' @param include_writers Logical whether to include writer proteins
#' @param include_readers Logical whether to include reader proteins
#' @param include_erasers Logical whether to include eraser proteins
#' @param organism Character string specifying organism (default: "Homo sapiens")
#' @param data_source Character string specifying data source filter (default: "Curated_database")
#' @param use_built_in Logical whether to use built-in data frame (default: TRUE)
#' @param file_path Character string specifying path to user-provided CSV file (required if use_built_in = FALSE)
#'
#' @return A data frame containing RNA modification proteins with columns:
#'   When use_built_in = TRUE: gene_symbol, modification_type, functional_role, evidence_source, pmid, reference, organism
#'   When use_built_in = FALSE: gene_symbol, modification_type, functional_role, organism
#'
#' @details
#' This function can use either the built-in comprehensive database or a user-provided CSV file.
#' 
#' **Built-in database columns**:
#'   * ensembl_id, gene_id, gene_symbol, modification_type, functional_role, evidence_source, pmid, reference, organism
#'
#' **User file required columns**:
#'   * gene_symbol, modification_type, functional_role, organism
#'
#' @export
#'
#' @examples
#' # Use built-in database
#' m6a_proteins <- get_rna_mod_proteins(
#'   modification_types = "m6A",
#'   include_writers = TRUE,
#'   include_readers = TRUE,
#'   include_erasers = TRUE,
#'   use_built_in = TRUE
#' )
#' 
#' # Use user-provided file
#' user_proteins <- get_rna_mod_proteins(
#'   modification_types = "m6A",
#'   include_writers = TRUE,
#'   use_built_in = FALSE,
#'   file_path = "my_proteins.csv"
#' )
#' 
#' # Get only human m5C writers from built-in
#' m5c_writers <- get_rna_mod_proteins(
#'   modification_types = "m5C",
#'   include_writers = TRUE,
#'   include_readers = FALSE,
#'   include_erasers = FALSE,
#'   organism = "Homo sapiens",
#'   use_built_in = TRUE
#' )
get_rna_mod_proteins <- function(modification_types = c("m6A", "m5C", "Psi", "m1A", "m7G", "A-to-I"),
                                 include_writers = TRUE,
                                 include_readers = TRUE,
                                 include_erasers = TRUE,
                                 organism = "Homo sapiens",
                                 data_source = "Curated_database",
                                 use_built_in = TRUE,
                                 file_path = NULL) {
  
  if (!use_built_in) {
    return(.get_user_rna_mod_proteins(
      modification_types = modification_types,
      include_writers = include_writers,
      include_readers = include_readers,
      include_erasers = include_erasers,
      organism = organism,
      file_path = file_path
    ))
  }
  
  # Use built-in database
  return(.get_built_in_rna_mod_proteins(
    modification_types = modification_types,
    include_writers = include_writers,
    include_readers = include_readers,
    include_erasers = include_erasers,
    organism = organism,
    data_source = data_source
  ))
}

#' Get RNA Modification Proteins from User File
#' @noRd
.get_user_rna_mod_proteins <- function(modification_types, include_writers, include_readers, 
                                   include_erasers, organism, file_path) {
  
  if (is.null(file_path)) {
    stop("file_path must be provided when use_built_in = FALSE")
  }
  
  if (!file.exists(file_path)) {
    stop("File not found: ", file_path)
  }
  
  # Read user-provided file
  proteins <- utils::read.csv(file_path, stringsAsFactors = FALSE)
  
  # Check required columns
  required_cols <- c("gene_symbol", "modification_type", "functional_role", "organism")
  missing_cols <- setdiff(required_cols, colnames(proteins))
  
  if (length(missing_cols) > 0) {
    stop("Missing required columns: ", paste(missing_cols, collapse = ", "))
  }
  
  # Filter by organism
  if (!is.null(organism) && organism != "all") {
    proteins <- proteins[proteins$organism == organism, ]
  }
  
  # Filter by modification types
  if (!is.null(modification_types)) {
    proteins <- proteins[proteins$modification_type %in% modification_types, ]
  }
  
  # Filter by functional roles
  roles_to_include <- c()
  if (include_writers) roles_to_include <- c(roles_to_include, "writer")
  if (include_readers) roles_to_include <- c(roles_to_include, "reader")
  if (include_erasers) roles_to_include <- c(roles_to_include, "eraser")
  
  if (length(roles_to_include) > 0) {
    proteins <- proteins[proteins$functional_role %in% roles_to_include, ]
  }
  
  # Return only the required columns for user-provided data
  result <- proteins[, c("gene_symbol", "modification_type", "functional_role", "organism")]
  rownames(result) <- NULL
  return(unique(result))
}

#' Get RNA Modification Proteins from Built-in Database
#' @noRd
.get_built_in_rna_mod_proteins <- function(modification_types, include_writers, include_readers, 
                                        include_erasers, organism, data_source) {
  
  # Read the built-in database
  protein_data <- system.file("extdata", "rna_modification_proteins.csv", package = "MPAGE")
  
  if (!file.exists(protein_data)) {
    stop("Built-in RNA modification protein database not found")
  }
  
  proteins <- utils::read.csv(protein_data, stringsAsFactors = FALSE)
  
  # Filter by organism
  if (!is.null(organism) && organism != "all") {
    proteins <- proteins[proteins$organism == organism, ]
  }
  
  # Filter by data source
  if (!is.null(data_source) && data_source != "all") {
    proteins <- proteins[proteins$evidence_source == data_source, ]
  }
  
  # Filter by modification types
  if (!is.null(modification_types)) {
    proteins <- proteins[proteins$modification_type %in% modification_types, ]
  }
  
  # Filter by functional roles
  roles_to_include <- c()
  if (include_writers) roles_to_include <- c(roles_to_include, "writer")
  if (include_readers) roles_to_include <- c(roles_to_include, "reader")
  if (include_erasers) roles_to_include <- c(roles_to_include, "eraser")
  
  if (length(roles_to_include) > 0) {
    proteins <- proteins[proteins$functional_role %in% roles_to_include, ]
  }
  
  # Return all columns from built-in database
  result <- proteins[, c("gene_symbol", "modification_type", "functional_role", 
                          "evidence_source", "pmid", "reference", "organism")]
  rownames(result) <- NULL
  return(unique(result))
}

#' Get RNA Modification References
#'
#' Retrieve literature references for RNA modification proteins from built-in database
#'
#' @param reference_filter Character vector of reference names to filter
#' @param organism Character string specifying organism (default: "Homo sapiens")
#' @param data_source Character string specifying data source filter (default: "Curated_database")
#'
#' @return Data frame with unique references and their PMIDs
#'
#' @export
#'
#' @examples
#' refs <- get_rna_mod_references()
#' print(refs)
get_rna_mod_references <- function(reference_filter = NULL, organism = "Homo sapiens", 
                                  data_source = "Curated_database") {
  proteins <- .get_built_in_rna_mod_proteins(
    modification_types = NULL,
    include_writers = TRUE,
    include_readers = TRUE,
    include_erasers = TRUE,
    organism = organism,
    data_source = data_source
  )
  
  refs <- unique(proteins[, c("reference", "pmid")])
  refs <- refs[order(refs$reference), ]
  
  if (!is.null(reference_filter)) {
    refs <- refs[refs$reference %in% reference_filter, ]
  }
  
  return(refs)
}

#' List Available Modification Types
#'
#' List all available RNA modification types in the built-in database
#'
#' @param organism Character string specifying organism (default: "Homo sapiens")
#' @param data_source Character string specifying data source filter (default: "Curated_database")
#'
#' @return Character vector of modification types
#'
#' @export
#'
#' @examples
#' mod_types <- list_modification_types()
#' print(mod_types)
list_modification_types <- function(organism = "Homo sapiens", 
                                   data_source = "Curated_database") {
  proteins <- .get_built_in_rna_mod_proteins(
    modification_types = NULL,
    include_writers = TRUE,
    include_readers = TRUE,
    include_erasers = TRUE,
    organism = organism,
    data_source = data_source
  )
  unique(proteins$modification_type)
}

#' List Available Organisms
#'
#' List all available organisms in the built-in database
#'
#' @param data_source Character string specifying data source filter (default: "Curated_database")
#'
#' @return Character vector of organisms
#'
#' @export
#'
#' @examples
#' organisms <- list_organisms()
#' print(organisms)
list_organisms <- function(data_source = "Curated_database") {
  protein_data <- system.file("extdata", "rna_modification_proteins.csv", package = "MPAGE")
  
  if (!file.exists(protein_data)) {
    return("Homo sapiens")
  }
  
  proteins <- utils::read.csv(protein_data, stringsAsFactors = FALSE)
  
  if (!is.null(data_source) && data_source != "all") {
    proteins <- proteins[proteins$evidence_source == data_source, ]
  }
  
  unique(proteins$organism)
}

#' Validate User Protein File
#'
#' Validate that a user-provided protein file has the correct format
#'
#' @param file_path Character string specifying path to the CSV file
#'
#' @return TRUE if file is valid, FALSE otherwise
#'
#' @export
#'
#' @examples
#' # Create a sample file
#' sample_data <- data.frame(
#'   gene_symbol = c("METTL3", "YTHDF1"),
#'   modification_type = c("m6A", "m6A"),
#'   functional_role = c("writer", "reader"),
#'   organism = c("Homo sapiens", "Homo sapiens")
#' )
#' write.csv(sample_data, "sample_proteins.csv", row.names = FALSE)
#' validate_protein_file("sample_proteins.csv")
validate_protein_file <- function(file_path) {
  if (!file.exists(file_path)) {
    stop("File not found: ", file_path)
  }
  
  tryCatch({
    proteins <- utils::read.csv(file_path, stringsAsFactors = FALSE)
    
    required_cols <- c("gene_symbol", "modification_type", "functional_role", "organism")
    missing_cols <- setdiff(required_cols, colnames(proteins))
    
    if (length(missing_cols) > 0) {
      warning("Missing required columns: ", paste(missing_cols, collapse = ", "))
      return(FALSE)
    }
    
    # Check data types
    if (!all(sapply(proteins[required_cols], is.character))) {
      warning("All columns should be character type")
      return(FALSE)
    }
    
    # Check for empty values
    if (any(is.na(proteins[required_cols]) | proteins[required_cols] == "")) {
      warning("Empty values found in required columns")
      return(FALSE)
    }
    
    return(TRUE)
  }, error = function(e) {
    warning("Error reading file: ", e$message)
    return(FALSE)
  })
}