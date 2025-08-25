#' Functional Enrichment Analysis
#'
#' Perform functional enrichment analysis using KEGG, GO-Biological Process, and WikiPathways
#'
#' @param gene_set Vector of gene symbols for enrichment analysis
#' @param key_proteins Vector of key proteins to check for presence in enriched pathways
#' @param organism Organism identifier (default: "hsa" for human)
#' @param p_value_cutoff P-value cutoff for significance (default: 0.05)
#' @param q_value_cutoff Q-value/FDR cutoff for significance (default: 0.2)
#' @param min_genes Minimum number of genes required in a pathway (default: 5)
#' @param max_genes Maximum number of genes allowed in a pathway (default: 500)
#'
#' @return List containing enrichment results for KEGG, GO-BP, and WikiPathways
#'   Each result includes pathway information, statistics, and key protein presence
#'
#' @export
#'
#' @examples
#' \dontrun{
#' # Basic enrichment analysis
#' genes <- c("TP53", "BRCA1", "BRCA2", "ATM", "CHEK2")
#' results <- functional_enrichment(gene_set = genes)
#'
#' # Enrichment with key proteins
#' key_prots <- c("TP53", "BRCA1")
#' results <- functional_enrichment(gene_set = genes, key_proteins = key_prots)
#'
#' # Custom parameters
#' results <- functional_enrichment(
#'   gene_set = genes,
#'   key_proteins = key_prots,
#'   p_value_cutoff = 0.01,
#'   min_genes = 10
#' )
#' }
functional_enrichment <- function(gene_set,
                                  key_proteins = NULL,
                                  organism = "hsa",
                                  p_value_cutoff = 0.05,
                                  q_value_cutoff = 0.2,
                                  min_genes = 5,
                                  max_genes = 500) {
  # Check required packages
  if (!requireNamespace("clusterProfiler", quietly = TRUE)) {
    stop("clusterProfiler package is required. Please install it with: install.packages('clusterProfiler')")
  }

  if (!requireNamespace("org.Hs.eg.db", quietly = TRUE)) {
    stop("org.Hs.eg.db package is required. Please install it with: BiocManager::install('org.Hs.eg.db')")
  }

  # Convert gene symbols to Entrez IDs
  gene_symbols <- toupper(unique(gene_set))

  # Map gene symbols to Entrez IDs
  org_db <- org.Hs.eg.db::org.Hs.eg.db
  entrez_ids <- AnnotationDbi::mapIds(
    org_db,
    keys = gene_symbols,
    keytype = "SYMBOL",
    column = "ENTREZID",
    multiVals = "first"
  )

  # Remove NAs and get unique IDs
  entrez_ids <- na.omit(unique(as.character(entrez_ids)))

  if (length(entrez_ids) == 0) {
    stop("No valid gene IDs found for enrichment analysis")
  }

  # Background universe
  universe <- AnnotationDbi::mapIds(
    org_db,
    keys = keys(org_db, keytype = "SYMBOL"),
    keytype = "SYMBOL",
    column = "ENTREZID",
    multiVals = "first"
  )
  universe <- na.omit(unique(as.character(universe)))

  # Prepare key proteins mapping
  key_protein_check <- !is.null(key_proteins)
  if (key_protein_check) {
    key_proteins_upper <- toupper(key_proteins)
  }

  # Function to check key protein presence in pathway
  check_key_proteins_in_pathway <- function(pathway_genes) {
    if (!key_protein_check) {
      return(list(has_key = FALSE, key_proteins = character(0)))
    }

    # Convert pathway genes to uppercase for matching
    pathway_genes_upper <- toupper(pathway_genes)

    # Find intersection with key proteins
    found_keys <- key_proteins_upper[key_proteins_upper %in% pathway_genes_upper]

    list(
      has_key = length(found_keys) > 0,
      key_proteins = found_keys
    )
  }

  # Function to run enrichment and process results
  run_enrichment <- function(OrgDb, gene_list, universe, database, p_cutoff, q_cutoff, min_g, max_g) {
    tryCatch(
      {
        enrich_result <- clusterProfiler::enrichGO(
          gene = gene_list,
          universe = universe,
          OrgDb = OrgDb,
          keyType = "ENTREZID",
          ont = database,
          pAdjustMethod = "BH",
          pvalueCutoff = p_cutoff,
          qvalueCutoff = q_cutoff,
          minGSSize = min_g,
          maxGSSize = max_g
        )

        if (is.null(enrich_result) || nrow(enrich_result@result) == 0) {
          return(NULL)
        }

        # Get detailed pathway information
        results <- as.data.frame(enrich_result)

        # Add key protein information
        if (key_protein_check) {
          key_info <- lapply(results$geneID, function(genes_str) {
            genes <- unlist(strsplit(genes_str, "/"))
            check_key_proteins_in_pathway(genes)
          })

          results$has_key_proteins <- sapply(key_info, function(x) x$has_key)
          results$key_proteins_found <- sapply(key_info, function(x) paste(x$key_proteins, collapse = ", "))
        } else {
          results$has_key_proteins <- FALSE
          results$key_proteins_found <- ""
        }

        return(results)
      },
      error = function(e) {
        warning(paste("Error in", database, "enrichment:", e$message))
        return(NULL)
      }
    )
  }

  # Function to run KEGG enrichment
  run_kegg_enrichment <- function(gene_list, organism, p_cutoff, q_cutoff, min_g, max_g) {
    tryCatch(
      {
        kegg_result <- clusterProfiler::enrichKEGG(
          gene = gene_list,
          organism = organism,
          pAdjustMethod = "BH",
          pvalueCutoff = p_cutoff,
          qvalueCutoff = q_cutoff,
          minGSSize = min_g,
          maxGSSize = max_g
        )

        if (is.null(kegg_result) || nrow(kegg_result@result) == 0) {
          return(NULL)
        }

        results <- as.data.frame(kegg_result)

        # Add key protein information
        if (key_protein_check) {
          key_info <- lapply(results$geneID, function(genes_str) {
            genes <- unlist(strsplit(genes_str, "/"))
            check_key_proteins_in_pathway(genes)
          })

          results$has_key_proteins <- sapply(key_info, function(x) x$has_key)
          results$key_proteins_found <- sapply(key_info, function(x) paste(x$key_proteins, collapse = ", "))
        } else {
          results$has_key_proteins <- FALSE
          results$key_proteins_found <- ""
        }

        return(results)
      },
      error = function(e) {
        warning(paste("Error in KEGG enrichment:", e$message))
        return(NULL)
      }
    )
  }

  # Function to run WikiPathways enrichment using enrichWP
  run_wikipathways_enrichment <- function(gene_list, organism, p_cutoff, q_cutoff, min_g, max_g) {
    tryCatch(
      {
        # Use enrichWP from clusterProfiler for WikiPathways
        if (organism == "hsa") {
          organism <- "Homo sapiens"
        } else {
          stop("only support hsa, the homo sapiens")
        }
        wp_result <- clusterProfiler::enrichWP(
          gene = gene_list,
          organism = organism,
          pvalueCutoff = p_cutoff,
          pAdjustMethod = "BH",
          qvalueCutoff = q_cutoff,
          minGSSize = min_g,
          maxGSSize = max_g
        )

        if (is.null(wp_result) || nrow(wp_result@result) == 0) {
          return(NULL)
        }

        results <- as.data.frame(wp_result)

        # Add key protein information
        if (key_protein_check) {
          key_info <- lapply(results$geneID, function(genes_str) {
            genes <- unlist(strsplit(genes_str, "/"))
            check_key_proteins_in_pathway(genes)
          })

          results$has_key_proteins <- sapply(key_info, function(x) x$has_key)
          results$key_proteins_found <- sapply(key_info, function(x) paste(x$key_proteins, collapse = ", "))
        } else {
          results$has_key_proteins <- FALSE
          results$key_proteins_found <- ""
        }

        return(results)
      },
      error = function(e) {
        warning(paste("Error in WikiPathways enrichment:", e$message))
        return(NULL)
      }
    )
  }

  # Run all enrichment analyses
  message("Running functional enrichment analysis...")

  # KEGG enrichment
  message("Running KEGG enrichment...")
  kegg_results <- run_kegg_enrichment(
    entrez_ids, organism, p_value_cutoff, q_value_cutoff, min_genes, max_genes
  )

  # GO-Biological Process enrichment
  message("Running GO-Biological Process enrichment...")
  go_bp_results <- run_enrichment(
    org_db, entrez_ids, universe, "BP", p_value_cutoff, q_value_cutoff, min_genes, max_genes
  )

  # WikiPathways enrichment
  message("Running WikiPathways enrichment...")
  wp_results <- run_wikipathways_enrichment(
    entrez_ids, organism, p_value_cutoff, q_value_cutoff, min_genes, max_genes
  )

  # Prepare summary
  summary_stats <- list(
    input_genes = length(gene_set),
    mapped_genes = length(entrez_ids),
    organism = organism,
    parameters = list(
      p_value_cutoff = p_value_cutoff,
      q_value_cutoff = q_value_cutoff,
      min_genes = min_genes,
      max_genes = max_genes
    ),
    results_count = list(
      kegg = ifelse(is.null(kegg_results), 0, nrow(kegg_results)),
      go_bp = ifelse(is.null(go_bp_results), 0, nrow(go_bp_results)),
      wikipathways = ifelse(is.null(wp_results), 0, nrow(wp_results))
    )
  )

  # Return comprehensive results
  results <- list(
    summary = summary_stats,
    kegg = kegg_results,
    go_bp = go_bp_results,
    wikipathways = wp_results,
    key_proteins = key_proteins,
    input_genes = gene_symbols,
    mapped_genes = entrez_ids
  )

  message("Enrichment analysis completed!")
  return(results)
}

#' Print enrichment results summary
#'
#' @param enrichment_results Results from functional_enrichment function
#' @param top_n Number of top results to display (default: 10)
#' @export
print_enrichment_summary <- function(enrichment_results, top_n = 10) {
  if (is.null(enrichment_results)) {
    cat("No enrichment results available.\n")
    return(invisible(NULL))
  }

  cat("=== Functional Enrichment Analysis Summary ===\n")
  cat(sprintf("Input genes: %d\n", enrichment_results$summary$input_genes))
  cat(sprintf("Mapped genes: %d\n", enrichment_results$summary$mapped_genes))
  cat(sprintf("Organism: %s\n", enrichment_results$summary$organism))
  cat("\n")

  databases <- c("kegg", "go_bp", "wikipathways")

  for (db in databases) {
    results <- enrichment_results[[db]]
    if (!is.null(results) && nrow(results) > 0) {
      cat(sprintf(
        "=== %s Enrichment (%d significant pathways) ===\n",
        toupper(db), nrow(results)
      ))

      # Show top results
      top_results <- head(results[order(results$pvalue), ], top_n)

      for (i in 1:nrow(top_results)) {
        row <- top_results[i, ]
        cat(sprintf(
          "%d. %s (p=%.2e, q=%.2e, genes=%d)",
          i, row$Description, row$pvalue, row$p.adjust, row$Count
        ))

        if (enrichment_results$summary$results_count[[db]] > 0 &&
          "has_key_proteins" %in% colnames(top_results) && row$has_key_proteins) {
          cat(sprintf(" [Key proteins: %s]", row$key_proteins_found))
        }
        cat("\n")
      }
      cat("\n")
    } else {
      cat(sprintf("=== %s Enrichment (No significant pathways) ===\n", toupper(db)))
    }
  }
}

#' Plot enrichment results
#'
#' @param enrichment_results Results from functional_enrichment function
#' @param database Database to plot ("kegg", "go_bp", or "wikipathways")
#' @param top_n Number of top pathways to display (default: 20)
#' @export
plot_enrichment_results <- function(enrichment_results, database = "kegg", top_n = 20) {
  if (!database %in% c("kegg", "go_bp", "wikipathways")) {
    stop("Database must be one of: kegg, go_bp, wikipathways")
  }

  results <- enrichment_results[[database]]
  if (is.null(results) || nrow(results) == 0) {
    warning(paste("No results available for", database))
    return(NULL)
  }

  if (!requireNamespace("ggplot2", quietly = TRUE)) {
    stop("ggplot2 package is required for plotting")
  }

  # Get top results
  top_results <- head(results[order(results$pvalue), ], top_n)
  top_results$Description <- factor(top_results$Description,
    levels = rev(unique(top_results$Description))
  )

  # Create dot plot
  p <- ggplot2::ggplot(top_results, ggplot2::aes(x = Count, y = Description)) +
    ggplot2::geom_point(ggplot2::aes(size = Count, color = p.adjust)) +
    ggplot2::scale_color_gradient(low = "blue", high = "red") +
    ggplot2::labs(
      title = paste(toupper(database), "Enrichment Results"),
      x = "Gene Count",
      y = "Pathway Description",
      color = "Adjusted P-value",
      size = "Gene Count"
    ) +
    ggplot2::theme_minimal() +
    ggplot2::theme(
      axis.text.y = ggplot2::element_text(size = 8),
      plot.title = ggplot2::element_text(hjust = 0.5)
    )

  return(p)
}
