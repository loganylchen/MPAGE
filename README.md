# mpage: Module-based Pair-wise Analysis of Gene Expression

<!-- badges: start -->
[![R-CMD-check](https://github.com/loganylchen/MPAGE/workflows/R-CMD-check/badge.svg)](https://github.com/loganylchen/MPAGE/actions)
[![Codecov test coverage](https://codecov.io/gh/loganylchen/MPAGE/branch/main/graph/badge.svg)](https://app.codecov.io/gh/loganylchen/MPAGE?branch=main)
[![License: MIT](https://img.shields.io/badge/License-MIT-yellow.svg)](https://opensource.org/licenses/MIT)
<!-- badges: end -->

`mpage` is the R reference implementation of **R3M** (RNA Modification
Module Model), an analytical framework that turns sparse RMP
(RNA modification protein) annotation into RMP-anchored, biologically
coherent network modules whose activity can be scored in any
whole-blood transcriptome.

It takes RMPs (writers / erasers / readers of m6A, m5C, m1A, psi, ...)
as anchors, integrates them into a high-confidence PPI network from
STRING / BioGRID / IntAct, recursively identifies RMP-anchored
modules under a size constraint, and provides downstream tooling for
module functional enrichment, three scales of pair-wise discriminative
features (intra-module gene pair / inter-module gene pair /
inter-module activity pair), and module-activity scoring
(GSVA, ssGSEA, z-score).

The companion paper, the full analysis pipeline, all figure-generation
scripts and the per-figure result objects live in the **mpage-paper**
repository:
[github.com/loganylchen/mpage-paper](https://github.com/loganylchen/mpage-paper).

## Table of contents

- [What `mpage` does](#what-mpage-does)
- [Installation](#installation)
- [One-shot pipeline: `run_mpage()`](#one-shot-pipeline-run_mpage)
- [Five-step pipeline (manual)](#five-step-pipeline-manual)
- [Output structure](#output-structure)
- [Reproducibility](#reproducibility)
- [Companion paper and Zenodo archive](#companion-paper-and-zenodo-archive)
- [Citation](#citation)
- [License](#license)
- [Contributing and issues](#contributing-and-issues)

## What `mpage` does

The R3M framework solves a specific problem: most RNA-modification
proteins (RMPs) are **sparsely annotated** in pathway databases. In
our curation only 12 of 85 RMPs are annotated across all three of
GO-BP, KEGG, and WikiPathways, which makes single-protein,
single-database enrichment uninformative for most of them. `mpage`
infers function from **network context** instead: each module is
anchored by one or more RMPs and inherits its functional
interpretation from its co-network-neighbour, expression-coordinated
gene set.

In one workflow `mpage` lets you:

1. **Curate RMPs** (built-in curation of 85 writers / erasers /
   readers across 8 modification types, or load your own CSV).
2. **Build per-database PPI networks** from STRING + BioGRID + IntAct
   and **merge** them into an integrated network.
3. **Recursively identify RMP-anchored modules** under a size
   constraint, anchored on a user-supplied target-protein list.
4. **Classify and enrich modules** by modification specificity
   (single vs cross-modification) and functional enrichment
   (GO / KEGG / WikiPathways).
5. **Score module activity per sample** (GSVA, ssGSEA, z-score) and
   build pair-wise discriminative features at three scales
   (intra-module gene pairs, inter-module gene pairs, inter-module
   activity pairs).

## Installation

`mpage` depends on several Bioconductor packages. Install those first,
then install `mpage` from GitHub:

```r
# 1. Bioconductor dependencies
if (!requireNamespace("BiocManager", quietly = TRUE))
  install.packages("BiocManager")
BiocManager::install(c(
  "STRINGdb", "clusterProfiler", "GSEABase", "GSVA",
  "Biobase", "org.Hs.eg.db", "AnnotationDbi"
))

# 2. mpage from GitHub
if (!requireNamespace("devtools", quietly = TRUE))
  install.packages("devtools")
devtools::install_github("loganylchen/MPAGE")
```

System requirements:

- R >= 4.0.0 (tested on R 4.3.3).
- `build_ppi_network()` reads pre-processed STRING / BioGRID / IntAct
  snapshots that ship inside the package at `inst/extdata/*.rds`, so
  no internet access is required to run it. To rebuild those
  snapshots from upstream you can use the helpers in `R/preprocessing_*.R`.

Optional but recommended for visualizations:

```r
install.packages(c("ggplot2", "dplyr", "tidyr", "pheatmap",
                   "RColorBrewer", "visNetwork"))
```

## One-shot pipeline: `run_mpage()`

If you already have a gene-expression matrix with case/control labels
and a set of modules (either the built-in 42 or your own), use
`run_mpage()`. It runs all three pair-wise scales (intra-module gene,
inter-module gene, inter-module activity) at once and returns a
single result object:

```r
library(mpage)

# expression_data : genes x samples (matrix or data.frame)
# sample_condition: named character vector with values "case" / "control"

results <- run_mpage(
  expression_data  = expr_mat,
  sample_condition = conditions,
  modules          = NULL,                       # NULL => built-in 42 RMP-anchored modules
  module_pairs     = NULL,                       # NULL => all inter-module pairs
  min_samples      = 3,
  activity_methods = c("gsva", "ssgsea", "zscore", "plage"),
  temp_dir         = "./temp"
)

names(results)
# -> "intra_module_results", "inter_module_results",
#    "module_activity_results", "summary", "metadata",
#    "main_analysis_metadata"
```

Use `run_mpage()` when you want the **whole pipeline in one call**.
Use the step-by-step pipeline below when you want to swap in a
different PPI network, change the module-identification parameters,
or stop after module discovery.

## Five-step pipeline (manual)

```r
library(mpage)

# Step 1. Curate RMPs (built-in curation of 85 RMPs, or your own CSV)
rmps <- get_rna_mod_proteins(use_built_in = TRUE)
# -> data.frame with columns gene_symbol / modification_type / functional_role
#
# Custom: get_rna_mod_proteins(use_built_in = FALSE,
#                              file_path = "my_proteins.csv")

# Step 2. Build per-source networks, then merge into the integrated PPI
networks <- build_ppi_network(
  data_sources  = c("STRING", "BIOGRID", "INTACT"),
  filters       = list(string = 700, intact = 0.7, biogrid = "physical"),
  versions      = list(string  = "12",
                       intact  = "2025-03-28",
                       biogrid = "4.4.248"),
  processed_dir = "./ppi_cache",
  species       = "9606"
)
# -> list of igraph objects keyed by source

ppi <- merge_ppi_networks(networks,
                          merge_method     = "union",
                          add_source_labels = TRUE)
# -> single igraph object representing the integrated PPI

# Step 3. Identify RMP-anchored modules
modules <- identify_modules_iterative(
  ppi_network     = ppi,
  min_module_size = 8,
  max_module_size = 95,
  target_proteins = rmps$gene_symbol
)
# -> list of modules; each has $module_id, $genes, $target_proteins, $size, ...

# Step 4. Classify modules and run functional enrichment
modules <- classify_modules(
  modules,
  classification = "modification_type",
  rna_proteins   = rmps
)

# functional_enrichment() runs on ONE gene set at a time (a character
# vector of gene symbols), so loop over your modules:
enrichment <- lapply(modules, function(m) {
  functional_enrichment(
    gene_set       = m$genes,
    key_proteins   = m$target_proteins,
    organism       = "hsa",
    p_value_cutoff = 0.05,
    q_value_cutoff = 0.2
  )
})
print_enrichment_summary(enrichment[[1]], top_n = 10)

# Step 5. Score module activity and compare case/control
activity <- score_modules(
  expression_data = expr_mat,
  modules         = modules,
  methods         = c("ssGSEA", "GSVA", "Zscore")
)
comparison <- compare_module_scores(
  module_scores = activity,
  sample_groups = conditions,
  method        = "ssGSEA"
)
```

Helper visualizations. `plot_network_summary()`,
`plot_module_features()`, `plot_module_scores()` and
`visualize_ppi_network()` write image files via their `save_path` or
`output_dir` argument; `plot_enrichment_results()` returns a `ggplot`
object that you save with `ggsave()`:

```r
plot_network_summary(ppi, save_path = "ppi_summary.png")        # writes images
plot_module_features(analyze_module_features(modules, ppi),
                     output_dir = "module_plots")               # writes images
plot_module_scores(comparison, output_dir = "score_plots")      # writes images
visualize_ppi_network(ppi, save_path = "ppi.png",
                      width = 8, height = 6)                    # writes image

g_enr <- plot_enrichment_results(enrichment[[1]],
                                  database = "kegg", top_n = 20) # ggplot
ggplot2::ggsave("enr.png", g_enr, width = 8, height = 6)
```

## Output structure

`run_mpage()` returns a list with the following slots:

| Slot | Description |
|---|---|
| `intra_module_results` | data.frame of intra-module gene-pair discriminative statistics (Fisher's exact test on the binary relative-ordering encoding, BH-FDR, selected pairs) |
| `inter_module_results` | data.frame of inter-module gene-pair statistics (same encoding and test) |
| `module_activity_results` | data.frame of per-method (GSVA / ssGSEA / z-score / PLAGE) module activity case-vs-control test results |
| `summary` | counts of significant pairs per scale and per method |
| `metadata` | run parameters carried through from `final_integrated_analysis()` |
| `main_analysis_metadata` | top-level metadata (timestamp, R version, package version, input parameters) |

The individual step functions return their own structured outputs
(see `?build_ppi_network`, `?identify_modules_iterative`,
`?functional_enrichment`, `?score_modules`). All exported functions
are documented; `?mpage` opens the package index.

## Reproducibility

- Randomized operations (cv.glmnet folds, random gene-set nulls) use
  `set.seed(422)` by default in the companion analysis pipeline; for
  your own runs set the seed at the start of your script.
- The companion **mpage-paper** repository runs the exact pipeline
  reported in the paper inside a pinned Docker image
  (`btrspg/env-mpage-paper:20260302`); see the
  [mpage-paper README](https://github.com/loganylchen/mpage-paper) for
  the per-figure entry points.

## Companion paper and Zenodo archive

| Resource | URL |
|---|---|
| Manuscript source (`paper/manuscript.md`) | [github.com/loganylchen/mpage-manuscript](https://github.com/loganylchen/mpage-manuscript) |
| Full analysis pipeline + figure generators | [github.com/loganylchen/mpage-paper](https://github.com/loganylchen/mpage-paper) |
| `mpage` R package source (this repo) | [github.com/loganylchen/MPAGE](https://github.com/loganylchen/MPAGE) |
| `mpage` R package frozen release on Zenodo | DOI: 10.5281/zenodo.XXXXXXX *(populated after the next tagged release)* |
| Paper analysis bundle on Zenodo | DOI: 10.5281/zenodo.XXXXXXX *(populated alongside the paper release)* |

## Citation

If you use `mpage` in published work, please cite both the paper and
the package:

```
Chen, L. et al. R3M: A Network-Prior Framework that Turns Sparse RMP
Annotation into a Platform-Robust Module-Activity Coordinate System.
[journal/year TBD]. DOI: <paper DOI>.

Chen, L. (2026). mpage: Module-based Pair-wise Analysis of Gene
Expression (Version 1.0.0) [Computer software].
DOI: 10.5281/zenodo.XXXXXXX.
```

After installation you can also retrieve a machine-readable citation:

```r
citation("mpage")
```

## License

`mpage` is released under the [MIT License](LICENSE).

## Contributing and issues

Bug reports and feature requests are very welcome at
[github.com/loganylchen/MPAGE/issues](https://github.com/loganylchen/MPAGE/issues).
Pull-request guidelines live in
[CONTRIBUTING.md](CONTRIBUTING.md).
