# MPAGE: Module-based Pair-wise Analysis of Gene Expression

<!-- badges: start -->
[![R-CMD-check](https://github.com/loganylchen/MPAGE/workflows/R-CMD-check/badge.svg)](https://github.com/loganylchen/MPAGE/actions)
[![Codecov test coverage](https://codecov.io/gh/loganylchen/MPAGE/branch/main/graph/badge.svg)](https://app.codecov.io/gh/loganylchen/MPAGE?branch=main)
<!-- badges: end --

The MPAGE package provides a complete workflow for analyzing RNA modification modules and their dynamic changes across samples. It implements a module-based pair-wise analysis of gene expression, integrating PPI network construction, module identification, and various quantification methods for module activity and gene pair relationships.

## Installation

You can install the development version of MPAGE from [GitHub](https://github.com/) with:

``` r
# install.packages("devtools")
devtools::install_github("loganylchen/MPAGE")
```

## Quick Start

```r
library(MPAGE)

# Step 1: Get RNA modification proteins
rna_mod_proteins <- get_rna_mod_proteins(
  modification_types = c("m6A", "m5C", "Psi"),
  include_writers = TRUE,
  include_readers = TRUE,
  include_erasers = TRUE
)

# Step 2: Build PPI network
ppi_network <- build_ppi_network(
  proteins = rna_mod_proteins$gene_symbol,
  data_sources = c("STRING", "BioGRID"),
  min_confidence = 0.7
)

# Step 3: Identify modules
modules <- identify_modules(
  ppi_network = ppi_network,
  algorithms = c("FASTGREEDY", "MCODE"),
  min_module_size = 5
)

# Step 4: Filter RNA modules
rna_modules <- filter_rna_modules(
  modules,
  rna_proteins = rna_mod_proteins$gene_symbol,
  min_rna_proteins = 2
)

# Step 5: Annotate modules
annotated_modules <- annotate_modules(
  rna_modules,
  ppi_network = ppi_network,
  organism = "hsa"
)
```

## Documentation

For detailed documentation and tutorials, see the package [vignettes](https://loganylchen.github.io/MPAGE/articles/MPAGE.html).

## Contributing

Contributions are welcome! Please see [CONTRIBUTING.md](CONTRIBUTING.md) for guidelines.

## License

This package is released under the [MIT License](LICENSE).