# MPAGE 0.1.0 (2026-03-05)

## Initial Release

This is the first development release of MPAGE (Module-based Pair-wise Analysis of Gene Expression).

### Features
* `plot_modification_volcano()`: Create volcano plots and MA plots for RNA modification differential analysis
  * Support for both volcano and MA plot types
  * Customizable thresholds and colors
  * Automatic labeling of significant genes
  * Publication-quality output using ggplot2

* `plot_modification_heatmap()`: Visualize RNA modification patterns across samples
  * Interactive clustering options
  * Customizable color palettes
  * Support for complex experimental designs

* `consensus_modifications()`: Identify consensus RNA modifications across multiple samples
  * Support for multiple consensus algorithms
  * Flexible threshold settings
  * Integration with downstream analysis tools

### Infrastructure
* Set up package structure with roxygen2 documentation
* Add comprehensive test suite using testthat
* Configure GitHub Actions for CI/CD
* Add code coverage tracking with codecov

### Documentation
* Initial README with installation instructions and quick start guide
* Function-level documentation with examples
* Contributing guidelines

## Future Plans

### Upcoming Features (v0.2.0)
* `get_rna_mod_proteins()`: Retrieve RNA modification protein databases
* `build_ppi_network()`: Construct protein-protein interaction networks
* `identify_modules()`: Identify functional modules from PPI networks
* `filter_rna_modules()`: Filter modules for RNA modification relevance

### Known Issues
* See GitHub Issues for current limitations and planned improvements
