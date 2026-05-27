# MPAGE News

## Version 1.0.0 (2026-03-02)

### New Features

* **Initial Release**: First stable version of MPAGE package
* **RNA Modification Protein Database**: Curated database of RNA modification-related proteins (writers, readers, erasers) with support for m6A, m5C, and pseudouridine modifications
* **PPI Network Construction**: Build protein-protein interaction networks from multiple data sources (STRING, BioGRID, IntAct)
* **Module Identification**: Identify functional modules from PPI networks using graph clustering algorithms (Louvain, Fast Greedy)
* **Gene Pair Analysis**: 
  - Intra-module gene pair analysis
  - Inter-module gene pair analysis  
  - Module activity comparison between conditions
* **Scoring Methods**: Multiple module scoring methods (GSVA, ssGSEA, Z-score, PLAGE)
* **Visualization Tools**: 
  - PPI network visualization (static and interactive)
  - Module score comparison plots
  - Enrichment result visualization
* **Functional Enrichment**: GO, KEGG, and WikiPathways enrichment analysis with clusterProfiler

### Documentation

* Comprehensive vignette with workflow examples
* Function-level documentation with examples
* AGENTS.md with coding guidelines

### Bug Fixes

None (initial release)

## Future Plans

* Support for additional RNA modification types
* Integration with more PPI databases
* Enhanced visualization options
* Parallel processing support for large datasets

---

For bug reports and feature requests, please visit:
https://github.com/loganylchen/MPAGE/issues
