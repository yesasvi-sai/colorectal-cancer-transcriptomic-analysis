# Changelog

All notable changes to this project will be documented in this file.

The format is based on [Keep a Changelog](https://keepachangelog.com/en/1.0.0/),
and this project adheres to [Semantic Versioning](https://semver.org/spec/v2.0.0.html).

## [1.0.0] - 2026-02-19

### Added
- Initial release of CRC transcriptomics analysis
- Complete analysis pipeline for GSE33126 dataset
- Paired differential expression analysis using limma
- Pathway enrichment analysis (ORA and GSEA)
- Quality control visualizations
- Comprehensive documentation
- Example results and figures
- MIT License
- Contributing guidelines

### Features
- **Data Processing**:
  - Automated GEO data download
  - Log2 transformation and quantile normalization
  - Expression-based filtering
  - Quality control checks

- **Statistical Analysis**:
  - limma with paired design (patient blocking)
  - Empirical Bayes moderation
  - Benjamini-Hochberg FDR correction
  - Dual thresholds (adj. P < 0.05, |logFC| > 1)

- **Pathway Enrichment**:
  - Over-representation analysis (upregulated genes)
  - Over-representation analysis (downregulated genes)
  - Gene Set Enrichment Analysis (GSEA)
  - Reactome pathway database

- **Visualizations**:
  - Boxplots of normalized expression
  - Sample correlation heatmaps
  - PCA plots
  - Volcano plots
  - MA plots
  - Pathway enrichment dot plots
  - GSEA enrichment plots

- **Documentation**:
  - Comprehensive README with methods and results
  - Detailed methods documentation
  - Code comments throughout
  - Contributing guidelines
  - Example outputs

### Results
- 1,473 significantly dysregulated genes identified
- 678 upregulated, 795 downregulated
- Strong proliferation pathway enrichment (NES up to 2.74)
- Metabolic pathway suppression (p-values < 10⁻⁷)
- Validation of core CRC biology

## [Unreleased]

### Planned
- Additional visualization options
- Integration with clinical data
- Survival analysis module
- Subtype classification
- Comparison with other CRC datasets
- RNA-seq adaptation
- Unit tests
- Continuous integration setup

---

## Version History Summary

- **v1.0.0** (2026-02-19): Initial release with complete analysis pipeline
