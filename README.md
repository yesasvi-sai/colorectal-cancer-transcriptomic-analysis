# Colorectal Cancer Transcriptomic Analysis

[![R Version](https://img.shields.io/badge/R-4.3%2B-blue)](https://www.r-project.org/)
[![License](https://img.shields.io/badge/license-MIT-green)](LICENSE)
[![DOI](https://img.shields.io/badge/DOI-GSE33126-orange)](https://www.ncbi.nlm.nih.gov/geo/query/acc.cgi?acc=GSE33126)

## Overview

This repository contains a comprehensive transcriptomic analysis of paired tumor and normal colon tissues from colorectal cancer (CRC) patients. Using the GSE33126 microarray dataset, we perform rigorous differential expression analysis with paired statistical design and pathway enrichment analysis to validate established CRC molecular signatures.

**Key Finding**: This analysis confirms that CRC shows robust, reproducible transcriptomic changes characterized by proliferation pathway upregulation (G1/S transition NES=2.61, rRNA processing NES=2.74) and metabolic pathway suppression (biological oxidations p=2.0×10⁻⁷), validating core cancer biology principles.

## Table of Contents

- [Background](#background)
- [Dataset](#dataset)
- [Methods](#methods)
- [Results](#results)
- [Repository Structure](#repository-structure)
- [Installation](#installation)
- [Usage](#usage)
- [Key Findings](#key-findings)
- [Limitations](#limitations)
- [Citation](#citation)
- [License](#license)
- [Contact](#contact)

## Background

Colorectal cancer is one of the best-characterized solid tumors at the molecular level. This project validates established CRC pathway alterations using:

- **Paired experimental design**: Tumor and normal tissues from the same 9 patients
- **Rigorous statistics**: limma with patient blocking to eliminate inter-individual variability  
- **Comprehensive pathway analysis**: Over-representation analysis (ORA) and Gene Set Enrichment Analysis (GSEA)

### Scientific Contribution

This is a **validation study** that confirms core CRC biology is:
- ✅ Robust across different cohorts
- ✅ Reproducible with appropriate statistical methods
- ✅ Consistent with established cancer hallmarks
- ✅ Detectable even with modest sample sizes when using paired designs

**This is NOT a discovery study** - we confirm known biology rather than identify novel mechanisms.

## Dataset

**Source**: [GSE33126](https://www.ncbi.nlm.nih.gov/geo/query/acc.cgi?acc=GSE33126) from NCBI Gene Expression Omnibus (GEO)

**Platform**: Illumina HumanHT-12 V3.0 expression beadchip

**Samples**: 
- 18 total samples (9 paired tumor-normal)
- 9 colorectal adenocarcinomas
- 9 matched normal colon mucosa (resected ≥20 cm from tumor)

**Reference**: Callari M, et al. Comparison of microarray platforms for measuring differential microRNA expression in paired normal/cancer colon tissues. PLoS One. 2012;7(9):e45105.

## Methods

### Data Processing Pipeline

1. **Quality Control**
   - Boxplots for normalization assessment
   - Sample-sample correlation with hierarchical clustering
   - Principal component analysis (PCA)

2. **Preprocessing**
   - Log₂ transformation
   - Quantile normalization (`limma::normalizeBetweenArrays`)
   - Filtering: retained probes expressed above median in >2 samples
   - **Result**: 18,535 probes retained from 48,803

3. **Differential Expression Analysis**
   - Method: `limma` with paired design
   - Model: `Expression ~ Patient + Group`
   - Patient blocking eliminates inter-individual variability
   - Empirical Bayes moderation for variance stabilization
   - **Significance criteria**: adj. P < 0.05 (Benjamini-Hochberg) AND |logFC| > 1

4. **Pathway Enrichment Analysis**
   - Database: Reactome pathways
   - Methods:
     - Over-representation analysis (ORA) - upregulated genes
     - Over-representation analysis (ORA) - downregulated genes
     - Gene Set Enrichment Analysis (GSEA) - all genes ranked by t-statistic
   - Background: All genes with valid Entrez IDs (n=15,407)

### Software & Packages

```r
# Core packages
library(GEOquery)      # v2.60.0 - GEO data retrieval
library(limma)         # v3.48.0 - Differential expression
library(ReactomePA)    # v1.36.0 - Pathway enrichment
library(enrichplot)    # Visualization
library(org.Hs.eg.db)  # Gene annotation

# Data manipulation & visualization
library(dplyr)
library(tidyr)
library(ggplot2)
library(pheatmap)
library(ggrepel)
```

**R version**: 4.3.0 or higher

## Results

### 1. Quality Control Results

| Metric | Result | Interpretation |
|--------|--------|----------------|
| Samples passing QC | 18/18 (100%) | All samples high quality |
| Median correlation (within tissue) | 0.96 | High technical reproducibility |
| Median correlation (tumor-normal) | 0.91 | Clear biological separation |
| PC1 variance explained | ~28% | Captures tumor-normal difference |

**Conclusion**: Data quality is excellent with clear tumor-normal separation.

### 2. Differential Expression Summary

| Category | Count | Percentage |
|----------|-------|------------|
| Total genes tested | 18,535 | 100% |
| Significantly dysregulated | 1,473 | 7.9% |
| Upregulated in tumor | 678 | 46.0% |
| Downregulated in tumor | 795 | 54.0% |

**Top 5 Most Significant Genes**:

| Probe ID | Gene Symbol | logFC | Adj. P-Value | Direction |
|----------|-------------|-------|--------------|-----------|
| ILMN_1704294 | - | 3.27 | 9.83×10⁻⁶ | UP |
| ILMN_1659990 | - | 1.26 | 2.61×10⁻⁵ | UP |
| ILMN_1652431 | - | -6.00 | 3.15×10⁻⁵ | DOWN |
| ILMN_1724686 | - | 3.63 | 8.36×10⁻⁵ | UP |
| ILMN_1813704 | - | 4.75 | 8.36×10⁻⁵ | UP |

### 3. Pathway Enrichment Results

#### Upregulated Pathways (Top 5)

| Pathway | Genes | P-value (adj) | Biological Significance |
|---------|-------|---------------|------------------------|
| rRNA processing in nucleus/cytosol | 13/200 | 0.028 | ⭐⭐⭐ Increased protein synthesis capacity |
| G1/S Transition | 10/200 | 0.028 | ⭐⭐⭐ Loss of cell cycle checkpoints |
| Response to amino acid deficiency | 9/200 | 0.028 | ⭐⭐ Metabolic stress response |
| IL-4 and IL-13 signaling | 9/200 | 0.028 | ⭐⭐ Immune modulation |
| Collagen degradation | 8/200 | 0.011 | ⭐⭐⭐ Invasive phenotype |

#### Downregulated Pathways (Top 5)

| Pathway | Genes | P-value (adj) | Biological Significance |
|---------|-------|---------------|------------------------|
| Metallothioneins bind metals | 8/273 | 1.15×10⁻⁸ | ⭐⭐⭐ Metal homeostasis loss |
| Response to metal ions | 8/273 | 6.04×10⁻⁸ | ⭐⭐⭐ Oxidative stress response |
| Biological oxidations | 24/273 | 2.00×10⁻⁷ | ⭐⭐⭐ Phase I metabolism loss |
| Phase II conjugation | 13/273 | 3.18×10⁻⁴ | ⭐⭐⭐ Detoxification suppression |
| Drug ADME | 12/273 | 5.81×10⁻⁴ | ⭐⭐ Altered pharmacokinetics |

#### GSEA Results (Top 5)

| Pathway | NES | P-value (adj) | Interpretation |
|---------|-----|---------------|----------------|
| rRNA processing | 2.74 | 3.7×10⁻⁹ | Extremely strong upregulation |
| G1/S Transition | 2.61 | 3.7×10⁻⁹ | Strong coordinated activation |
| Translation | 2.47 | 3.7×10⁻⁹ | Enhanced protein synthesis |
| S Phase | 2.45 | 3.7×10⁻⁹ | Accelerated cell cycle |
| Mitotic G1 | 2.53 | 3.7×10⁻⁹ | Checkpoint dysregulation |

**NES Interpretation**: NES > 2.0 indicates very strong coordinated pathway regulation

## Repository Structure

```
crc-transcriptomics-analysis/
│
├── README.md                          # This file
├── LICENSE                            # MIT License
├── environment.yml                    # Conda environment specification
├── requirements.txt                   # Python requirements (if any)
│
├── scripts/                           # Analysis scripts
│   ├── 01_download_data.R            # GEO data retrieval
│   ├── 02_preprocessing.R            # QC and normalization
│   ├── 03_differential_expression.R  # limma paired analysis
│   ├── 04_pathway_enrichment.R       # ORA and GSEA
│   └── main_analysis.R               # Complete workflow
│
├── data/                              # Data directory
│   ├── raw/                          # Raw GEO data (gitignored)
│   └── processed/                    # Normalized, filtered data
│
├── results/                           # Analysis outputs
│   ├── DEG_results.csv               # Differential expression table
│   ├── ORA_upregulated.csv           # ORA results (up genes)
│   ├── ORA_downregulated.csv         # ORA results (down genes)
│   └── GSEA_results.csv              # GSEA results
│
├── figures/                           # Publication-quality figures
│   ├── qc/                           # Quality control plots
│   │   ├── boxplot_normalized.png
│   │   ├── correlation_heatmap.png
│   │   └── pca_plot.png
│   ├── de/                           # Differential expression plots
│   │   ├── volcano_plot.png
│   │   └── ma_plot.png
│   └── pathways/                     # Pathway analysis plots
│       ├── ora_upregulated.png
│       ├── ora_downregulated.png
│       ├── gsea_g1s_transition.png
│       └── gsea_muscle_contraction.png
│
└── docs/                              # Documentation
    ├── methods.md                     # Detailed methods
    ├── results_interpretation.md      # Biological interpretation
    └── CRC_Analysis_Report.pdf        # Full scientific report
```

## Installation

### Option 1: Using Conda (Recommended)

```bash
# Clone repository
git clone https://github.com/yourusername/crc-transcriptomics-analysis.git
cd crc-transcriptomics-analysis

# Create conda environment
conda env create -f environment.yml
conda activate crc-analysis

# Launch R
R
```

### Option 2: Manual Installation in R

```r
# Install Bioconductor if needed
if (!require("BiocManager", quietly = TRUE))
    install.packages("BiocManager")

# Core Bioconductor packages
BiocManager::install(c(
    "GEOquery",
    "limma",
    "ReactomePA",
    "enrichplot",
    "org.Hs.eg.db"
))

# CRAN packages
install.packages(c(
    "dplyr",
    "tidyr",
    "ggplot2",
    "pheatmap",
    "ggrepel",
    "stringr",
    "tibble"
))
```

## Usage

### Quick Start

Run the complete analysis pipeline:

```r
# Set working directory
setwd("path/to/crc-transcriptomics-analysis")

# Run complete analysis
source("scripts/main_analysis.R")
```

### Step-by-Step Analysis

```r
# 1. Download and load data
source("scripts/01_download_data.R")

# 2. Quality control and preprocessing
source("scripts/02_preprocessing.R")

# 3. Differential expression analysis
source("scripts/03_differential_expression.R")

# 4. Pathway enrichment analysis
source("scripts/04_pathway_enrichment.R")
```

### Expected Runtime

| Step | Time | Memory |
|------|------|--------|
| Data download | ~2 min | <500 MB |
| Preprocessing | ~1 min | <1 GB |
| Differential expression | ~30 sec | <1 GB |
| Pathway enrichment | ~2 min | <2 GB |
| **Total** | **~6 min** | **<2 GB** |

*Tested on: Intel i7, 16GB RAM*

## Key Findings

### 1. Proliferation Signature is Dominant

- **13 ribosome biogenesis genes** upregulated (NES=2.74, p<3.7×10⁻⁹)
- **10 G1/S transition genes** upregulated (NES=2.61, p<3.7×10⁻⁹)
- Validates accelerated cell cycle and increased protein synthesis capacity

### 2. Metabolic Functions are Profoundly Suppressed

- **24 Phase I metabolism genes** downregulated (p=2.0×10⁻⁷)
- **13 Phase II conjugation genes** downregulated (p=3.2×10⁻⁴)
- **8 metallothionein genes** downregulated (p=1.2×10⁻⁸)
- Indicates loss of normal colon detoxification functions

### 3. Paired Design Provides Strong Statistical Power

- Despite small sample size (n=9), achieved p-values <10⁻⁸ for key pathways
- PCA shows tissue type dominates over patient variability
- Demonstrates efficiency of paired experimental designs

### 4. Findings are Highly Consistent with Prior Literature

- All major pathways identified match TCGA CRC results
- Validates reproducibility of core cancer biology
- Confirms these are robust, fundamental CRC features

## Limitations

### Sample Size
- **n=9 patients** is modest but appropriate for paired validation study
- Cannot stratify by molecular subtypes (CMS1-4)
- Limited power for detecting rare or subtype-specific changes

### Study Design
- **Cross-sectional** - cannot determine temporal sequence of events
- Cannot distinguish driver from passenger alterations
- Requires functional validation to establish causation

### Technology
- **Microarray** has lower sensitivity and dynamic range than RNA-seq
- Cannot detect novel transcripts or isoform-specific changes
- Limited to well-annotated genes

### Clinical Data
- **No outcome data** (survival, metastasis, treatment response)
- Cannot assess prognostic value of signatures
- Cannot correlate with clinical behavior

### Scope
- **Validation study** - confirms known biology, does not discover novel mechanisms
- Pathways identified are well-established CRC features
- Not identifying novel therapeutic targets

## Biological Interpretation

### The Proliferation-Metabolism Axis

This analysis reveals the fundamental trade-off in cancer biology:

**Resource Reallocation**: Cancer cells cannot simultaneously maximize proliferation and maintain specialized functions. They solve this by:

1. **Upregulating growth infrastructure**
   - Ribosome biogenesis (NES=2.74)
   - Cell cycle progression (NES=2.61)
   - Translation machinery (NES=2.47)

2. **Downregulating differentiated functions**
   - Detoxification enzymes (p<10⁻⁷)
   - Metal homeostasis (p<10⁻⁸)
   - Specialized metabolism (p<10⁻⁴)

### Clinical Implications

While this study confirms rather than discovers biology, the findings support:

1. **Therapeutic rationale** for CDK4/6 inhibitors targeting G1/S transition
2. **Biomarker potential** of ribosome signature for proliferation assessment
3. **Pharmacokinetic considerations** due to altered drug metabolism enzymes

However, these applications require:
- Larger validation cohorts
- Clinical outcome correlation
- Functional validation of specific targets

## Citation

If you use this analysis or code, please cite:

```bibtex
@article{callari2012comparison,
  title={Comparison of microarray platforms for measuring differential microRNA expression in paired normal/cancer colon tissues},
  author={Callari, Maurizio and Dugo, Matteo and Musella, Vera and others},
  journal={PLoS One},
  volume={7},
  number={9},
  pages={e45105},
  year={2012},
  publisher={Public Library of Science}
}
```

**Dataset**: GSE33126 - https://www.ncbi.nlm.nih.gov/geo/query/acc.cgi?acc=GSE33126

## Reproducibility

All analyses are fully reproducible:

- ✅ Complete code provided
- ✅ Public dataset (GEO)
- ✅ Specific package versions documented
- ✅ Random seeds set for permutation tests
- ✅ Step-by-step instructions

## Contributing

Contributions are welcome! Please:

1. Fork the repository
2. Create a feature branch (`git checkout -b feature/improvement`)
3. Commit changes (`git commit -am 'Add improvement'`)
4. Push to branch (`git push origin feature/improvement`)
5. Open a Pull Request

## License

This project is licensed under the MIT License - see [LICENSE](LICENSE) file for details.

## Acknowledgments

- Original dataset: Callari M, et al. (2012)
- Bioconductor core team for excellent R packages
- TCGA for establishing CRC molecular classification framework
- Reactome for curated pathway annotations

## Contact

**Project Maintainer**: [Your Name]
- GitHub: [@yourusername](https://github.com/yourusername)
- Email: your.email@example.com

**Issues**: Please report bugs or request features via [GitHub Issues](https://github.com/yourusername/crc-transcriptomics-analysis/issues)

---

**Keywords**: colorectal cancer, transcriptomics, differential expression, pathway enrichment, paired analysis, limma, Reactome, GSEA, microarray, CRC, validation study

**Last Updated**: February 2026
