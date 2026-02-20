# Quick Start Guide

Get up and running with the CRC transcriptomics analysis in 5 minutes.

## Prerequisites

- R version 4.3 or higher
- 8 GB RAM recommended
- Internet connection (for data download)

## Installation (5 minutes)

### Option 1: Using Conda (Recommended)

```bash
# 1. Clone repository
git clone https://github.com/yourusername/crc-transcriptomics-analysis.git
cd crc-transcriptomics-analysis

# 2. Create conda environment
conda env create -f environment.yml

# 3. Activate environment
conda activate crc-analysis

# 4. Launch R
R
```

### Option 2: Manual R Installation

Open R and run:

```r
# Install BiocManager
if (!require("BiocManager", quietly = TRUE))
    install.packages("BiocManager")

# Install Bioconductor packages
BiocManager::install(c("GEOquery", "limma", "ReactomePA", 
                       "enrichplot", "org.Hs.eg.db"))

# Install CRAN packages
install.packages(c("dplyr", "tidyr", "ggplot2", 
                   "pheatmap", "ggrepel", "stringr", "tibble"))
```

## Run Analysis (5 minutes)

### Complete Pipeline

```r
# Set working directory
setwd("path/to/crc-transcriptomics-analysis")

# Run complete analysis
source("scripts/main_analysis.R")
```

The script will:
1. ✅ Download GSE33126 from GEO (~2 min)
2. ✅ Perform QC and preprocessing (~1 min)
3. ✅ Run differential expression analysis (~30 sec)
4. ✅ Perform pathway enrichment (~2 min)
5. ✅ Generate all figures
6. ✅ Save results to CSV files

### Expected Output

```
results/
├── DEG_results.csv                # All differential expression results
├── ORA_upregulated.csv            # Enriched pathways (upregulated genes)
├── ORA_downregulated.csv          # Enriched pathways (downregulated genes)
└── GSEA_results.csv               # Gene set enrichment analysis

figures/
├── boxplot_normalized.png         # QC: normalized expression
├── correlation_heatmap.png        # QC: sample clustering
├── pca_plot.png                   # QC: principal components
├── volcano_plot.png               # DE: volcano plot
├── ma_plot.png                    # DE: MA plot
├── ora_upregulated.png            # Pathways: upregulated
├── ora_downregulated.png          # Pathways: downregulated (if generated)
└── gsea_*.png                     # GSEA: top pathways
```

## Examine Results

### In R

```r
# Load results
deg <- read.csv("results/DEG_results.csv")
ora_up <- read.csv("results/ORA_upregulated.csv")
gsea <- read.csv("results/GSEA_results.csv")

# Summary statistics
table(deg$Direction)
nrow(deg[deg$Significant, ])

# Top pathways
head(ora_up[, c("Description", "p.adjust", "geneID")], 10)
head(gsea[order(gsea$p.adjust), c("Description", "NES", "p.adjust")], 10)
```

### Key Findings

You should see:

1. **~1,473 significant genes** (adj. P < 0.05, |logFC| > 1)
   - ~678 upregulated
   - ~795 downregulated

2. **Upregulated pathways**:
   - rRNA processing (NES ≈ 2.74)
   - G1/S Transition (NES ≈ 2.61)
   - Cell cycle progression

3. **Downregulated pathways**:
   - Biological oxidations (p ≈ 10⁻⁷)
   - Phase II conjugation (p ≈ 10⁻⁴)
   - Metal homeostasis (p ≈ 10⁻⁸)

## Troubleshooting

### "Cannot download GSE33126"

**Problem**: Network timeout or GEO server unavailable

**Solution**:
```r
# Try with longer timeout
options(timeout = 300)
gse <- getGEO("GSE33126", GSEMatrix = TRUE)
```

### "Package not found"

**Problem**: Missing dependencies

**Solution**:
```r
# Check which packages are installed
installed.packages()[, c("Package", "Version")]

# Reinstall specific package
BiocManager::install("limma", force = TRUE)
```

### "Out of memory"

**Problem**: Insufficient RAM

**Solution**:
- Close other programs
- Use 64-bit R
- Increase R memory limit:
  ```r
  memory.limit(size = 8000)  # Windows only
  ```

### "No enriched pathways found"

**Problem**: Either too stringent cutoffs or analysis issue

**Check**:
```r
# How many significant genes?
sum(deg2$Significant)  # Should be ~1,473

# Try relaxing pathway cutoff
ora_up <- enrichPathway(gene = deg_up_ids, 
                        universe = universe_ids,
                        pvalueCutoff = 0.1)  # Less stringent
```

## Next Steps

### Explore Results

1. **Examine specific genes**:
   ```r
   # Top upregulated gene
   deg[deg$Direction == "Up", ][1, ]
   
   # Find specific gene
   deg[grep("TP53", deg$Symbol), ]
   ```

2. **Investigate pathways**:
   ```r
   # Genes in a pathway
   pathway_genes <- strsplit(ora_up$geneID[1], "/")[[1]]
   deg[deg$Symbol %in% pathway_genes, ]
   ```

3. **Custom visualizations**:
   ```r
   # Heatmap of top 50 genes
   top50 <- deg[deg$Significant, ][1:50, "probe"]
   pheatmap(expr_filt[top50, ], 
            annotation_col = sampleInfo,
            scale = "row")
   ```

### Modify Analysis

1. **Different significance thresholds**:
   ```r
   # More stringent
   deg_strict <- deg2 %>% 
     filter(adj.P.Val < 0.01, abs(logFC) > 2)
   ```

2. **Focus on specific pathways**:
   ```r
   # Only cell cycle pathways
   cc_pathways <- ora_up@result %>%
     filter(grepl("cell cycle|G1|G2|mitotic", Description, 
                  ignore.case = TRUE))
   ```

3. **Additional plots**:
   ```r
   # Expression of specific gene across samples
   gene_of_interest <- "ILMN_1704294"
   boxplot(expr_filt[gene_of_interest, ] ~ sampleInfo$group,
           main = gene_of_interest,
           ylab = "Expression")
   ```

## Learn More

- **Full documentation**: See [README.md](README.md)
- **Detailed methods**: See [docs/methods.md](docs/methods.md)
- **Scientific report**: See [docs/CRC_Analysis_Report.docx](docs/CRC_Analysis_Report.docx)
- **Contributing**: See [CONTRIBUTING.md](CONTRIBUTING.md)

## Get Help

- **Questions**: Open a GitHub Discussion
- **Bugs**: Create an Issue
- **Email**: your.email@example.com

---

**Estimated time to results**: 10 minutes
**Difficulty**: Beginner-friendly
**Data size**: ~50 MB download
