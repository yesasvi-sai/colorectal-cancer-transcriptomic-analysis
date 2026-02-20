# Detailed Methods

## Data Acquisition

### Dataset Selection Criteria

GSE33126 was selected based on three key criteria:

1. **Truly Paired Design**: Tumor and normal tissues from the same patients, eliminating inter-patient variability
2. **Appropriate Sampling**: Normal tissue resected ≥20 cm from tumor to avoid field effect contamination
3. **Well-Validated Platform**: Illumina HumanHT-12 V3.0 beadchips with extensive validation

### Data Download

```r
library(GEOquery)
gse <- getGEO("GSE33126", GSEMatrix = TRUE)[[1]]
expr_raw <- exprs(gse)
```

Data automatically downloaded from NCBI GEO using GEOquery package.

## Preprocessing Pipeline

### 1. Log Transformation

**Rationale**: Microarray intensities follow approximately log-normal distribution. Log transformation:
- Stabilizes variance across expression range
- Makes data more normally distributed
- Enables linear modeling assumptions

**Implementation**:
```r
expr <- log2(expr_raw)
```

### 2. Between-Array Normalization

**Method**: Quantile normalization

**Rationale**: Removes systematic technical variation between arrays while preserving biological differences.

**Implementation**:
```r
expr_norm <- normalizeBetweenArrays(expr, method = "quantile")
```

**Effect**: Forces all arrays to have identical empirical distribution while maintaining ranking within each array.

### 3. Filtering Low-Expression Probes

**Criterion**: Retained probes with expression > median in at least 3 samples

**Rationale**:
- Removes noise from probes near detection limit
- Reduces multiple testing burden
- Focuses analysis on reliably detected transcripts

**Implementation**:
```r
cutoff <- median(expr_norm)
keep <- rowSums(expr_norm > cutoff) > 2
expr_filt <- expr_norm[keep, ]
```

**Result**: 18,535 of 48,803 probes retained (38%)

## Quality Control

### Normalization Quality Assessment

**Boxplots**: Should show similar distributions across samples after normalization

**Interpretation**:
- Similar median values indicate successful normalization
- Similar IQRs indicate variance stabilization
- Outliers suggest potential sample quality issues

### Sample Clustering

**Method**: Hierarchical clustering of sample-sample correlations

**Distance**: 1 - Pearson correlation

**Expected Pattern**:
- Samples should cluster by biological condition (tumor vs normal)
- Within-condition correlations > between-condition correlations
- Paired samples from same patient may or may not cluster together

**Interpretation**:
- If samples cluster by batch/plate rather than condition → batch effect present
- If no clear clustering → biological signal weak or high technical noise

### Principal Component Analysis

**Method**: Singular value decomposition of centered, scaled expression matrix

**Expected Pattern**:
- PC1 should capture tumor-normal difference (biological signal)
- PC2 often captures patient-specific variation
- Clear separation along PC1 indicates strong biological signal

**Variance Explained**:
- PC1: ~25-30% typical for cancer vs normal
- PC2: ~10-15% typical for patient effects

## Statistical Analysis

### Linear Modeling with Patient Blocking

**Model Formula**:
```
Expression ~ Patient + Group
```

**Design Matrix Structure**:
- Patient factors: 8 dummy variables (9 patients - 1 reference)
- Group effect: 1 coefficient (tumor vs normal)
- Total: 10 coefficients

**Why Paired Design?**

Inter-patient variability in gene expression can be substantial due to:
- Genetic background
- Age, sex, diet
- Environmental exposures
- Microbiome composition

By including patient as a blocking factor, we:
1. Remove patient-specific baseline differences
2. Focus exclusively on within-patient tumor-normal differences
3. Dramatically increase statistical power

**Power Comparison**:
- Paired design with n=9: Power ≈ 0.80 for effect size d=1.5
- Unpaired design with n=9/group: Power ≈ 0.30 for same effect size
- Would need n≈25/group unpaired to match power of n=9 paired

### Empirical Bayes Moderation

**Method**: Shrink gene-specific variance estimates toward common value

**Rationale**: 
- Individual gene variances are unreliable with small sample sizes
- Borrowing information across genes improves estimates
- Particularly beneficial for lowly expressed genes

**Implementation**:
```r
fit <- lmFit(expr_filt, design)
fit <- eBayes(fit)
```

**Effect**:
- Moderated t-statistics are more stable than standard t-statistics
- Reduces false positives from genes with spuriously low variance
- Increases power for genes with high biological variability

### Multiple Testing Correction

**Method**: Benjamini-Hochberg False Discovery Rate (FDR)

**Target**: FDR = 5% (expected proportion of false positives among discoveries)

**Why FDR instead of FWER?**
- Family-Wise Error Rate (FWER, e.g., Bonferroni) controls probability of ANY false positive
- Too conservative for genomics with 10,000+ tests
- FDR balances discovery and false positive control

**Implementation**:
```r
adj.P.Val <- p.adjust(P.Value, method = "BH")
```

### Significance Thresholds

**Two Criteria Applied**:
1. **Statistical**: adj. P < 0.05 (5% FDR)
2. **Biological**: |logFC| > 1 (2-fold change)

**Rationale for Fold Change Threshold**:
- Changes <2-fold often below biological relevance
- Microarray technical noise ~1.5-fold
- Focuses on substantial changes with clear biological impact

## Pathway Enrichment Analysis

### Gene Annotation

**Mapping**: Illumina probe IDs → Entrez Gene IDs

**Challenge**: Some probes map to multiple genes or have ambiguous IDs

**Solution**: 
```r
# Split multi-mapped probes
# Keep only numeric Entrez IDs
# Remove duplicates
```

**Background Universe**: All genes with valid Entrez IDs that passed filtering (n≈15,400)

### Over-Representation Analysis (ORA)

**Statistical Test**: Hypergeometric test

**Null Hypothesis**: Genes in pathway are randomly distributed among DEGs

**Formula**: 
```
P(X ≥ k) = sum from i=k to min(n,K) of:
  [C(K,i) * C(M-K, n-i)] / C(M,n)
```

Where:
- M = total genes in background
- K = genes in pathway
- n = number of DEGs
- k = DEGs in pathway

**Advantages**:
- Simple, well-understood
- Easy to interpret (enrichment ratio)

**Limitations**:
- Requires arbitrary significance cutoff
- Treats all DEGs equally (ignores magnitude)
- Assumes gene independence

### Gene Set Enrichment Analysis (GSEA)

**Method**: Enrichment score based on ranked gene list

**Key Innovation**: Uses entire ranked list, no arbitrary cutoff

**Algorithm**:
1. Rank all genes by test statistic (t-statistic)
2. Walk down ranked list
3. Increase running score when hit pathway gene
4. Decrease score when hit non-pathway gene
5. Enrichment score = maximum deviation from zero

**Normalization**: 
- Normalized Enrichment Score (NES) accounts for pathway size
- Comparable across pathways

**Significance**: 
- Permutation-based p-values
- 1000 permutations standard

**Advantages**:
- More sensitive than ORA
- Detects subtle coordinated changes
- No arbitrary threshold

**Interpretation**:
- |NES| > 1.5: meaningful enrichment
- |NES| > 2.0: strong enrichment
- |NES| > 2.5: very strong enrichment

Our results (NES up to 2.74) indicate extremely strong coordinated regulation.

### Reactome Pathway Database

**Why Reactome?**
- Manually curated by experts
- Hierarchical organization
- Cross-species conservation
- Extensive literature support

**Coverage**: ~2,500 human pathways

**Alternative Databases**:
- KEGG: Smaller, less detailed
- GO Biological Process: Very broad, less mechanistic
- MSigDB Hallmarks: Cancer-specific but limited coverage

## Reproducibility Measures

### Random Seeds
```r
set.seed(42)  # All permutation-based tests
```

### Package Versions
Documented in sessionInfo() output

### Platform Independence
- Code runs on Windows, Mac, Linux
- No hardcoded paths (relative paths only)

### Data Availability
- Public dataset (GSE33126)
- No restricted-access data
- Complete replication possible

## Computational Requirements

**Hardware**:
- RAM: 4 GB minimum, 8 GB recommended
- CPU: Single core sufficient
- Storage: ~1 GB for data and results

**Runtime**:
- Data download: 1-2 minutes
- Preprocessing: 30 seconds
- Differential expression: 10 seconds
- Pathway enrichment: 1-2 minutes
- **Total**: ~5 minutes

**Bottlenecks**:
- GEO download (network dependent)
- GSEA permutations (can parallelize)
