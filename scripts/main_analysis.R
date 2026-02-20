## =========================================================
## GSE33126 Colorectal Cancer Transcriptomic Analysis
## Complete Analysis Pipeline
## =========================================================
## 
## This script performs a comprehensive paired differential expression
## and pathway enrichment analysis of colorectal cancer vs normal tissues
##
## Author: [Your Name]
## Date: February 2026
## Dataset: GSE33126 (9 paired tumor-normal samples)
## =========================================================

## ---- Setup and Configuration ----

# Set random seed for reproducibility
set.seed(42)

# Create output directories if they don't exist
dir.create("data/raw", recursive = TRUE, showWarnings = FALSE)
dir.create("data/processed", recursive = TRUE, showWarnings = FALSE)
dir.create("results", recursive = TRUE, showWarnings = FALSE)
dir.create("figures/qc", recursive = TRUE, showWarnings = FALSE)
dir.create("figures/de", recursive = TRUE, showWarnings = FALSE)
dir.create("figures/pathways", recursive = TRUE, showWarnings = FALSE)

## ---- Load Required Libraries ----

cat("Loading required libraries...\n")

suppressPackageStartupMessages({
  # Core Bioconductor packages
  library(GEOquery)       # GEO data retrieval
  library(limma)          # Differential expression
  library(ReactomePA)     # Pathway enrichment
  library(enrichplot)     # Enrichment visualization
  library(org.Hs.eg.db)   # Gene annotation
  
  # Data manipulation
  library(dplyr)
  library(tidyr)
  library(stringr)
  library(tibble)
  
  # Visualization
  library(ggplot2)
  library(ggrepel)
  library(pheatmap)
})

cat("✓ All libraries loaded successfully\n\n")

## =========================================================
## STEP 1: DATA DOWNLOAD AND LOADING
## =========================================================

cat("=" , rep("=", 50), "=\n", sep="")
cat("STEP 1: Downloading GSE33126 from GEO\n")
cat("=" , rep("=", 50), "=\n\n", sep="")

# Download GEO dataset
my_id <- "GSE33126"
cat(paste0("Downloading ", my_id, "...\n"))

gse <- getGEO(my_id, GSEMatrix = TRUE, destdir = "data/raw")[[1]]

cat("✓ Dataset downloaded successfully\n")
cat(paste0("  - Platform: ", annotation(gse), "\n"))
cat(paste0("  - Samples: ", ncol(gse), "\n"))
cat(paste0("  - Features: ", nrow(gse), "\n\n"))

# Extract expression data
expr_raw <- exprs(gse)

cat("Expression data summary:\n")
print(summary(expr_raw[,1:3]))  # Show first 3 samples
cat("\n")

## =========================================================
## STEP 2: PREPROCESSING AND QUALITY CONTROL
## =========================================================

cat("=" , rep("=", 50), "=\n", sep="")
cat("STEP 2: Preprocessing and Quality Control\n")
cat("=" , rep("=", 50), "=\n\n", sep="")

## ---- Log2 transformation ----
cat("Applying log2 transformation...\n")
expr <- log2(expr_raw)

## ---- Normalization ----
cat("Normalizing between arrays...\n")
expr_norm <- normalizeBetweenArrays(expr, method = "quantile")
cat("✓ Quantile normalization complete\n\n")

## ---- Sample annotation ----
cat("Extracting sample annotations...\n")
sampleInfo <- pData(gse) %>%
  dplyr::select(source_name_ch1, characteristics_ch1.1) %>%
  dplyr::rename(group = source_name_ch1,
                patient = characteristics_ch1.1)

# Clean patient IDs
sampleInfo$patient <- gsub("patient: ", "", sampleInfo$patient)
sampleInfo$patient <- factor(sampleInfo$patient)

# Standardize group labels
sampleInfo$group <- tolower(sampleInfo$group)
sampleInfo$group <- gsub("tumour", "tumor", sampleInfo$group)
sampleInfo$group <- factor(sampleInfo$group, levels = c("normal", "tumor"))

cat("✓ Sample annotations:\n")
print(table(sampleInfo$group))
cat("\n")

## ---- Filtering ----
cat("Filtering low-expression probes...\n")
cutoff <- median(expr_norm)
keep <- rowSums(expr_norm > cutoff) > 2
expr_filt <- expr_norm[keep, ]

cat(paste0("✓ Retained ", nrow(expr_filt), " of ", nrow(expr_norm), 
           " probes (", round(100*nrow(expr_filt)/nrow(expr_norm), 1), "%)\n\n"))

## ---- QC Plots ----

cat("Generating quality control plots...\n")

# 1. Boxplot
png("figures/qc/boxplot_normalized.png", width=800, height=600, res=100)
boxplot(expr_norm, outline = FALSE, las = 2, 
        main = "Log2 Expression (Normalized)", 
        ylab = "Log2 Expression",
        cex.axis = 0.7)
dev.off()
cat("  ✓ Boxplot saved\n")

# 2. Correlation heatmap
corMatrix <- cor(expr_filt, use = "pairwise.complete.obs")
rownames(sampleInfo) <- colnames(corMatrix)

png("figures/qc/correlation_heatmap.png", width=800, height=800, res=100)
pheatmap(corMatrix, 
         annotation_col = sampleInfo,
         main = "Sample Correlation Heatmap",
         fontsize = 8)
dev.off()
cat("  ✓ Correlation heatmap saved\n")

# 3. PCA plot
pca <- prcomp(t(expr_filt), scale. = FALSE)
pca_df <- cbind(sampleInfo, as.data.frame(pca$x))
pca_var <- round(100 * summary(pca)$importance[2, 1:2], 1)

p_pca <- ggplot(pca_df, aes(PC1, PC2, color = group, label = paste("P", patient))) +
  geom_point(size = 4, alpha = 0.8) +
  geom_text_repel(size = 3, max.overlaps = 50) +
  theme_minimal(base_size = 12) +
  labs(title = "PCA: Tumor vs Normal Samples",
       x = paste0("PC1 (", pca_var[1], "% variance)"),
       y = paste0("PC2 (", pca_var[2], "% variance)"),
       color = "Tissue Type") +
  scale_color_manual(values = c("normal" = "#E74C3C", "tumor" = "#3498DB"))

ggsave("figures/qc/pca_plot.png", p_pca, width=8, height=6, dpi=300)
cat("  ✓ PCA plot saved\n\n")

## =========================================================
## STEP 3: DIFFERENTIAL EXPRESSION ANALYSIS
## =========================================================

cat("=" , rep("=", 50), "=\n", sep="")
cat("STEP 3: Differential Expression Analysis (Paired Design)\n")
cat("=" , rep("=", 50), "=\n\n", sep="")

## ---- Design matrix with patient blocking ----
cat("Creating design matrix with patient blocking...\n")
design <- model.matrix(~ patient + group, data = sampleInfo)
cat("✓ Design matrix created\n")
cat("  Model: Expression ~ Patient + Group\n")
cat(paste0("  Coefficients: ", ncol(design), "\n\n"))

## ---- Fit linear model ----
cat("Fitting linear model...\n")
fit <- lmFit(expr_filt, design)
fit <- eBayes(fit)
cat("✓ Model fitted with empirical Bayes moderation\n\n")

## ---- Extract results ----
cat("Extracting differential expression results...\n")
deg <- topTable(fit, coef = "grouptumor", number = Inf, sort.by = "P")

# Add probe IDs and calculate derived statistics
deg2 <- deg %>%
  tibble::rownames_to_column("probe") %>%
  mutate(
    Significant = (adj.P.Val < 0.05) & (abs(logFC) > 1),
    Direction = case_when(
      Significant & logFC > 0 ~ "Up",
      Significant & logFC < 0 ~ "Down",
      TRUE ~ "Not Sig"
    ),
    negLogAdjP = -log10(adj.P.Val)
  )

# Summary statistics
n_sig <- sum(deg2$Significant)
n_up <- sum(deg2$Direction == "Up")
n_down <- sum(deg2$Direction == "Down")

cat("✓ Results summary:\n")
cat(paste0("  - Total genes tested: ", nrow(deg2), "\n"))
cat(paste0("  - Significantly dysregulated: ", n_sig, " (", 
           round(100*n_sig/nrow(deg2), 1), "%)\n"))
cat(paste0("  - Upregulated: ", n_up, "\n"))
cat(paste0("  - Downregulated: ", n_down, "\n\n"))

cat("Top 10 upregulated genes:\n")
print(deg2 %>% filter(Direction == "Up") %>% 
        select(probe, logFC, P.Value, adj.P.Val) %>% head(10))
cat("\n")

cat("Top 10 downregulated genes:\n")
print(deg2 %>% filter(Direction == "Down") %>% 
        select(probe, logFC, P.Value, adj.P.Val) %>% head(10))
cat("\n")

## ---- DE Visualizations ----

cat("Generating differential expression plots...\n")

# Volcano plot
p_volcano <- ggplot(deg2, aes(x = logFC, y = negLogAdjP, color = Direction)) +
  geom_point(alpha = 0.5, size = 1) +
  scale_color_manual(values = c("Up" = "#E74C3C", "Down" = "#3498DB", 
                                  "Not Sig" = "grey70")) +
  theme_minimal(base_size = 12) +
  labs(title = "Volcano Plot: Tumor vs Normal (Paired Analysis)",
       x = "Log2 Fold Change",
       y = "-Log10(Adjusted P-value)",
       color = "Regulation") +
  geom_hline(yintercept = -log10(0.05), linetype = "dashed", alpha = 0.5) +
  geom_vline(xintercept = c(-1, 1), linetype = "dashed", alpha = 0.5)

ggsave("figures/de/volcano_plot.png", p_volcano, width=8, height=6, dpi=300)
cat("  ✓ Volcano plot saved\n")

# MA plot
p_ma <- ggplot(deg2, aes(x = AveExpr, y = logFC, color = Direction)) +
  geom_point(alpha = 0.5, size = 1) +
  scale_color_manual(values = c("Up" = "#E74C3C", "Down" = "#3498DB", 
                                  "Not Sig" = "grey70")) +
  theme_minimal(base_size = 12) +
  labs(title = "MA Plot: Tumor vs Normal (Paired Analysis)",
       x = "Average Expression (log2)",
       y = "Log2 Fold Change",
       color = "Regulation") +
  geom_hline(yintercept = c(-1, 0, 1), linetype = "dashed", alpha = 0.5)

ggsave("figures/de/ma_plot.png", p_ma, width=8, height=6, dpi=300)
cat("  ✓ MA plot saved\n\n")

## =========================================================
## STEP 4: PATHWAY ENRICHMENT ANALYSIS
## =========================================================

cat("=" , rep("=", 50), "=\n", sep="")
cat("STEP 4: Pathway Enrichment Analysis\n")
cat("=" , rep("=", 50), "=\n\n", sep="")

## ---- Prepare gene lists ----

cat("Preparing gene lists with Entrez IDs...\n")

# Get gene annotations
anno <- as.data.frame(fData(gse)) %>% 
  tibble::rownames_to_column("probe")

# Merge with DE results
deg_annot <- deg2 %>% left_join(anno, by = "probe")

# Function to clean Entrez IDs
clean_entrez <- function(x) {
  x %>%
    as.character() %>%
    na.omit() %>%
    str_split("[; ,/]+") %>% 
    unlist() %>%
    str_trim() %>%
    .[grepl("^[0-9]+$", .)] %>%
    unique()
}

# Extract gene lists
universe_ids <- clean_entrez(deg_annot$Entrez_Gene_ID)
deg_up_ids <- deg_annot %>% 
  filter(Significant, logFC > 0) %>% 
  pull(Entrez_Gene_ID) %>% 
  clean_entrez()
deg_down_ids <- deg_annot %>% 
  filter(Significant, logFC < 0) %>% 
  pull(Entrez_Gene_ID) %>% 
  clean_entrez()

cat(paste0("  ✓ Universe: ", length(universe_ids), " genes\n"))
cat(paste0("  ✓ Upregulated: ", length(deg_up_ids), " genes\n"))
cat(paste0("  ✓ Downregulated: ", length(deg_down_ids), " genes\n\n"))

## ---- Over-Representation Analysis ----

cat("Running over-representation analysis...\n")

# ORA for upregulated genes
cat("  - Analyzing upregulated genes...\n")
ora_up <- enrichPathway(
  gene = deg_up_ids,
  universe = universe_ids,
  pvalueCutoff = 0.05,
  pAdjustMethod = "BH",
  readable = TRUE,
  minGSSize = 5,
  maxGSSize = 500
)

if (!is.null(ora_up) && nrow(ora_up@result) > 0) {
  cat(paste0("    ✓ Found ", nrow(ora_up@result), " enriched pathways\n"))
  cat("    Top 5 pathways:\n")
  print(ora_up@result %>% 
          select(Description, GeneRatio, p.adjust) %>% 
          head(5))
} else {
  cat("    ⚠ No significant pathways found\n")
}
cat("\n")

# ORA for downregulated genes
cat("  - Analyzing downregulated genes...\n")
ora_down <- enrichPathway(
  gene = deg_down_ids,
  universe = universe_ids,
  pvalueCutoff = 0.05,
  pAdjustMethod = "BH",
  readable = TRUE,
  minGSSize = 5,
  maxGSSize = 500
)

if (!is.null(ora_down) && nrow(ora_down@result) > 0) {
  cat(paste0("    ✓ Found ", nrow(ora_down@result), " enriched pathways\n"))
  cat("    Top 5 pathways:\n")
  print(ora_down@result %>% 
          select(Description, GeneRatio, p.adjust) %>% 
          head(5))
} else {
  cat("    ⚠ No significant pathways found\n")
}
cat("\n")

## ---- ORA Visualizations ----

cat("Generating ORA plots...\n")

if (!is.null(ora_up) && nrow(ora_up@result) > 0) {
  p_ora_up <- dotplot(ora_up, showCategory = 10, font.size = 10) +
    ggtitle("Reactome ORA: Upregulated in Tumor")
  ggsave("figures/pathways/ora_upregulated.png", p_ora_up, 
         width=10, height=6, dpi=300)
  cat("  ✓ ORA upregulated plot saved\n")
}

if (!is.null(ora_down) && nrow(ora_down@result) > 0) {
  p_ora_down <- dotplot(ora_down, showCategory = 10, font.size = 10) +
    ggtitle("Reactome ORA: Downregulated in Tumor")
  ggsave("figures/pathways/ora_downregulated.png", p_ora_down, 
         width=10, height=6, dpi=300)
  cat("  ✓ ORA downregulated plot saved\n")
}
cat("\n")

## ---- Gene Set Enrichment Analysis ----

cat("Running Gene Set Enrichment Analysis (GSEA)...\n")

# Create ranked gene list
rank_df <- deg_annot %>%
  mutate(Entrez = as.character(Entrez_Gene_ID)) %>%
  filter(!is.na(Entrez)) %>%
  mutate(Entrez = str_split(Entrez, "[; ,/]+")) %>%
  tidyr::unnest(Entrez) %>%
  mutate(Entrez = str_trim(Entrez)) %>%
  filter(grepl("^[0-9]+$", Entrez)) %>%
  group_by(Entrez) %>%
  summarise(t = max(t, na.rm = TRUE), .groups = "drop")

geneList <- rank_df$t
names(geneList) <- rank_df$Entrez
geneList <- sort(geneList, decreasing = TRUE)

cat(paste0("  ✓ Ranked gene list: ", length(geneList), " genes\n"))

# Run GSEA
gsea_reactome <- gsePathway(
  geneList = geneList,
  pvalueCutoff = 0.05,
  pAdjustMethod = "BH",
  verbose = FALSE,
  minGSSize = 15,
  maxGSSize = 500
)

if (!is.null(gsea_reactome) && nrow(gsea_reactome@result) > 0) {
  cat(paste0("  ✓ Found ", nrow(gsea_reactome@result), " enriched pathways\n"))
  cat("  Top 10 pathways (by NES):\n")
  print(gsea_reactome@result %>% 
          select(Description, NES, p.adjust) %>%
          arrange(p.adjust, desc(abs(NES))) %>%
          head(10))
} else {
  cat("  ⚠ No significant pathways found\n")
}
cat("\n")

## ---- GSEA Visualizations ----

cat("Generating GSEA plots...\n")

if (!is.null(gsea_reactome) && nrow(gsea_reactome@result) > 0) {
  # Get top pathway IDs
  top_pathways <- gsea_reactome@result %>%
    arrange(p.adjust) %>%
    head(5) %>%
    pull(ID)
  
  # Plot first two pathways
  if (length(top_pathways) >= 1) {
    p_gsea1 <- gseaplot2(gsea_reactome, 
                         geneSetID = top_pathways[1],
                         title = gsea_reactome@result %>% 
                           filter(ID == top_pathways[1]) %>% 
                           pull(Description))
    ggsave("figures/pathways/gsea_top1.png", p_gsea1, 
           width=8, height=6, dpi=300)
    cat("  ✓ GSEA plot 1 saved\n")
  }
  
  if (length(top_pathways) >= 2) {
    p_gsea2 <- gseaplot2(gsea_reactome, 
                         geneSetID = top_pathways[2],
                         title = gsea_reactome@result %>% 
                           filter(ID == top_pathways[2]) %>% 
                           pull(Description))
    ggsave("figures/pathways/gsea_top2.png", p_gsea2, 
           width=8, height=6, dpi=300)
    cat("  ✓ GSEA plot 2 saved\n")
  }
}
cat("\n")

## =========================================================
## STEP 5: SAVE RESULTS
## =========================================================

cat("=" , rep("=", 50), "=\n", sep="")
cat("STEP 5: Saving Results\n")
cat("=" , rep("=", 50), "=\n\n", sep="")

# Save differential expression results
write.csv(deg_annot, "results/DEG_results.csv", row.names = FALSE)
cat("✓ Differential expression results saved\n")

# Save pathway enrichment results
if (!is.null(ora_up)) {
  write.csv(as.data.frame(ora_up), "results/ORA_upregulated.csv", 
            row.names = FALSE)
  cat("✓ ORA upregulated results saved\n")
}

if (!is.null(ora_down)) {
  write.csv(as.data.frame(ora_down), "results/ORA_downregulated.csv", 
            row.names = FALSE)
  cat("✓ ORA downregulated results saved\n")
}

if (!is.null(gsea_reactome)) {
  write.csv(as.data.frame(gsea_reactome), "results/GSEA_results.csv", 
            row.names = FALSE)
  cat("✓ GSEA results saved\n")
}

# Save processed data
saveRDS(list(
  expr_raw = expr_raw,
  expr_norm = expr_norm,
  expr_filt = expr_filt,
  sampleInfo = sampleInfo,
  deg_results = deg2
), "data/processed/processed_data.rds")
cat("✓ Processed data saved\n\n")

## =========================================================
## ANALYSIS COMPLETE
## =========================================================

cat("=" , rep("=", 50), "=\n", sep="")
cat("ANALYSIS COMPLETE!\n")
cat("=" , rep("=", 50), "=\n\n", sep="")

cat("Summary of findings:\n")
cat(paste0("  • ", n_sig, " genes significantly dysregulated\n"))
cat(paste0("  • ", n_up, " genes upregulated in tumor\n"))
cat(paste0("  • ", n_down, " genes downregulated in tumor\n"))

if (!is.null(ora_up)) {
  cat(paste0("  • ", nrow(ora_up@result), " pathways enriched (upregulated)\n"))
}
if (!is.null(ora_down)) {
  cat(paste0("  • ", nrow(ora_down@result), " pathways enriched (downregulated)\n"))
}
if (!is.null(gsea_reactome)) {
  cat(paste0("  • ", nrow(gsea_reactome@result), " pathways enriched (GSEA)\n"))
}

cat("\nOutput locations:\n")
cat("  • Figures: figures/\n")
cat("  • Results: results/\n")
cat("  • Processed data: data/processed/\n\n")

cat("Session info:\n")
print(sessionInfo())
