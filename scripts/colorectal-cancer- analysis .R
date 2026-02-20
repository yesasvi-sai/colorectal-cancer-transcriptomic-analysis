
suppressPackageStartupMessages({
  library(GEOquery)
  library(limma)
  library(dplyr)
  library(tidyr)
  library(stringr)
  library(tibble)
  library(ggplot2)
  library(ggrepel)
  library(pheatmap)
  
  library(ReactomePA)
  library(enrichplot)
  library(org.Hs.eg.db) 
})

## ---- a/b) Download GSE33126 + expression summary ----
my_id <- "GSE33126"
gse <- getGEO(my_id, GSEMatrix = TRUE)[[1]]

expr_raw <- exprs(gse)
summary(expr_raw)

## ---- c/d) Log-normalize + boxplot ----
expr <- log2(expr_raw)   # dataset is usually already normalized; assignment asks log-normalize
boxplot(expr, outline = FALSE, las = 2, main = "Log2 expression (all samples)")

## ---- e) Phenotype data (group + patient) ----
sampleInfo <- pData(gse) %>%
  dplyr::select(source_name_ch1, characteristics_ch1.1) %>%
  dplyr::rename(group = source_name_ch1,
                patient = characteristics_ch1.1)

# clean: "patient: 1" -> "1"
sampleInfo$patient <- gsub("patient: ", "", sampleInfo$patient)
sampleInfo$patient <- factor(sampleInfo$patient)

# ensure group labels match your design
# check unique(sampleInfo$group)
# set levels to "normal" then "tumor" (adjust if your GEO uses different spelling)
sampleInfo$group <- tolower(sampleInfo$group)
sampleInfo$group <- factor(sampleInfo$group, levels = c("normal", "tumour", "tumor"))
# If your dataset uses "Tumour" spelling:
levels(sampleInfo$group) <- gsub("tumour", "tumor", levels(sampleInfo$group))
sampleInfo$group <- factor(sampleInfo$group, levels = c("normal", "tumor"))

sampleInfo

## ---- f) Sample clustering heatmap (correlation) ----
corMatrix <- cor(expr, use = "pairwise.complete.obs")
pheatmap(corMatrix, main = "Sample–sample correlation")

## ---- g) Heatmap with annotations ----
rownames(sampleInfo) <- colnames(corMatrix)
pheatmap(corMatrix,
         annotation_col = sampleInfo,
         main = "Correlation heatmap (annotated)")

## ---- h) PCA scatter plot ----
pca <- prcomp(t(expr), scale. = FALSE)
pca_df <- cbind(sampleInfo, as.data.frame(pca$x))

ggplot(pca_df, aes(PC1, PC2, color = group, label = paste("Patient", patient))) +
  geom_point(size = 3) +
  geom_text_repel(max.overlaps = 50) +
  theme_minimal() +
  labs(title = "PCA of samples", x = "PC1", y = "PC2")

## ---- i) Differential expression (paired limma) ----
# Better normalization for microarray intensities (recommended)
expr_norm <- normalizeBetweenArrays(expr)

# simple filter (keep genes expressed above median in >2 samples)
cutoff <- median(expr_norm)
keep <- rowSums(expr_norm > cutoff) > 2
expr_filt <- expr_norm[keep, ]

# paired model: patient blocking + group effect
design <- model.matrix(~ patient + group, data = sampleInfo)
fit <- lmFit(expr_filt, design)
fit <- eBayes(fit)

# tumor vs normal coefficient
deg <- topTable(fit, coef = "grouptumor", number = Inf, sort.by = "P")

deg2 <- deg %>%
  tibble::rownames_to_column("probe") %>%
  mutate(Significant = (adj.P.Val < 0.05) & (abs(logFC) > 1),
         negLogAdjP = -log10(adj.P.Val))

head(deg2[, c("probe","logFC","P.Value","adj.P.Val")], 20)

## ---- j) Volcano plot (correct y-axis) ----
ggplot(deg2, aes(x = logFC, y = negLogAdjP, color = Significant)) +
  geom_point(alpha = 0.6) +
  theme_minimal() +
  labs(title = "Volcano: Tumor vs Normal (paired)",
       x = "log2 Fold Change", y = "-log10(adj P-value)")

## ---- k) MA plot (correct axes: AveExpr vs logFC) ----
ggplot(deg2, aes(x = AveExpr, y = logFC, color = Significant)) +
  geom_point(alpha = 0.6) +
  theme_minimal() +
  labs(title = "MA plot: Tumor vs Normal (paired)",
       x = "Average Expression", y = "log2 Fold Change")

## ---- l) Reactome ORA (Up/Down) using Entrez from annotation ----
anno <- as.data.frame(fData(gse)) %>% tibble::rownames_to_column("probe")
deg_annot <- deg2 %>% left_join(anno, by = "probe")

# helper: clean Entrez IDs
clean_entrez <- function(x) {
  x %>%
    as.character() %>%
    na.omit() %>%
    str_split("[; ,/]+") %>% unlist() %>%
    str_trim() %>%
    .[grepl("^[0-9]+$", .)] %>%
    unique()
}

# universe = all tested genes that have Entrez IDs
universe_ids <- clean_entrez(deg_annot$Entrez_Gene_ID)

# DEG hits (up/down)
deg_up_ids <- deg_annot %>% filter(Significant, logFC > 0) %>% pull(Entrez_Gene_ID) %>% clean_entrez()
deg_down_ids <- deg_annot %>% filter(Significant, logFC < 0) %>% pull(Entrez_Gene_ID) %>% clean_entrez()

ora_up <- enrichPathway(gene = deg_up_ids, universe = universe_ids,
                        pvalueCutoff = 0.05, pAdjustMethod = "BH",
                        readable = TRUE)

ora_down <- enrichPathway(gene = deg_down_ids, universe = universe_ids,
                          pvalueCutoff = 0.05, pAdjustMethod = "BH",
                          readable = TRUE)

# quick tables
ora_up@result %>% dplyr::select(Description, GeneRatio, p.adjust) %>% head(10)
ora_down@result %>% dplyr::select(Description, GeneRatio, p.adjust) %>% head(10)

# plots
dotplot(ora_up, showCategory = 10) + ggtitle("Reactome ORA (UP in Tumor)")
dotplot(ora_down, showCategory = 10) + ggtitle("Reactome ORA (DOWN in Tumor)")

## ---- m) Reactome GSEA (ranked by limma t-stat) ----
# Build ranked list across ALL tested genes
rank_df <- deg_annot %>%
  mutate(Entrez = as.character(Entrez_Gene_ID)) %>%
  filter(!is.na(Entrez)) %>%
  mutate(Entrez = str_split(Entrez, "[; ,/]+")) %>%
  tidyr::unnest(Entrez) %>%
  mutate(Entrez = str_trim(Entrez)) %>%
  filter(grepl("^[0-9]+$", Entrez)) %>%
  group_by(Entrez) %>%
  summarise(t = max(t, na.rm = TRUE), .groups = "drop")  # collapse duplicates

geneList <- rank_df$t
names(geneList) <- rank_df$Entrez
geneList <- sort(geneList, decreasing = TRUE)

gsea_reactome <- gsePathway(geneList = geneList,
                            pvalueCutoff = 0.05,
                            pAdjustMethod = "BH",
                            verbose = FALSE)

# top 10
gsea_reactome@result %>%
  dplyr::select(Description, NES, p.adjust) %>%
  dplyr::arrange(p.adjust, dplyr::desc(NES)) %>%
  head(10)

## ---- n) Example pathway visualization (optional) ----
# If you want a GSEA enrichment curve (recommended for report figures):
gseaplot2(gsea_reactome, geneSetID = "R-HSA-69206", title = "G1/S Transition")
gseaplot2(gsea_reactome, geneSetID = "R-HSA-397014", title = "Muscle contraction")

## ---- Save outputs (optional but nice) ----
write.csv(deg_annot, "Tumor_vs_Normal_paired_limma_DEG_annotated.csv", row.names = FALSE)
write.csv(as.data.frame(ora_up), "Reactome_ORA_UP.csv", row.names = FALSE)
write.csv(as.data.frame(ora_down), "Reactome_ORA_DOWN.csv", row.names = FALSE)
write.csv(as.data.frame(gsea_reactome), "Reactome_GSEA.csv", row.names = FALSE)

outdir <- "~/Desktop/GSE33126_outputs"
dir.create(outdir, showWarnings = FALSE, recursive = TRUE)

write.csv(deg_annot, file.path(outdir, "Tumor_vs_Normal_paired_limma_DEG_annotated.csv"), row.names = FALSE)
write.csv(as.data.frame(ora_up), file.path(outdir, "Reactome_ORA_UP.csv"), row.names = FALSE)
write.csv(as.data.frame(ora_down), file.path(outdir, "Reactome_ORA_DOWN.csv"), row.names = FALSE)
write.csv(as.data.frame(gsea_reactome), file.path(outdir, "Reactome_GSEA.csv"), row.names = FALSE)
