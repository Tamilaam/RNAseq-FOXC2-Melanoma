# === Load required libraries ===
library(tidyverse)
library(tximport)
library(biomaRt)
library(edgeR)
library(limma)
library(ggplot2)
library(clusterProfiler)
library(org.Mm.eg.db)

# === Set working directory ===
setwd("~/Desktop/ABioinformatics/assignment")

# === Prepare output folder ===
output_dir <- "last_sec_results"
dir.create(output_dir, showWarnings = FALSE)

# === Read sample metadata ===
samples <- read_csv("samples.csv")

# === Build paths to Kallisto quantifications ===
files <- file.path("data/kallisto_output", samples$sample, 
"abundance.tsv")
names(files) <- samples$sample
stopifnot(all(file.exists(files)))

# === Get transcript-to-gene mapping ===
mart <- useMart(biomart = "ensembl", dataset = "mmusculus_gene_ensembl")
tx2gene <- getBM(attributes = c("ensembl_transcript_id", 
"ensembl_gene_id"), mart = mart)
colnames(tx2gene) <- c("target_id", "gene_id")

# === Import transcript quantification ===
txi <- tximport(files,
                type = "kallisto",
                tx2gene = tx2gene,
                countsFromAbundance = "lengthScaledTPM",
                ignoreTxVersion = TRUE)

# === Create DGEList and normalize ===
myDGEList <- DGEList(counts = txi$counts)
myDGEList <- calcNormFactors(myDGEList)

# === Experimental design ===
samples$condition <- factor(samples$condition)
design <- model.matrix(~0 + samples$condition)
colnames(design) <- levels(samples$condition)

# === Voom transformation and plot ===
png(file.path(output_dir, "voom_plot.png"), width = 800, height = 600)
v <- voom(myDGEList, design, plot = TRUE)
dev.off()

# === Linear model fitting and contrasts ===
fit <- lmFit(v, design)
contrast <- makeContrasts(delta_vs_WT = delta - WT, levels = design)
fit2 <- contrasts.fit(fit, contrast)
fit2 <- eBayes(fit2)

# === Extract and save DE results ===
top_genes <- topTable(fit2, adjust.method = "BH", number = Inf, sort.by = 
"logFC") %>%
  rownames_to_column(var = "gene_id")
write.csv(top_genes, file.path(output_dir, "DE_results_full.csv"), 
row.names = FALSE)

# === Volcano plot ===
volcano <- top_genes %>%
  mutate(sig = ifelse(adj.P.Val < 0.05 & abs(logFC) > 1, "yes", "no"))
ggplot(volcano, aes(x = logFC, y = -log10(adj.P.Val), color = sig)) +
  geom_point(alpha = 0.6) +
  scale_color_manual(values = c("no" = "grey", "yes" = "red")) +
  theme_minimal() +
  labs(title = "Volcano Plot", x = "log2 Fold Change", y = "-log10 
Adjusted P-value")
ggsave(file.path(output_dir, "volcano_plot.png"))

# === PCA plot ===
pca <- prcomp(t(v$E))
pca_data <- data.frame(PC1 = pca$x[, 1], PC2 = pca$x[, 2], condition = 
samples$condition)
ggplot(pca_data, aes(x = PC1, y = PC2, color = condition)) +
  geom_point(size = 4) +
  theme_minimal() +
  labs(title = "PCA of Samples")
ggsave(file.path(output_dir, "pca_plot.png"))

# === Hierarchical clustering ===
dist_mat <- dist(t(v$E))
hc <- hclust(dist_mat, method = "complete")
png(file.path(output_dir, "hierarchical_clustering.png"), width = 800, 
height = 600)
plot(hc, labels = samples$sample, main = "Hierarchical Clustering 
Dendrogram")
dev.off()

# === Enrichment analysis ===
upregulated_genes <- top_genes %>%
  filter(adj.P.Val < 0.05 & logFC > 1) %>%
  pull(gene_id) %>%
  unique()
downregulated_genes <- top_genes %>%
  filter(adj.P.Val < 0.05 & logFC < -1) %>%
  pull(gene_id) %>%
  unique()

# GO for upregulated genes
ego_up <- enrichGO(gene = upregulated_genes,
                   OrgDb = org.Mm.eg.db,
                   keyType = "ENSEMBL",
                   ont = "BP",
                   pAdjustMethod = "BH",
                   pvalueCutoff = 0.05,
                   qvalueCutoff = 0.2,
                   readable = TRUE)
write.csv(ego_up@result, file.path(output_dir, "GO_upregulated.csv"))
barplot_up <- barplot(ego_up, showCategory = 10)
ggsave(file.path(output_dir, "go_enrichment_upregulated.png"), plot = 
barplot_up)

# GO for downregulated genes
ego_down <- enrichGO(gene = downregulated_genes,
                     OrgDb = org.Mm.eg.db,
                     keyType = "ENSEMBL",
                     ont = "BP",
                     pAdjustMethod = "BH",
                     pvalueCutoff = 0.05,
                     qvalueCutoff = 0.2,
                     readable = TRUE)
write.csv(ego_down@result, file.path(output_dir, "GO_downregulated.csv"))
barplot_down <- barplot(ego_down, showCategory = 10)
ggsave(file.path(output_dir, "go_enrichment_downregulated.png"), plot = 
barplot_down)

