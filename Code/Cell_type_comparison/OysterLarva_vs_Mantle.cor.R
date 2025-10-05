library(Seurat)
library(future)
library(tidyverse)
library(pheatmap)
library(philentropy)
library(RColorBrewer)
library(viridis)

# Increase memory limit for parallelization
options(future.globals.maxSize = 50 * 1024^3)

# Load larval and mantle single-cell datasets
load(file = "data/Cgigas/larva/OysterLarva.scRNAseq.rdata")
larva_seob <- seob
load(file = "data/Cgigas/mantle/OysterMantle.scRNAseq.rdata")
mantle_seob <- seob

# Identify genes shared between larval and mantle datasets
common_genes <- intersect(rownames(mantle_seob), rownames(larva_seob))

# Compute average expression per cell type for shared genes
mantle_avg <- AverageExpression(mantle_seob, assays = "SCT", features = common_genes, group.by = "cell_type")$SCT
larvae_avg <- AverageExpression(larva_seob, assays = "SCT", features = common_genes, group.by = "cell_type")$SCT
expr_mantle <- as.matrix(mantle_avg)
expr_larva  <- as.matrix(larvae_avg)

# Extract cell type names
ct_mantle <- colnames(expr_mantle)
ct_larva  <- colnames(expr_larva)

# Initialize empty correlation matrix
cor_matrix <- matrix(NA, nrow = length(ct_mantle), ncol = length(ct_larva))
rownames(cor_matrix) <- ct_mantle
colnames(cor_matrix) <- ct_larva

# Pearson correlation
# Compute Pearson correlation between each mantle and larval cell type
for (i in ct_mantle) {
  for (j in ct_larva) {
    cor_matrix[i, j] <- cor(expr_mantle[, i], expr_larva[, j], method = "pearson")
  }
}

# Define the plotting order for cell types
mantle_order <- rev(c(
  "Shell-secreting epithelial cells1", "Shell-secreting epithelial cells2",
  "Shell-secreting epithelial cells3", "Shell-secreting epithelial cells4",
  "Shell-secreting epithelial cells5", "Proliferative cells",
  "Neurons", "Hemocytes (blast-like cells)", "Hemocytes (hyalinocytes)",
  "Hemocytes (small granule cells)", "Hemocytes (vesicular cells)",
  "Myocytes", "Myocytes + Hemocytes", "Cilia",
  "Cilia + Hemocytes", "Unknown", "Vesicular connective tissue cells"
))

larva_order <- c(
  "Shell gland1", "Shell gland2", "Shell gland3", "Shell gland4",
  "Proliferative cells", "Proliferative cells + Myocytes", "Neurons",
  "Neurons + Myocytes", "Hemocytes", "Hemocytes + Myocytes",
  "Myocytes", "Cilia", "Cilia + Myocytes", "Unknown"
)

# Reshape correlation matrix into long format for ggplot
cor_df <- as_tibble(cor_matrix, rownames = "mantle_celltype") %>%
  pivot_longer(-mantle_celltype, names_to = "larva_celltype", values_to = "correlation")

# Set factor levels for ordered plotting
cor_df_pearson <- cor_df %>%
  mutate(
    mantle_celltype = factor(mantle_celltype, levels = mantle_order),
    larva_celltype = factor(larva_celltype, levels = larva_order)
  )

# Plot Pearson correlation heatmap
p_pearson <- ggplot(cor_df_pearson, aes(x = larva_celltype, y = mantle_celltype, fill = correlation)) +
  geom_tile(color = "grey", size = 0.25) +
  scale_fill_gradient2(
    low = "white", high = "#d73027",
    midpoint = 0.2, limits = c(0, 0.8), name = "Pearson\nCorrelation"
  ) +
  theme_minimal() +
  theme(
    axis.text.x = element_text(angle = 45, hjust = 1, size = 9, color = "black"),
    axis.text.y = element_text(size = 9, color = "black"),
    axis.title.x = element_text(color = "black"),
    axis.title.y = element_text(color = "black"),
    panel.grid = element_line(0),
    legend.title = element_text(size = 10, color = "black"),
    legend.text = element_text(size = 9, color = "black"),
    plot.title = element_text(color = "black", size = 12, face = "bold")
  ) +
  labs(title = "Cell Type Correlation: Mantle vs Larva",
       x = "Larval Cell Types",
       y = "Mantle Cell Types",
       fill = "Pearson\nCorrelation") +
  coord_fixed()

# Save figure and matrix results
ggsave(filename = "Cgi.larva_mantle.pearson.pdf", plot = p_pearson, width = 8, height = 6)
write.csv(cor_df_pearson, "Cgi.larva_mantle.pearson.csv", row.names = FALSE)
write.csv(cor_matrix, file = "Cgi.larva_mantle.pearson.csv", quote = FALSE)

# Spearman correlation
# Initialize empty matrix for Spearman correlation
cor_matrix_spearman <- matrix(NA, nrow = length(ct_mantle), ncol = length(ct_larva))
rownames(cor_matrix_spearman) <- ct_mantle
colnames(cor_matrix_spearman) <- ct_larva

# Compute Spearman correlation for all pairs
for (i in ct_mantle) {
  for (j in ct_larva) {
    cor_matrix_spearman[i, j] <- cor(expr_mantle[, i], expr_larva[, j], method = "spearman")
  }
}

# Reshape matrix for plotting
cor_df_spearman <- as_tibble(cor_matrix_spearman, rownames = "mantle_celltype") %>%
  pivot_longer(cols = -mantle_celltype, names_to = "larva_celltype", values_to = "correlation") %>%
  mutate(
    mantle_celltype = factor(mantle_celltype, levels = mantle_order),
    larva_celltype = factor(larva_celltype, levels = larva_order)
  )

# Plot Spearman correlation heatmap
p_spearman <- ggplot(cor_df_spearman, aes(x = larva_celltype, y = mantle_celltype, fill = correlation)) +
  geom_tile(color = "grey", size = 0.25) +
  scale_fill_gradient(low = "white", high = "#327bb1", limits = c(0.3, 0.8), name = "Spearman\nCorrelation") +
  theme_minimal() +
  theme(
    axis.text.x = element_text(angle = 45, hjust = 1, size = 9, color = "black"),
    axis.text.y = element_text(size = 9, color = "black"),
    axis.title.x = element_text(color = "black"),
    axis.title.y = element_text(color = "black"),
    panel.grid = element_line(0),
    legend.title = element_text(size = 10, color = "black"),
    legend.text = element_text(size = 9, color = "black"),
    plot.title = element_text(color = "black", size = 12, face = "bold")
  ) +
  labs(title = "Cell Type Correlation: Mantle vs Larva",
       x = "Larval Cell Types",
       y = "Mantle Cell Types",
       fill = "Spearman\nCorrelation") +
  coord_fixed()
ggsave(filename = "Cgi.larva_mantle.spearman.pdf", plot = p_spearman, width = 8, height = 6)

# Jensen–Shannon divergence (JSD)
# Convert expression matrix to probability distributions (column normalized)
expr_mantle_prob <- apply(expr_mantle, 2, function(x) x / sum(x))
expr_larva_prob <- apply(expr_larva, 2, function(x) x / sum(x))

# Initialize empty matrix for JSD similarity
jsd_matrix <- matrix(NA, nrow = length(ct_mantle), ncol = length(ct_larva))
rownames(jsd_matrix) <- ct_mantle
colnames(jsd_matrix) <- ct_larva

# Compute JSD for all pairs and convert to similarity (1 - JSD)
for (i in ct_mantle) {
  for (j in ct_larva) {
    P <- expr_mantle_prob[, i]
    Q <- expr_larva_prob[, j]
    mat <- rbind(P, Q)
    jsd_val <- JSD(mat, unit = "log2")
    jsd_matrix[i, j] <- 1 - jsd_val
  }
}

# Convert to long format for plotting
cor_df_jsd <- as_tibble(jsd_matrix, rownames = "mantle_celltype") %>%
  pivot_longer(cols = -mantle_celltype, names_to = "larva_celltype", values_to = "similarity") %>%
  mutate(
    mantle_celltype = factor(mantle_celltype, levels = mantle_order),
    larva_celltype = factor(larva_celltype, levels = larva_order)
  )

# Plot JSD-based similarity heatmap
p_jsd <- ggplot(cor_df_jsd, aes(x = larva_celltype, y = mantle_celltype, fill = similarity)) +
  geom_tile(color = "grey", size = 0.25) +
  scale_fill_gradient(low = "white", high = "#015b34", limits = c(0.3, 0.8), name = "JSD\nSimilarity") +
  theme_minimal() +
  theme(
    axis.text.x = element_text(angle = 45, hjust = 1, size = 9, color = "black"),
    axis.text.y = element_text(size = 9, color = "black"),
    axis.title.x = element_text(color = "black"),
    axis.title.y = element_text(color = "black"),
    legend.title = element_text(size = 10, color = "black"),
    legend.text = element_text(size = 9, color = "black"),
    plot.title = element_text(color = "black", size = 12, face = "bold")
  ) +
  labs(title = "Cell Type Similarity (JSD): Mantle vs Larva",
       x = "Larval Cell Types",
       y = "Mantle Cell Types") +
  coord_fixed()
ggsave(filename = "Cgi.larva_mantle.jsd.pdf", plot = p_jsd, width = 8, height = 6)

# Load cell-type-level marker gene tables
mantle_markers <- read.csv("data/Cgigas/mantle/Cg_cell_markers.csv")
larva_markers <- read.csv("data/Cgigas/larva/Cg_cell_markers.csv")

# Extract positively expressed marker genes per cluster
mantle_marker_sets <- mantle_markers %>%
  filter(avg_log2FC > 0) %>%
  group_by(cluster) %>%
  summarise(genes = list(gene), .groups = 'drop')

larva_marker_sets <- larva_markers %>%
  filter(avg_log2FC > 0) %>%
  group_by(cluster) %>%
  summarise(genes = list(gene), .groups = 'drop')

# Define Jaccard similarity function
jaccard <- function(set1, set2) {
  if (length(set1) == 0 | length(set2) == 0) return(0)
  inter <- length(intersect(set1, set2))
  union <- length(union(set1, set2))
  return(inter / union)
}

# Build Jaccard similarity matrix
ct_mantle <- mantle_marker_sets$cluster
ct_larva  <- larva_marker_sets$cluster

jaccard_matrix <- matrix(NA, nrow = length(ct_mantle), ncol = length(ct_larva))
rownames(jaccard_matrix) <- ct_mantle
colnames(jaccard_matrix) <- ct_larva

for (i in seq_along(ct_mantle)) {
  for (j in seq_along(ct_larva)) {
    set1 <- unlist(mantle_marker_sets$genes[i])
    set2 <- unlist(larva_marker_sets$genes[j])
    jaccard_matrix[i, j] <- jaccard(set1, set2)
  }
}

# Plot Jaccard index heatmap
jaccard_df <- as_tibble(jaccard_matrix, rownames = "mantle_celltype") %>%
  pivot_longer(cols = -mantle_celltype, names_to = "larva_celltype", values_to = "jaccard_index")

jaccard_df$mantle_celltype <- factor(jaccard_df$mantle_celltype, levels = mantle_order)
jaccard_df$larva_celltype <- factor(jaccard_df$larva_celltype, levels = larva_order)

# Plot jaccard heatmap
p_jaccard <- ggplot(jaccard_df, aes(x = larva_celltype, y = mantle_celltype, fill = jaccard_index)) +
  geom_tile(color = "grey", size = 0.25) +
  scale_fill_gradient(low = "white", high = "#db4003", limits = c(0, 0.25), name = "Jaccard\nIndex") +
  theme_minimal() +
  theme(
    axis.text.x = element_text(angle = 45, hjust = 1, size = 9, color = "black"),
    axis.text.y = element_text(size = 9, color = "black"),
    axis.title.x = element_text(color = "black"),
    axis.title.y = element_text(color = "black"),
    panel.grid = element_blank(),
    legend.title = element_text(size = 10, color = "black"),
    legend.text = element_text(size = 9, color = "black"),
    plot.title = element_text(color = "black", size = 12, face = "bold")
  ) +
  labs(title = "Cell Type Correlation: Mantle vs Larva",
       x = "Larval Cell Types",
       y = "Mantle Cell Types",
       fill = "Jaccard\nIndex") +
  coord_fixed()
ggsave(filename = "Cgi.larva_mantle.jaccard.pdf", plot = p_jaccard, width = 8, height = 6)

# Hypergeometric test for marker overlap
# Define optimized hypergeometric test function
test_overlap <- function(genes_A, genes_B, background_genes) {
  overlap <- length(intersect(genes_A, genes_B))
  m <- length(genes_A)
  n <- length(genes_B)
  N <- length(background_genes)
  if (overlap == 0) {
    p_val <- 1
    or <- NA
  } else {
    p_val <- phyper(overlap - 1, m, N - m, n, lower.tail = FALSE)
    or <- (overlap / n) / ((m - overlap) / (N - n))
  }
  return(list(overlap = overlap, p_val = p_val, or = or))
}

# Use all expressed genes as background
background_genes <- union(mantle_markers$gene, larva_markers$gene)

# Extract upregulated marker genes per cluster
mantle_marker_sets <- mantle_markers %>%
  filter(avg_log2FC > 0) %>%
  group_by(cluster) %>%
  summarise(genes = list(gene), .groups = "drop")

larva_marker_sets <- larva_markers %>%
  filter(avg_log2FC > 0) %>%
  group_by(cluster) %>%
  summarise(genes = list(gene), .groups = "drop")

ct_mantle <- mantle_marker_sets$cluster
ct_larva <- larva_marker_sets$cluster

# Initialize matrices for p-values and overlaps
p_list <- c()
p_matrix <- matrix(NA, nrow = length(ct_mantle), ncol = length(ct_larva),
                   dimnames = list(ct_mantle, ct_larva))
overlap_matrix <- matrix(NA, nrow = length(ct_mantle), ncol = length(ct_larva),
                         dimnames = list(ct_mantle, ct_larva))

# Perform hypergeometric tests for all cell type pairs
for (i in seq_along(ct_mantle)) {
  genes_A <- unlist(mantle_marker_sets$genes[i])
  for (j in seq_along(ct_larva)) {
    genes_B <- unlist(larva_marker_sets$genes[j])
    res <- test_overlap(genes_A, genes_B, background_genes)
    p_matrix[i, j] <- res$p_val
    overlap_matrix[i, j] <- res$overlap
    p_list <- c(p_list, res$p_val)
  }
}

# Multiple testing correction using Benjamini-Hochberg (BH) method
p_adj_list <- p.adjust(p_list, method = "BH")

# Convert the adjusted p-values back into a matrix (same dimensions as p_matrix)
p_adj_matrix <- matrix(p_adj_list, nrow = length(ct_mantle), ncol = length(ct_larva),
                       dimnames = list(ct_mantle, ct_larva), byrow = TRUE)

# Convert adjusted p-values to -log10 scale (retain Inf values)
logp_matrix <- -log10(p_adj_matrix)

# Cap maximum -log10 values for visualization
cap_value <- 20
logp_matrix_capped <- pmin(logp_matrix, cap_value)

# Define significance threshold (default: p < 0.01)
sig_threshold <- -log10(0.01)

# Convert matrix to long-format dataframe for ggplot2 visualization
logp_df <- as_tibble(logp_matrix_capped, rownames = "mantle_celltype") %>%
  pivot_longer(-mantle_celltype, names_to = "larva_celltype", values_to = "log10_adj_pvalue") %>%
  mutate(
    mantle_celltype = factor(mantle_celltype, levels = mantle_order),
    larva_celltype = factor(larva_celltype, levels = larva_order)
  )

# Plot heatmap of adjusted hypergeometric test results
p <- ggplot(logp_df, aes(x = larva_celltype, y = mantle_celltype, fill = log10_adj_pvalue)) +
  geom_tile(color = "grey", size = 0.25) +
  scale_fill_gradient(low = "white", high = "#e71f19", name = expression(-log[10]~adjusted~p)) +
  geom_text(
    data = logp_df %>% filter(log10_adj_pvalue > sig_threshold),
    aes(label = "**"),
    color = "black",
    size = 3
  ) +
  theme_minimal() +
  theme(
    axis.text.x = element_text(angle = 45, hjust = 1, size = 9, color = "black"),
    axis.text.y = element_text(size = 9, color = "black"),
    axis.title.x = element_text(color = "black"),
    axis.title.y = element_text(color = "black"),
    legend.title = element_text(size = 10, color = "black"),
    legend.text = element_text(size = 9, color = "black"),
    panel.grid = element_blank(),
    plot.title = element_text(color = "black", size = 12, face = "bold")
  ) +
  labs(title = "Adjusted Hypergeometric Overlap Significance",
       x = "Larval Cell Types",
       y = "Mantle Cell Types",
       fill = expression(-log[10](adjusted~p))) +
  coord_fixed()
ggsave("mantle_vs_larva_hypergeom_adjusted_logp.pdf", p, width = 8, height = 6)

# Export raw and adjusted p-value matrices and overlap counts
write.csv(p_matrix, file = "mantle_vs_larva_raw_pvalue_matrix.csv", quote = FALSE)
write.csv(p_adj_matrix, file = "mantle_vs_larva_adjusted_pvalue_matrix.csv", quote = FALSE)
write.csv(overlap_matrix, file = "mantle_vs_larva_overlap_matrix.csv", quote = FALSE)
