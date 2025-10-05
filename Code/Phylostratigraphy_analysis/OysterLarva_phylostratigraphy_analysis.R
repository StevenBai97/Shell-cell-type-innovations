# Load libraries
library(dplyr)
library(ggplot2)
library(tidyr)
library(readr)
library(pheatmap)
library(myTAI)
library(tibble)

# Define custom colors for cell types
my_colors <- c(
  "Shell gland1" = "#c5362b",
  "Shell gland2" = "#e55442",
  "Shell gland3" = "#f08a74",
  "Shell gland4" = "#f5c7b9",
  "Proliferative cells" = "#f6c461",
  "Neurons" = "#ccc45d",
  "Hemocytes" = "#36729d",
  "Cilia" = "#9dbb65",
  "Myocytes" = "#72a87a",
  "Hemocytes + Myocytes" = "#5d998a",
  "Proliferative cells + Myocytes" = "#529192",
  "Neurons + Myocytes" = "#4b899b",
  "Cilia + Myocytes" = "#6cbacc",
  "Unknown" = "#add3d0"
)

# Load cell type average expression data
load(file = "Cg_average.rdata")

# Load phylostratum (gene age) information
phylostrata_Cg <- read_delim("Cgi.gene_ages.tsv", delim = "\t", escape_double = FALSE, trim_ws = TRUE)

# Remove missing entries
phylostrata_Cg <- na.omit(phylostrata_Cg)

# Keep only genes present in both average expression and phylostratum data
Cg_average <- subset(Cg_average, subset = Cg_average$GeneID %in% phylostrata_Cg$GeneID)

# Match expression matrix with phylostrata
Phylo_expr_Cg <- MatchMap(phylostrata_Cg[, c("Phylostratum", "GeneID")], Cg_average)

# Compute TAI (Transcriptome Age Index) per cell type
TAI <- as.data.frame(TAI(Phylo_expr_Cg))
TAI <- rownames_to_column(TAI, var = "cell_type")
colnames(TAI) <- c("cell_type", "TAI")

# Plot TAI per cell type
Oyster_TAI <- ggplot(TAI, aes(y = reorder(cell_type, TAI), x = TAI, color = cell_type)) +
  geom_point(size = 5) +
  scale_color_manual(values = my_colors) +
  theme_bw() +
  theme(
    axis.text.x = element_text(color = "black"),
    axis.text.y = element_text(color = "black"),
    axis.title = element_text(color = "black"),
    axis.ticks = element_line(color = "black"),
    plot.title = element_text(color = "black", hjust = 0.5),
    legend.position = "none"
  ) +
  labs(y = "Cell Type", x = "TAI") +
  scale_y_discrete(position = "right") +
  scale_x_continuous(expand = expansion(mult = c(0.05, 0.05)))
ggsave("Oyster_larva_TAI_celltype.pdf", plot = Oyster_TAI, width = 6, height = 8)

# Phylostratum enrichment of cell type marker genes
load(file = "Cg_cell_markers.rdata")
Cg_cell_markers <- rownames_to_column(as.data.frame(Cg_cell_markers), "GeneID")

# Subset marker genes present in phylostrata data
markers_phylo_cg <- subset(phylostrata_Cg, subset = phylostrata_Cg$GeneID %in% Cg_cell_markers$GeneID)

# Merge marker genes with phylostrata info
markers_phylo_cg <- left_join(markers_phylo_cg, Cg_cell_markers, by = "GeneID")

# Remove missing data and select relevant columns
oyster <- na.omit(select(markers_phylo_cg, cluster, PS_name, GeneID, Phylostratum))

# Hypergeometric test for enrichment of each phylostratum per cluster
Cg_phyp <- oyster %>%
  group_by(cluster) %>% mutate(n_cluster = n()) %>%
  group_by(PS_name) %>% mutate(n_strata = n_distinct(GeneID)) %>%
  ungroup() %>% mutate(n_genes = n_distinct(GeneID)) %>%
  group_by(cluster, PS_name) %>% mutate(n = n()) %>%
  select(-GeneID) %>% unique() %>% rowwise() %>%
  mutate(phyp = phyper(n, n_strata, n_genes - n_strata, n_cluster, lower.tail = FALSE))

# Transform enrichment results into matrix for heatmap
Cg_phyp_mat <- Cg_phyp %>%
  select(cluster, Phylostratum, phyp) %>% 
  pivot_wider(names_from = cluster, id_cols = Phylostratum, values_from = phyp) %>%
  left_join(unique(select(oyster, PS_name, Phylostratum)))

Cg_phyp_mat[is.na(Cg_phyp_mat)] <- 1  # replace missing values with 1

# Convert to data frame and sort by Phylostratum
Cg_phyp_mat <- as.data.frame(Cg_phyp_mat)
Cg_phyp_mat <- Cg_phyp_mat[order(as.numeric(Cg_phyp_mat$Phylostratum)), ]

# Set rownames
rownames(Cg_phyp_mat) <- Cg_phyp_mat$PS_name
Cg_phyp_mat$PS_name <- NULL
Cg_phyp_mat$Phylostratum <- NULL

# Plot pheatmap for phylostratum enrichment (-log10 scale)
colors <- colorRampPalette(c("#f7f4f9", "#cab2d6", "#6a3d9a"))(100)
Oyster_larva_pheatmap_celltype <- pheatmap(
  -log10(as.matrix(Cg_phyp_mat + 1e-10)),
  breaks = seq(0, 5, length.out = 100),
  color = colors,
  cluster_rows = FALSE,
  cluster_cols = TRUE,
  fontsize = 10,
  fontsize_row = 8,
  fontsize_col = 8,
  border_color = NA,
  show_rownames = TRUE,
  show_colnames = TRUE,
  legend = TRUE
)
ggsave("Oyster_larva_pheatmap_celltype.pdf", plot = Oyster_larva_pheatmap_celltype, width = 8, height = 10)

# Save matrix as CSV
write.csv(Cg_phyp_mat, file = "Cg_phyp_mat.csv", row.names = TRUE)
