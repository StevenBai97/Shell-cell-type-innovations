# =============================================================================
# PART I. Transcriptome age and phylostratum enrichment
# =============================================================================

# Load libraries
library(tidyverse)
library(ggplot2)
library(tidyr)
library(readr)
library(pheatmap)
library(myTAI)
library(tibble)

# Define custom colors for cell types
my_colors <- c(
  "Shell gland1" = "#c5362b",
  "Shell gland2" = "#E86A60",
  "Shell gland3" = "#f08a74",
  "Shell gland4" = "#f5c7b9",
  "Proliferative cells" = "#f6c461",
  "Neurons" = "#ccc45d",
  "Hemocytes" = "#36729d",
  "Cilia" = "#9dbb65",
  "Myocytes" = "#72a87a",
  "Unknown" = "#b0d8d5"
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
  group_by(cluster) %>%
  mutate(
    n_cluster = n_distinct(GeneID)
  ) %>%
  group_by(PS_name) %>%
  mutate(
    n_strata = n_distinct(GeneID)
  ) %>%
  ungroup() %>%
  mutate(
    n_genes = n_distinct(GeneID)
  ) %>%
  group_by(cluster, PS_name, Phylostratum) %>%
  mutate(
    n = n_distinct(GeneID)
  ) %>%
  select(
    cluster,
    PS_name,
    Phylostratum,
    n,
    n_cluster,
    n_strata,
    n_genes
  ) %>%
  distinct() %>%
  rowwise() %>%
  mutate(
    phyp = phyper(
      n - 1,
      n_strata,
      n_genes - n_strata,
      n_cluster,
      lower.tail = FALSE
    )
  ) %>%
  ungroup() %>%
  mutate(
    FDR = p.adjust(phyp, method = "BH")
  )

# Transform enrichment results into matrix for heatmap
Cg_phyp_mat <- Cg_phyp %>%
  select(
    cluster,
    PS_name,
    Phylostratum,
    phyp
  ) %>%
  pivot_wider(
    names_from = cluster,
    values_from = phyp,
    values_fill = 1
  ) %>%
  arrange(as.numeric(Phylostratum))

PS_info <- Cg_phyp_mat %>%
  select(PS_name, Phylostratum)

Cg_phyp_mat2 <- Cg_phyp_mat %>%
  select(-PS_name, -Phylostratum) %>%
  as.data.frame()

rownames(Cg_phyp_mat2) <- paste0(
  PS_info$PS_name
)

Cg_phyp_logP <- -log10(
  as.matrix(Cg_phyp_mat2) + 1e-300
)
Cg_phyp_logP[Cg_phyp_logP > 5] <- 5

# Plot pheatmap for phylostratum enrichment
Oyster_larva_pheatmap_celltype <- pheatmap(
  Cg_phyp_logP,
  color = colorRampPalette(c("white", "#6a3d9a"))(100),
  breaks = seq(0, 5, length.out = 101),
  cluster_rows = FALSE,
  cluster_cols = TRUE,
  fontsize = 10,
  fontsize_row = 8,
  fontsize_col = 8,
  border_color = NA,
  show_rownames = TRUE,
  show_colnames = TRUE,
  legend = TRUE,
  number_color = "black",
  fontsize_number = 10
)

ggsave("Oyster_larva_pheatmap_celltype.pdf", plot = Oyster_larva_pheatmap_celltype, width = 8, height = 10)

# Save matrix as CSV
write.csv(Cg_phyp_mat, file = "Cg_phyp_mat.csv", row.names = TRUE)

# ============================================================
# PART II. Duplication enrichment in larval cell-type markers
# ============================================================

# Read paralog pairs and ID mapping
para <- read.delim(
  "cgi_paralogs.txt",
  header = TRUE,
  stringsAsFactors = FALSE,
  check.names = FALSE
)

idmap <- read.delim(
  "Cgi.pep2gene",
  header = TRUE,
  stringsAsFactors = FALSE,
  check.names = FALSE
)

# Convert paralog pairs to gene-level duplication table
dup_gene <- bind_rows(
  para %>%
    transmute(
      XP_ID = gene1,
      duplication_node,
      Orthogroup
    ),
  para %>%
    transmute(
      XP_ID = gene2,
      duplication_node,
      Orthogroup
    )
) %>%
  distinct() %>%
  left_join(
    idmap,
    by = "XP_ID"
  ) %>%
  filter(
    !is.na(GeneID),
    GeneID != ""
  ) %>%
  select(
    GeneID,
    duplication_node,
    Orthogroup
  ) %>%
  distinct()

# Gene-node table used for enrichment
dup_gene_node <- dup_gene %>%
  distinct(
    GeneID,
    duplication_node
  )

# Load scRNA-seq marker and Seurat data
load("data/Cgigas/larva/Cg_cell_markers.rdata")
load("data/Cgigas/larva/OysterLarva.rdata")

# Prefer explicitly using RNA assay features
detected_genes <- rownames(seob)

orthofinder_genes <- idmap %>%
  filter(
    !is.na(GeneID),
    GeneID != ""
  ) %>%
  pull(GeneID) %>%
  unique()

universe <- intersect(
  detected_genes,
  orthofinder_genes
)

N <- length(universe)

# Prepare duplication gene sets
dup_gene_use <- dup_gene_node %>%
  filter(
    GeneID %in% universe
  )

dup_sets <- split(
  dup_gene_use$GeneID,
  dup_gene_use$duplication_node
) %>%
  lapply(unique)

# Number of duplicated genes per node after background filtering
dup_count <- tibble(
  duplication_node = names(dup_sets),
  N_genes = lengths(dup_sets)
) %>%
  arrange(desc(N_genes))

# Prepare marker sets
markers_use <- Cg_cell_markers %>%
  filter(
    gene %in% universe
  ) %>%
  distinct(
    cluster,
    gene
  )

marker_sets <- split(
  markers_use$gene,
  markers_use$cluster
) %>%
  lapply(unique)

# Hypergeometric enrichment
res <- expand_grid(
  cluster = names(marker_sets),
  duplication_node = names(dup_sets)
) %>%
  mutate(
    stats = map2(
      cluster,
      duplication_node,
      function(ct, node) {

        marker_genes <- marker_sets[[ct]]
        dup_genes <- dup_sets[[node]]

        n <- length(marker_genes)
        K <- length(dup_genes)
        k <- length(
          intersect(
            marker_genes,
            dup_genes
          )
        )

        expected <- n * K / N

        pvalue <- phyper(
          k - 1,
          K,
          N - K,
          n,
          lower.tail = FALSE
        )

        tibble(
          N_background = N,
          N_dup_genes = K,
          N_markers = n,
          N_overlap = k,
          Expected = expected,
          FoldEnrichment = ifelse(
            expected > 0,
            k / expected,
            NA_real_
          ),
          Pvalue = pvalue
        )
      }
    )
  ) %>%
  unnest(stats) %>%
  mutate(
    FDR = p.adjust(
      Pvalue,
      method = "BH"
    ),

    minus_log10_FDR = -log10(
      pmax(
        FDR,
        .Machine$double.xmin
      )
    )
  )

# Save complete enrichment results
write.table(res, file = "Cgi_larva_duplication_enrichment.tsv", sep = "\t", quote = FALSE, row.names = FALSE)

# Output duplicated marker genes with node and Orthogroup
marker_dup_detail <- Cg_cell_markers %>%
  select(
    gene,
    cluster,
    p_val,
    avg_log2FC,
    pct.1,
    pct.2,
    p_val_adj
  ) %>%
  inner_join(
    dup_gene,
    by = c("gene" = "GeneID"),
    relationship = "many-to-many"
  ) %>%
  filter(
    gene %in% universe
  ) %>%
  distinct(
    gene,
    cluster,
    p_val,
    avg_log2FC,
    pct.1,
    pct.2,
    p_val_adj,
    duplication_node,
    Orthogroup
  ) %>%
  arrange(
    cluster,
    duplication_node,
    Orthogroup,
    gene
  )

write.table(marker_dup_detail, file = "Cgi_larva_duplicated_marker_genes.tsv", sep = "\t", quote = FALSE, row.names = FALSE)

# Duplication node order
node_order <- c(
  "Metazoa",
  "Eumetazoa",
  "Bilateria",
  "Protostomia",
  "Lophotrochozoa",
  "Mollusca",
  "Bivalvia",
  "Pteriomorphia",
  "Pterioida",
  "Ostreidae",
  "Crassostrea",
  "C.gigas"
)

node_order <- intersect(
  node_order,
  unique(res$duplication_node)
)

# Build -log10(FDR) heatmap matrix
heatmap_mat <- res %>%
  mutate(
    duplication_node = factor(
      duplication_node,
      levels = node_order
    ),
    heat_value = ifelse(
      N_overlap == 0,
      0,
      minus_log10_FDR
    )
  ) %>%
  select(
    cluster,
    duplication_node,
    heat_value
  ) %>%
  pivot_wider(
    names_from = duplication_node,
    values_from = heat_value,
    values_fill = 0
  ) %>%
  column_to_rownames("cluster") %>%
  as.matrix()

# Force phylogenetic order
heatmap_mat <- heatmap_mat[
  ,
  node_order,
  drop = FALSE
]

# Plot heatmap
my_colors <- colorRampPalette(c("white", "#6a3d9a"))(100)

pdf("Cgi_larva_duplication_enrichment_heatmap.pdf", width = 8, height = 6)

pheatmap(
  heatmap_mat,
  cluster_rows = TRUE,
  cluster_cols = FALSE,
  color = my_colors,
  border_color = NA,
  angle_col = 45,
  fontsize = 10,
  fontsize_row = 10,
  fontsize_col = 9,
  legend = TRUE
)

dev.off()

# Bubble plot
cluster_order <- c(
  "Shell gland1",
  "Shell gland2",
  "Shell gland3",
  "Shell gland4",
  "Hemocytes",
  "Myocytes",
  "Cilia",
  "Neurons",
  "Proliferative cells",
  "Unknown"
)

plot_df <- res %>%
  mutate(
    duplication_node = factor(
      duplication_node,
      levels = node_order
    ),
    cluster = factor(
      cluster,
      levels = rev(cluster_order)
    ),
    log2FE = log2(FoldEnrichment),
    log2FE_plot = pmax(log2FE, 0),
    significance = case_when(
      FDR < 0.05 & FoldEnrichment > 1 ~ "Significant enrichment",
      TRUE ~ "Not significant"
    )
  ) %>%
  filter(
    N_overlap > 0
  )

p_bubble <- ggplot(
  plot_df,
  aes(
    x = duplication_node,
    y = cluster
  )
) +
  geom_point(
    aes(
      size = N_overlap,
      color = log2FE_plot
    )
  ) +
  scale_size_continuous(
    name = "Duplicated marker genes",
    range = c(2, 10)
  ) +
  scale_color_gradient(
    low = "white",
    high = "#6a3d9a",
    name = expression(log[2]("Fold enrichment"))
  ) +
  scale_y_discrete(
    position = "right"
  ) +
  labs(
    x = "Duplication age",
    y = NULL
  ) +
  theme_bw() +
  theme(
    axis.text.x = element_text(
      angle = 45,
      hjust = 1
    ),
    axis.text.y.left = element_blank(),
    axis.ticks.y.left = element_blank(),
    axis.text.y.right = element_text(
      size = 9
    ),
    axis.ticks.y.right = element_line(),
    legend.title = element_text(size = 9),
    legend.text = element_text(size = 8)
  )

pdf("Cgi_larva_duplication_enrichment_bubble.pdf", width = 8, height = 6)
p_bubble
dev.off()
