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
  "Chaetal sacs" = "#C5362B",
  "Piwi+ cells" = "#3677A5",
  "Epidermis" = "#479EB5",
  "Muscle" = "#36A7A4",
  "Gut" = "#4B899B",
  "Neurons" = "#529192",
  "Eleocytes" = "#5D998A",
  "Globin+ cells" = "#72A87A",
  "Polycystin cells" = "#88B859",
  "Lipoxygenase+ cells" = "#A1C265",
  "Vigilin+ cells" = "#CCC45D",
  "Lumbrokinase+ cells" = "#E1C45E",
  "Carbohydrate metabolic cells" = "#F6C461",
  "Secretory cells" = "#EB9C61",
  "Ldlrr+ cells" = "#C49B54",
  "Unknown" = "#ADD3D0")

# Load cell type average expression data
load(file = "Ple_average.rdata")

# Load phylostratum (gene age) information
phylostrata_Ple <- read_delim("gene_ages.tsv", delim = "\t", escape_double = FALSE, trim_ws = TRUE)

# Remove missing entries
phylostrata_Ple <- na.omit(phylostrata_Ple)

# Keep only genes present in both average expression and phylostratum data
Ple_average <- subset(Ple_average, subset = Ple_average$GeneID %in% phylostrata_Ple$GeneID)

# Match expression matrix with phylostrata
Phylo_expr_Ple <- MatchMap(phylostrata_Ple[,c("Phylostratum","GeneID")], Ple_average)

# Compute TAI (Transcriptome Age Index) per cell type
TAI <- as.data.frame(TAI(Phylo_expr_Ple))
TAI <- rownames_to_column(TAI, var= "cell_type")
colnames(TAI) <- c("cell_type","TAI")

# Plot TAI per cell type
Ple_TAI <- ggplot(TAI, aes(y = reorder(cell_type, TAI), x = TAI, color = cell_type)) +
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
ggsave("Ple_TAI_celltype.pdf", plot = Ple_TAI, width = 6, height = 6)

# Phylostratum enrichment of cell type marker genes
load(file = "Ple_cell_markers.rdata")
Ple_cell_markers <- rownames_to_column(as.data.frame(Ple_cell_markers), "GeneID")

# Subset marker genes present in phylostrata data
markers_phylo_Ple <- subset(phylostrata_Ple, subset = phylostrata_Ple$GeneID %in% Ple_cell_markers$GeneID)

# Merge marker genes with phylostrata info
markers_phylo_Ple <- left_join(markers_phylo_Ple, Ple_cell_markers, by= c("GeneID"="GeneID"))

# Remove missing data and select relevant columns
Ple <- na.omit(select(markers_phylo_Ple, cluster, PS_name, GeneID, Phylostratum))

# Hypergeometric test for enrichment of each phylostratum per cluster
Ple_phyp <- Ple %>%
  group_by(cluster) %>% mutate(n_cluster=n()) %>%
  group_by(PS_name) %>% mutate(n_strata=n_distinct(GeneID)) %>%
  ungroup %>% mutate(n_genes=n_distinct(GeneID)) %>%
  group_by(cluster, PS_name) %>% mutate(n=n()) %>%
  select(-GeneID) %>% unique() %>% rowwise() %>%
  mutate(phyp=phyper(n, n_strata, n_genes - n_strata, n_cluster, lower.tail = F))  -> Ple_phyp

# Transform enrichment results into matrix for heatmap
Ple_phyp_mat <- Ple_phyp %>%
  select(cluster, Phylostratum, phyp) %>% 
  pivot_wider(names_from = cluster, id_cols = Phylostratum ,values_from = phyp) %>%
  left_join(unique(select(Ple, PS_name, Phylostratum))) 

Ple_phyp_mat[is.na(Ple_phyp_mat)] <- 1 # replace missing values with 1

# Convert to data frame and sort by Phylostratum
Ple_phyp_mat <- as.data.frame(Ple_phyp_mat)
Ple_phyp_mat <- Ple_phyp_mat[order(as.numeric(Ple_phyp_mat$Phylostratum)), ]

# Set rownames
rownames(Ple_phyp_mat) <- Ple_phyp_mat$Phylostratum
Ple_phyp_mat$Phylostratum <- NULL
rownames(Ple_phyp_mat) <- Ple_phyp_mat$PS_name
Ple_phyp_mat$PS_name <- NULL

# Plot pheatmap for phylostratum enrichment (-log10 scale)
colors <- colorRampPalette(c("#fff5f0", "#fcbba1", "#fc9272", "#fb6a4a"))(100)
Ple_pheatmap_celltype <- pheatmap(
  -log10(as.matrix(Ple_phyp_mat + 1e-10)),
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
ggsave("Ple_pheatmap_celltype.pdf", plot = Ple_pheatmap_celltype, width = 8, height = 10)

# Save matrix as CSV
write.csv(Ple_phyp_mat, file = "Ple_phyp_mat.csv", row.names = TRUE)
