# Load required libraries
library(future)
library(tidyverse)
library(Seurat)
library(patchwork)
library(clustree)
library(ggplot2)
library(DoubletFinder)
library(shadowtext)

# Increase memory limit for parallelization
options(future.globals.maxSize = 50 * 1024^3)

# Initialize list to store Seurat objects per sample
seob_list <- list()

# Get sample names (directory names under data/Cgigas/mantle)
samples = basename(list.dirs("data/Cgigas/mantle", recursive = F))

# Loop through samples and create individual Seurat objects
for (sample in samples) {
  scrna_data <- Read10X(
    data.dir = str_c("data/Cgigas/mantle/", sample))
  seob <- CreateSeuratObject(
    counts = scrna_data,
    min.cells = 3,
    min.features = 300)
  seob[['samples']] <- sample
  seob_list[[sample]] = seob
}

# Merge all Seurat objects into one
seob <- merge(x = seob_list[[1]], y = seob_list[-1], add.cell.ids = names(seob_list))

# Calculate mitochondrial gene percentage (using oyster mt gene list)
seob[["percent.mt"]] <- PercentageFeatureSet(seob, features = c("CYTB", "COX2", "ATP6", "ND2", "ND4", "ND5", "ND6",  "ND3", "ND1", "ND4L", "COX1", "COX3"))

# QC violin plots before filtering
p0 <- VlnPlot(seob,
              features = c("nFeature_RNA", "nCount_RNA", "percent.mt"),
              assay = "RNA",
              layer = "counts",
              group.by  = "samples",
              log = F,
              pt.size = 0)

# Median ± 3 * MAD filtering
cells_filtered <- seob@meta.data %>%
  rownames_to_column(var = "cell_id") %>%
  group_by(samples) %>%
  mutate(
    nCount_median = median(nCount_RNA),
    nCount_mad = mad(nCount_RNA),
    nFeature_median = median(nFeature_RNA),
    nFeature_mad = mad(nFeature_RNA),
    percent_mt_median = median(percent.mt),
    percent_mt_mad = mad(percent.mt)
  ) %>%
  filter(
    nCount_RNA >= nCount_median - 3*nCount_mad,
    nCount_RNA <= nCount_median + 3*nCount_mad,
    nFeature_RNA >= nFeature_median - 3*nFeature_mad,
    nFeature_RNA <= nFeature_median + 3*nFeature_mad,
    percent.mt >= percent_mt_median - 3*percent_mt_mad,
    percent.mt <= percent_mt_median + 3*percent_mt_mad
  ) %>%
  pull(cell_id)

# Subset object with high-quality cells
seob_filter2 <- subset(seob, cells = cells_filtered)

# QC violin plots after filtering
p2 <- VlnPlot(seob_filter2,
              features = c("nFeature_RNA", "nCount_RNA", "percent.mt"),
              assay = "RNA",
              layer = "counts",
              group.by  = "samples",
              log = F,
              pt.size = 0)

# Save before/after filtering QC plots
pdf("data/Cgigas/mantle/Cgi.filter_plot.pdf", width = 10, height = 12)
p0 /p2
dev.off()

# Overwrite seob with filtered object
seob <- seob_filter2

# SCTransform normalization (no regression at this stage)
seob_sct_noregress <- SCTransform(seob)

## Read S phase and G2M phase gene IDs from file
s_genes.Cgi <- readLines("s_genes.Cgi.geneid")
g2m_genes.Cgi <- readLines("g2m_genes.Cgi.geneid")

# Assign cell cycle scores based on S and G2/M orthologs
seob_sct_noregress <- CellCycleScoring(
  seob_sct_noregress,
  s.features = s_genes.Cgi,
  g2m.features = g2m_genes.Cgi
)

# Re-run SCTransform, regressing out cell cycle effects (S and G2M scores)
seob_sct_regress <- SCTransform(seob_sct_noregress, vars.to.regress = c("S.Score", "G2M.Score"))

# Update object
seob <- seob_sct_regress

# Dimensionality reduction (PCA, UMAP)
seob <- RunPCA(seob, assay = "SCT")
seob <- RunUMAP(seob, assay = "SCT", dims = 1:30)

# PCA and UMAP plots grouped by sample
p_pca <- DimPlot(seob,
                 reduction = "pca",
                 group.by = "samples",
                 label = F)
p_umap <- DimPlot(seob,
                  reduction = "umap",
                  group.by = "samples",
                  label = F)
pdf("data/Cgigas/mantle/Cgi_PCA_UMAP.pdf", width = 12, height = 6)
print(p_pca + p_umap)
dev.off()

# Graph-based clustering
seob <- FindNeighbors(seob, reduction = "pca", dims = 1:30)

# Explore clustering resolutions (0.5 to 2.0)
resolutions <- seq(0.5, 2, by=0.1)
seob <- FindClusters(seob,
                     resolution = resolutions,
                     random.seed = 1)

# Clustree visualization to guide resolution choice
clustree_plot <- clustree(seob, prefix = "SCT_snn_res.")
ggsave("data/Cgigas/mantle/clustree_plot.pdf", plot = clustree_plot, width = 40, height = 20)

# Re-run clustering with optimal resolution (1.2 chosen here)
seob <- FindNeighbors(seob, reduction = "pca", dims = 1:30)
seob <- FindClusters(seob, resolution = 1.2, random.seed = 1)

# Plot UMAP colored by final clusters
cluster_plot <- DimPlot(seob,
                        reduction = "umap",
                        group.by  = c("seurat_clusters"),
                        label = T)
ggsave("data/Cgigas/mantle/cluster_plot.pdf", plot = cluster_plot, width = 11, height = 9)

# Doublet removal using DoubletFinder
## 1. Split by sample
seob_list <- SplitObject(seob, split.by = "samples")

for (i in 1:length(seob_list)) {
  s <- seob_list[[i]]
  sample_name <- names(seob_list)[i]
  
  ## 2. Estimate expected number of doublets
  nExp_poi <- round(0.08 * nrow(s@meta.data))
  homotypic.prop <- modelHomotypic(s$seurat_clusters)
  nExp_poi.adj <- round(nExp_poi * (1 - homotypic.prop))
  
  ## 3. Choose optimal pK
  sweep.res.list <- paramSweep(s, PCs = 1:30, sct = TRUE)
  sweep.stats <- summarizeSweep(sweep.res.list, GT = FALSE)
  bcmvn <- find.pK(sweep.stats)
  bcmvn$pK <- as.numeric(as.character(bcmvn$pK))
  pK_candidates <- bcmvn$pK[bcmvn$pK >= 0.01]
  BC_candidates <- bcmvn$BCmetric[bcmvn$pK >= 0.01]
  if (length(pK_candidates) == 0) {
    warning("No valid pK values >= 0.01. Skipping this sample.")
    next
  }
  pK <- pK_candidates[which.max(BC_candidates)]
  message("Sample: ", sample_name, " | Selected pK: ", pK)
  
  ## 4. Run DoubletFinder
  s <- doubletFinder(s, PCs = 1:30, pN = 0.25, pK = pK,
                     nExp = nExp_poi.adj, sct = TRUE)
  seob_list[[i]] <- s  # Update object
}

# Print classification summary for each sample
for (i in seq_along(seob_list)) {
  s <- seob_list[[i]]
  sample_name <- names(seob_list)[i]
  
  df_col <- grep("^DF.classifications_", colnames(s@meta.data), value = TRUE)
  if (length(df_col) == 0) {
    cat("Sample", sample_name, ": No DoubletFinder classification found.\n\n")
    next
  }
  
  cat("Sample", sample_name, "classification summary:\n")
  print(table(s@meta.data[[df_col[length(df_col)]]]))
  cat("\n")
}

# Extract all Singlet cell names across samples
Singlets <- unlist(lapply(seob_list, function(s) {
  df_col <- grep("^DF.classifications_", colnames(s@meta.data), value = TRUE)
  if (length(df_col) == 0) return(NULL)
  rownames(s@meta.data)[s@meta.data[[df_col[length(df_col)]]] == "Singlet"]
}))

# Subset the original Seurat object to retain only Singlets
seob_singlets <- subset(seob, cells = Singlets)

# Overwrite seob with singlet-only object for downstream analysis
seob <- seob_singlets

# Read eggNOG-mapper annotation
emapper <- read.delim("data/Cgigas/Cgi.pep.emapper_out.annotations", comment.char = "#", header = TRUE, sep = "\t")

# Replace missing "Preferred_name" with query ID
emapper$Preferred_name[emapper$Preferred_name == "-"] <- NA
emapper$Preferred_name[is.na(emapper$Preferred_name)] <- emapper$query[is.na(emapper$Preferred_name)]

# Create mapping from geneID to preferred gene name
gene2name <- unique(emapper[, c("geneID", "Preferred_name")])
gene_label_map <- setNames(gene2name$Preferred_name, gene2name$geneID)

# Define the list of marker genes
ordered_genes <- c(
  "LOC105329954", "LOC105335295", "LOC105325156", "LOC105324712", "LOC105344040", "LOC117688780",
  "LOC105336368", "LOC105336346", "LOC105347390", "Tyr", "LOC105334492", "LOC105326482", 
  "LOC105324352", "LOC109617045", "LOC105347705", "LOC105335177", "LOC105334556",
  "LOC105329634", "LOC105348570", "LOC105348490", "LOC105329481", "LOC105335713",
  "LOC105323195", "LOC105327706", "LOC105343071", "LOC105320177", "LOC105331920",
  "LOC105338883", "LOC105338950", "LOC105324772"， "LOC105326938", "LOC105328755",
  "LOC105329247", "LOC105319677", "LOC105323511", "LOC105328903", "LOC105333128",
  "LOC105337940", "LOC105339039", "LOC105341104"
)

# Extract expression data for DotPlot
marker_data <- DotPlot(seob, 
                   assay = "SCT", 
                   features = ordered_genes,
                   group.by = "cell_type")$data
marker_data <- marker_data[marker_data$pct.exp > 0, ]

# Replace gene IDs with preferred names for plotting
marker_data$features.plot <- ifelse(
  is.na(gene_label_map[as.character(marker_data$features.plot)]),
  as.character(marker_data$features.plot),
  gene_label_map[as.character(marker_data$features.plot)]
)

# Set factor levels to maintain the order of genes in the plot
feature_levels <- sapply(ordered_genes, function(gid) {
  if (!is.na(gene_label_map[gid])) {
    gene_label_map[gid]
  } else {
    gid
  }
})
marker_data$features.plot <- factor(marker_data$features.plot, levels = feature_levels)

# Define cluster order for plotting
cluster_order <- rev(c(22,18,27,29,33,21,1,10,32,8,17,24,9,25,31,4,7,11,34,5,30,12,14,23,0,2,3,6,13,15,16,19,20,26,28))
marker_data$id <- factor(marker_data$id, levels = as.character(cluster_order))

# Generate DotPlot using ggplot2
marker_plot <- ggplot(marker_data, aes(x = features.plot, y = id)) +
  geom_point(aes(size = pct.exp, color = avg.exp.scaled)) +
  scale_color_gradient2(low = "#4575b4", mid = "white", high = "#d73027", midpoint = 0) +
  scale_size(range = c(1, 8)) +
  theme_bw() +
  theme(
    axis.text.x = element_text(angle = 45, hjust = 1, vjust = 1, 
                               size = 10, color = "black"),
    axis.text.y = element_text(size = 10, color = "black"),
    axis.title = element_text(size = 10, color = "black"),
    panel.grid.major = element_blank(),
    panel.grid.minor = element_blank(),
    legend.title = element_text(size = 10, color = "black"),
    legend.text = element_text(size = 10, color = "black")
  )
pdf("data/Cgigas/mantle/marker_dotplot.pdf", width = 12, height = 15)
print(marker_plot)
dev.off()

# Define mapping from cluster ID to cell type
cluster2type <- c(
  "0" = "Unknown",
  "1" = "Hemocytes (small granule cells)",
  "2" = "Unknown",
  "3" = "Unknown",
  "4" = "Myocytes",
  "5" = "Cilia",
  "6" = "Unknown",
  "7" = "Myocytes",
  "8" = "Hemocytes (blast-like cells)",
  "9" = "Hemocytes (hyalinocytes)",
  "10" = "Hemocytes (small granule cells)",
  "11" = "Myocytes",
  "12" = "Neurons",
  "13" = "Unknown",
  "14" = "Vesicular connective tissue cells",
  "15" = "Unknown",
  "16" = "Unknown",
  "17" = "Hemocytes (blast-like cells)",
  "18" = "Shell-secreting epithelial cells2",
  "19" = "Unknown",
  "20" = "Unknown",
  "21" = "Shell-secreting epithelial cells5",
  "22" = "Shell-secreting epithelial cells1",
  "23" = "Proliferative cells",
  "24" = "Hemocytes (blast-like cells)",
  "25" = "Hemocytes (vesicular cells)",
  "26" = "Unknown",
  "27" = "Shell-secreting epithelial cells3",
  "28" = "Unknown",
  "29" = "Shell-secreting epithelial cells3",
  "30" = "Cilia",
  "31" = "Myocytes + Hemocytes",
  "32" = "Hemocytes (small granule cells)",
  "33" = "Shell-secreting epithelial cells4",
  "34" = "Cilia + Hemocytes"
)

# Add cell type annotation to Seurat metadata
seob@meta.data$cell_type <- cluster2type[as.character(seob@meta.data$seurat_clusters)]

# Create cell type label for plotting (omit cluster numbers)
seob@meta.data$cell_type_label <- seob@meta.data$cell_type

# Compute UMAP centroids for each cluster (for label placement)
umap_df <- as.data.frame(Embeddings(seob, "umap"))
colnames(umap_df) <- c("UMAP_1", "UMAP_2")
umap_df$cluster <- seob@meta.data$seurat_clusters

centroids <- umap_df %>%
  group_by(cluster) %>%
  summarise(
    UMAP_1 = median(UMAP_1),
    UMAP_2 = median(UMAP_2)
  )

# Define colors for each cell type (no cluster numbers in legend)
my_colors <- c(
  "Shell-secreting epithelial cells1" = "#B02E26",
  "Shell-secreting epithelial cells2" = "#C5362B",
  "Shell-secreting epithelial cells3" = "#E86A60",
  "Shell-secreting epithelial cells4" = "#f08a74",
  "Shell-secreting epithelial cells5" = "#F2A19A",
  "Hemocytes (small granule cells)" = "#357daf",
  "Hemocytes (blast-like cells)" = "#4185a8",
  "Hemocytes (hyalinocytes)" = "#4c8ea1",
  "Hemocytes (vesicular cells)" = "#549697",
  "Myocytes + Hemocytes" = "#5f9e8e",
  "Myocytes" = "#75ad7d",
  "Cilia + Hemocytes" = "#8dc05a",
  "Cilia" = "#a6c966",
  "Neurons" = "#d3cb5d",
  "Vesicular connective tissue cells" = "#e9ca5d",
  "Proliferative cells" = "#FFC961",
  "Unknown" = "#b0d8d5"
)

# Set factor levels to maintain consistent legend order
seob@meta.data$cell_type_label <- factor(seob@meta.data$cell_type_label, levels = names(my_colors))

# Plot UMAP with custom colors and cluster labels at centroids
cluster_plot <- DimPlot(
  seob,
  reduction = "umap",
  group.by = "cell_type_label",
  label = FALSE,
  cols = my_colors
) +
  geom_shadowtext(
    data = centroids,
    aes(x = UMAP_1, y = UMAP_2, label = cluster),
    color = "black",
    bg.color = "white",
    bg.r = 0.15,
    size = 5,
    fontface = "bold",
    inherit.aes = FALSE
  ) +
  theme_void(base_size = 12) +
  theme(legend.position = "right") +
  ggtitle(NULL)

# Add arrows to indicate UMAP axes
xrange <- ggplot_build(cluster_plot)$layout$panel_params[[1]]$x.range
yrange <- ggplot_build(cluster_plot)$layout$panel_params[[1]]$y.range
arrow_base_x <- xrange[1] + 0.1 * diff(xrange)
arrow_base_y <- yrange[1] + 0.1 * diff(yrange)

# Generate Plot
cluster_plot <- cluster_plot +
  annotate(
    "segment",
    x = arrow_base_x, xend = arrow_base_x + 0.1 * diff(xrange),
    y = arrow_base_y, yend = arrow_base_y,
    arrow = arrow(length = unit(0.2, "cm")), size = 0.5
  ) +
  annotate(
    "text",
    x = arrow_base_x + 0.1 * diff(xrange), y = arrow_base_y,
    label = "UMAP1", hjust = 0, vjust = 1, size = 3.5
  ) +
  annotate(
    "segment",
    x = arrow_base_x, xend = arrow_base_x,
    y = arrow_base_y, yend = arrow_base_y + 0.1 * diff(yrange),
    arrow = arrow(length = unit(0.2, "cm")), size = 0.5
  ) +
  annotate(
    "text",
    x = arrow_base_x + 0.05 * diff(xrange), y = arrow_base_y + 0.1 * diff(yrange),
    label = "UMAP2", hjust = 1, vjust = 0, size = 3.5
  ) +
  coord_fixed(ratio = 1)
pdf("data/Cgigas/mantle/cluster_plot.pdf", width = 15, height = 15)
print(cluster_plot)
dev.off()

# Identify marker genes for each cell type
Cg_cell_markers <- FindAllMarkers(seob, only.pos = TRUE, group.by = "cell_type", min.pct = 0.10, logfc.threshold = 0.25, assay = "SCT")
write.csv(Cg_cell_markers, file = "data/Cgigas/mantle/Cg_cell_markers.csv", row.names = FALSE)
save(Cg_cell_markers, file = "data/Cgigas/mantle/Cg_cell_markers.rdata")

# Compute average expression per cell type
Cg_average <- AverageExpression(seob, group.by = "cell_type")$SCT
Cg_average <- rownames_to_column(as.data.frame(Cg_average), "GeneID")
save(Cg_average, file = "data/Cgigas/mantle/Cg_average.rdata")

# Save the Seurat object
save(seob, file = "data/Cgigas/mantle/OysterMantle.scRNAseq.rdata")
