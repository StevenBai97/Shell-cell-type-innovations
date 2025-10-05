# Load required libraries
library(future)
library(tidyverse)
library(Seurat)
library(patchwork)
library(ggplot2)
library(DoubletFinder)

# Increase memory limit for parallelization
options(future.globals.maxSize = 50 * 1024^3)

# Subset Seurat object to retain only adult samples
seob <- readRDS("data/Pgotoi/Ch_v5_chim.RDS")
adult_cells <- rownames(seob@meta.data)[grepl("^Adult", seob$orig.ident)]
seob <- subset(seob, cells = adult_cells)

# QC violin plots before filtering
p0 <- VlnPlot(seob,
              features = c("nFeature_RNA", "nCount_RNA", "percent_mito"),
              assay = "RNA",
              layer = "counts",
              group.by  = "orig.ident",
              log = F,
              pt.size = 0)

# Median ± 3 * MAD filtering
cells_filtered <- seob@meta.data %>%
  rownames_to_column(var = "cell_id") %>%
  group_by(orig.ident) %>%
  mutate(
    nCount_median = median(nCount_RNA),
    nCount_mad = mad(nCount_RNA),
    nFeature_median = median(nFeature_RNA),
    nFeature_mad = mad(nFeature_RNA),
    percent_mt_median = median(percent_mito),
    percent_mt_mad = mad(percent_mito)
  ) %>%
  filter(
    nCount_RNA >= nCount_median - 3*nCount_mad,
    nCount_RNA <= nCount_median + 3*nCount_mad,
    nFeature_RNA >= nFeature_median - 3*nFeature_mad,
    nFeature_RNA <= nFeature_median + 3*nFeature_mad,
    percent_mito >= percent_mt_median - 3*percent_mt_mad,
    percent_mito <= percent_mt_median + 3*percent_mt_mad
  ) %>%
  pull(cell_id)

# Subset object with high-quality cells
seob_filter2 <- subset(seob, cells = cells_filtered)

# QC violin plots after filtering
p2 <- VlnPlot(seob_filter2,
              features = c("nFeature_RNA", "nCount_RNA", "percent_mito"),
              assay = "RNA",
              layer = "counts",
              group.by  = "orig.ident",
              log = F,
              pt.size = 0)

# Save before/after filtering QC plots
pdf("data/Pgotoi/Pgo.filter_plot.pdf", width = 10, height = 12)
p0 /p2
dev.off()

# Overwrite seob with filtered object
seob <- seob_filter2

# SCTransform normalization (no regression at this stage)
seob_sct_noregress <- SCTransform(seob)

# Read S phase and G2M phase gene IDs from file
s_genes.Pgo <- readLines("s_genes.Ple")
g2m_genes.Pgo <- readLines("g2m_genes.Ple")

# Assign cell cycle scores based on S and G2/M orthologs
seob_sct_noregress <- CellCycleScoring(
  seob_sct_noregress,
  s.features = s_genes.Pgo,
  g2m.features = g2m_genes.Pgo
)

# Re-run SCTransform, regressing out cell cycle effects (S and G2M scores)
seob_sct_regress <- SCTransform(seob_sct_noregress, vars.to.regress = c("S.Score", "G2M.Score"))

# Update object
seob <- seob_sct_regress

# Dimensionality reduction (PCA, UMAP)
seob <- RunPCA(seob, assay = "SCT")
seob <- RunUMAP(seob, assay = "SCT", dims = 1:50)

# PCA and UMAP plots grouped by sample
p_pca <- DimPlot(seob,
                 reduction = "pca",
                 group.by = "orig.ident",
                 label = F)
p_umap <- DimPlot(seob,
                  reduction = "umap",
                  group.by = "orig.ident",
                  label = F)
pdf("data/Pgotoi/Pgo_PCA_UMAP.pdf", width = 12, height = 6)
print(p_pca + p_umap)
dev.off()

# Graph-based clustering
seob <- FindNeighbors(seob, reduction = "pca", dims = 1:50)
seob <- FindClusters(seob, resolution = 1, random.seed = 1)

# Doublet removal using DoubletFinder
## 1. Split by sample
seob_list <- SplitObject(seob, split.by = "orig.ident")

for (i in 1:length(seob_list)) {
  s <- seob_list[[i]]
  samPgo_name <- names(seob_list)[i]
  
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
    warning("No valid pK values >= 0.01. Skipping this samPgo.")
    next
  }
  pK <- pK_candidates[which.max(BC_candidates)]
  message("SamPgo: ", samPgo_name, " | Selected pK: ", pK)
  
  ## 4. Run DoubletFinder
  s <- doubletFinder(s, PCs = 1:30, pN = 0.25, pK = pK,
                     nExp = nExp_poi.adj, sct = TRUE)
  seob_list[[i]] <- s  # Update object
}

# Print classification summary for each sample
for (i in seq_along(seob_list)) {
  s <- seob_list[[i]]
  samPgo_name <- names(seob_list)[i]
  
  df_col <- grep("^DF.classifications_", colnames(s@meta.data), value = TRUE)
  if (length(df_col) == 0) {
    cat("SamPgo", samPgo_name, ": No DoubletFinder classification found.\n\n")
    next
  }
  
  cat("SamPgo", samPgo_name, "classification summary:\n")
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

# Define the list of marker genes
ordered_genes <- c(
  "Chs2-5", "Chi10-4", "Dyh3-3", "Tbb6-1", "Tbb4-2", "Vwf-2", "Chs2-3", "Chit1-7",
  "Chit1-8", "Ifea-2", "Syt16-2", "Act", "Aqp10-2", "Fmar-5","Agrin-7",
  "Muc5a-2", "Aplp", "Tmps9-1", "Pros-1", "Titin-2", "Titin-3", "Titin-4", "Tnnt",
  "Tnni-3", "Tnni-2", "Act1-5", "Acts-1", "Cac1a-2", "Cac1a-3", "Tbx2", "Pxdn-9",
  "Syt15-3", "Piwl1-2", "Piwl1-5", "Ddx5-1", "Nano2-2", "Mael-2", "M10l1")

# Extract expression data for DotPlot
marker_data <- DotPlot(seob, 
                       assay = "SCT", 
                       features = ordered_genes,
                       group.by = "seurat_clusters")$data
marker_data <- marker_data[marker_data$pct.exp > 0, ]

# Define cluster order for plotting
cluster_order <- rev(c(22,10,14,20,25,5,12,23,7,13,3,11,27,30,6,17,21,18,19,16,31,2,4,8,15,28,33,0,1, 9, 24, 26, 29, 32))
marker_data$id <- factor(marker_data$id, levels = as.character(cluster_order))

# Generate DotPlot using ggplot2
marker_plot <- ggplot(marker_data, aes(x = features.plot, y = id)) +
  geom_point(aes(size = pct.exp, color = avg.exp.scaled)) +
  scale_color_gradient2(low = "#4575b4", mid = "white", high = "#d73027", midpoint = 0) +
  scale_size(range = c(1, 8)) +
  theme_bw() +
  theme(
    axis.text.x = element_text(angle = 90, hjust = 1, vjust = 1, 
                               size = 10, color = "black"),
    axis.text.y = element_text(size = 10, color = "black"),
    axis.title = element_text(size = 10, color = "black"),
    panel.grid.major = element_blank(),
    panel.grid.minor = element_blank(),
    legend.title = element_text(size = 10, color = "black"),
    legend.text = element_text(size = 10, color = "black")
  )
pdf("data/Pgotoi/marker_dotplot.pdf", width = 12, height = 15)
print(marker_plot)
dev.off()

# Define mapping from cluster ID to cell type
cluster2type <- c(
  "22" = "Chaete",
  "10" = "Cilia",
  "14" = "Cilia",
  "20" = "Cilia",
  "25" = "Cilia",
  "5" = "Epidermis",
  "12" = "Epidermis",
  "23" = "Epidermis",
  "7" = "Epidermis",
  "13" = "Fmar",
  "3" = "Gut",
  "11" = "Gut",
  "27" = "Gut",
  "30" = "Gut",
  "6" = "Hemocytes",
  "17" = "Muscle",
  "21" = "Muscle",
  "18" = "Neurons",
  "19" = "Neurons",
  "16" = "Papillae",
  "31" = "Papillae",
  "2" = "Stem cells",
  "4" = "Stem cells",
  "8" = "Stem cells",
  "15" = "Stem cells",
  "28" = "Stem cells",
  "33" = "Stem cells",
  "0" = "Unknown",
  "1" = "Unknown",
  "9" = "Unknown",
  "24" = "Unknown",
  "26" = "Unknown",
  "29" = "Unknown",
  "32" = "Unknown")

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
  summarise(UMAP_1 = median(UMAP_1),
            UMAP_2 = median(UMAP_2))

# Define colors for each cell type (no cluster numbers in legend)
my_colors <- c(
  "Chaete" = "#C5362B",
  "Cilia" = "#a1c265",
  "Epidermis" = "#479eb5",
  "Fmar" = "#72a87a",
  "Hemocytes" = "#3677a5",
  "Muscle" = "#36a7a4",
  "Neurons" = "#5d998a",
  "Papillae" = "#529192",
  "Gut" = "#4b899b",
  "Stem cells" = "#e1c45e",
  "Unknown" = "#add3d0")

# Set factor levels to maintain consistent legend order
seob@meta.data$cell_type_label <- factor(seob@meta.data$cell_type_label, levels = names(my_colors))

# Plot UMAP with custom colors and cluster labels at centroids
cluster_plot <- DimPlot(seob,
                        reduction = "umap",
                        group.by = "cell_type_label",
                        label = FALSE,
                        cols = my_colors) +
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
  annotate("segment",
           x = arrow_base_x, xend = arrow_base_x + 0.1 * diff(xrange),
           y = arrow_base_y, yend = arrow_base_y,
           arrow = arrow(length = unit(0.2, "cm")), size = 0.5) +
  annotate("text",
           x = arrow_base_x + 0.1 * diff(xrange), y = arrow_base_y,
           label = "UMAP1", hjust = 0, vjust = 1, size = 3.5) +
  annotate("segment",
           x = arrow_base_x, xend = arrow_base_x,
           y = arrow_base_y, yend = arrow_base_y + 0.1 * diff(yrange),
           arrow = arrow(length = unit(0.2, "cm")), size = 0.5) +
  annotate("text",
           x = arrow_base_x + 0.05 * diff(xrange), y = arrow_base_y + 0.1 * diff(yrange),
           label = "UMAP2", hjust = 1, vjust = 0, size = 3.5) +
  coord_fixed(ratio = 1)
pdf("data/Pgotoi/cluster_plot.pdf", width = 15, height = 15)
print(cluster_plot)
dev.off()

# Identify marker genes for each cell type
Pgo_cell_markers <- FindAllMarkers(seob, only.pos = TRUE, group.by = "cell_type", min.pct = 0.10, logfc.threshold = 0.25, assay = "SCT")
write.csv(Pgo_cell_markers, file = "data/Pgotoi/Pgo_cell_markers.csv", row.names = FALSE)
save(Pgo_cell_markers, file = "data/Pgotoi/Pgo_cell_markers.rdata")

# Compute average expression per cell type
Pgo_average <- AverageExpression(seob, group.by = "cell_type")$SCT
Pgo_average <- rownames_to_column(as.data.frame(Pgo_average), "GeneID")
save(Pgo_average, file = "data/Pgotoi/Pgo_average.rdata")

# Save the Seurat object
save(seob, file = "data/Pgotoi/Pgotoi.scRNAseq.rdata")
