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

# Get sample names (directory names under data/Pleidyi)
samples = basename(list.dirs("data/Pleidyi", recursive = F))

# Loop through samples and create individual Seurat objects
for (sample in samples) {
  scrna_data <- Read10X(
    data.dir = str_c("data/Pleidyi/", sample))
  seob <- CreateSeuratObject(
    counts = scrna_data,
    min.cells = 3,
    min.features = 50)
  seob[['samples']] <- sample
  seob_list[[sample]] = seob
}

# Merge all Seurat objects into one
seob <- merge(x = seob_list[[1]], y = seob_list[-1], add.cell.ids = names(seob_list))

# QC violin plots before filtering
p0 <- VlnPlot(seob, 
              features = c("nFeature_RNA", "nCount_RNA"), 
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
    nFeature_mad = mad(nFeature_RNA)) %>%
  filter(
    nCount_RNA >= nCount_median - 3*nCount_mad,
    nCount_RNA <= nCount_median + 3*nCount_mad,
    nFeature_RNA >= nFeature_median - 3*nFeature_mad,
    nFeature_RNA <= nFeature_median + 3*nFeature_mad) %>%
  pull(cell_id)

# Subset object with high-quality cells
seob_filter2 <- subset(seob, cells = cells_filtered)

# QC violin plots after filtering
p2 <- VlnPlot(seob_filter2, 
              features = c("nFeature_RNA", "nCount_RNA"), 
              assay = "RNA",
              layer = "counts",
              group.by  = "samples",
              log = F,
              pt.size = 0)

# Save before/after filtering QC plots
pdf("data/Pleidyi/Ple.filter_plot.pdf", width = 10, height = 12)
p0 /p2
dev.off()

# Overwrite seob with filtered object
seob <- seob_filter2

# SCTransform normalization (no regression at this stage)
seob_sct_noregress <- SCTransform(seob)

# Read S phase and G2M phase gene IDs from file
s_genes.Ple <- readLines("s_genes.Ple")
g2m_genes.Ple <- readLines("g2m_genes.Ple")

# Assign cell cycle scores based on S and G2/M orthologs
seob_sct_noregress <- CellCycleScoring(
  seob_sct_noregress,
  s.features = s_genes.Ple,
  g2m.features = g2m_genes.Ple
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
pdf("data/Pleidyi/Ple_PCA_UMAP.pdf", width = 12, height = 6)
print(p_pca + p_umap)
dev.off()

# Graph-based clustering
seob <- FindNeighbors(seob, reduction = "pca", dims = 1:30)

# Explore clustering resolutions (0.5 to 2.5)
resolutions <- seq(0.5, 2.5, by=0.1)
seob <- FindClusters(seob,
                     resolution = resolutions,
                     random.seed = 1)

# Clustree visualization to guide resolution choice
clustree_plot <- clustree(seob, prefix = "SCT_snn_res.")
ggsave("data/Pleidyi/clustree_plot.pdf", plot = clustree_plot, width = 40, height = 20)

# Re-run clustering with optimal resolution (1.2 chosen here)
seob <- FindNeighbors(seob, reduction = "pca", dims = 1:30)
seob <- FindClusters(seob, resolution = 1.5, random.seed = 1)

# Plot UMAP colored by final clusters
cluster_plot <- DimPlot(seob,
                        reduction = "umap",
                        group.by  = c("seurat_clusters"),
                        label = T)
ggsave("data/Pleidyi/cluster_plot.pdf", plot = cluster_plot, width = 11, height = 9)

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

# Define the list of marker genes
ordered_genes <- c(
  "PrileiEVm000573t1", "PrileiEVm000939t1", "PrileiEVm007502t1", "PrileiEVm003567t1", "PrileiEVm016887t1",
  "PrileiEVm016982t1", "PrileiEVm003521t1", "PrileiEVm022498t1", "PrileiEVm023936t1", "PrileiEVm008309t1",
  "PrileiEVm008287t1", "PrileiEVm011741t1", "PrileiEVm000300t1", "PrileiEVm019310t1", "PrileiEVm025675t1",
  "PrileiEVm020932t1", "PrileiEVm016172t1", "PrileiEVm015909t1", "PrileiEVm014978t1", "PrileiEVm018738t1",
  "PrileiEVm000447t1", "PrileiEVm005677t1", "PrileiEVm017310t1", "PrileiEVm026965t1", "PrileiEVm019805t1",
  "PrileiEVm010941t1", "PrileiEVm010132t1", "PrileiEVm000199t1", "PrileiEVm021761t1", "PrileiEVm001383t1",
  "PrileiEVm020317t1", "PrileiEVm022227t1", "PrileiEVm000558t1", "PrileiEVm014226t1", "PrileiEVm012030t1",
  "PrileiEVm013699t1", "PrileiEVm005681t1", "PrileiEVm009033t1", "PrileiEVm007557t1", "PrileiEVm018088t1",
  "PrileiEVm015446t1", "PrileiEVm020672t1", "PrileiEVm018061t1", "PrileiEVm020599t1", "PrileiEVm020147t1",
  "PrileiEVm005033t1", "PrileiEVm004079t1", "PrileiEVm000278t1", "PrileiEVm008087t1", "PrileiEVm008285t1",
  "PrileiEVm012391t1", "PrileiEVm016387t1", "PrileiEVm016330t1", "PrileiEVm016446t1", "PrileiEVm016633t1",
  "PrileiEVm001525t1", "PrileiEVm013994t1", "PrileiEVm026130t1", "PrileiEVm014911t1", "PrileiEVm008904t1",
  "PrileiEVm026983t1", "PrileiEVm002621t1", "PrileiEVm010163t1", "PrileiEVm014251t1", "PrileiEVm010669t1"
)

# Extract expression data for DotPlot
marker_data <- DotPlot(seob, 
                   assay = "SCT", 
                   features = ordered_genes,
                   group.by = "seurat_clusters")$data
marker_data <- marker_data[marker_data$pct.exp > 0, ]

# Define cluster order for plotting
cluster_order <- rev(c(36,49,26,2,4,14,18,6,0,5,9,22,23,27,46,47,38, 3,7,10,8,15,17,41,34,29,44,1,
                       12,24,30,35,13,19,20,16,33,11,21,31,43,25,39,40,45,37,28,32,42,48))
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
pdf("data/Pleidyi/marker_dotplot.pdf", width = 12, height = 15)
print(marker_plot)
dev.off()

# Define mapping from cluster ID to cell type
cluster2type <- c(
  "0" = "Epidermis",
  "1" = "Neurons",
  "2" = "Piwi+ cells",
  "3" = "Muscle",
  "4" = "Piwi+ cells",
  "5" = "Epidermis",
  "6" = "Piwi+ cells",
  "7" = "Muscle",
  "8" = "Gut",
  "9" = "Epidermis",
  "10" = "Muscle",
  "11" = "Polycystin cells",
  "12" = "Eleocytes",
  "13" = "Globin+ cells",
  "14" = "Piwi+ cells",
  "15" = "Gut",
  "16" = "Polycystin cells",
  "17" = "Gut",
  "18" = "Piwi+ cells",
  "19" = "Globin+ cells",
  "20" = "Globin+ cells",
  "21" = "Lipoxygenase+ cells",
  "22" = "Epidermis",
  "23" = "Epidermis",
  "24" = "Eleocytes",
  "25" = "Carbohydrate metabolic cells",
  "26" = "Chaetal sacs",
  "27" = "Epidermis",
  "28" = "Unknown",
  "29" = "Neurons",
  "30" = "Eleocytes",
  "31" = "Vigilin+ cells",
  "32" = "Unknown",
  "33" = "Polycystin cells",
  "34" = "Neurons",
  "35" = "Eleocytes",
  "36" = "Chaetal sacs",
  "37" = "Ldlrr+ cells",
  "38" = "Epidermis",
  "39" = "Secretory cells",
  "40" = "Secretory cells",
  "41" = "Gut",
  "42" = "Proliferative cells + Myocytes",
  "43" = "Lumbrokinase+ cells",
  "44" = "Neurons",
  "45" = "Secretory cells",
  "46" = "Epidermis",
  "47" = "Epidermis",
  "48" = "Unknown",
  "49" = "Chaetal sacs"
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
  "Unknown" = "#ADD3D0"
)

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
pdf("data/Pleidyi/cluster_plot.pdf", width = 15, height = 15)
print(cluster_plot)
dev.off()

# Identify marker genes for each cell type
Ple_cell_markers <- FindAllMarkers(seob, only.pos = TRUE, group.by = "cell_type", min.pct = 0.10, logfc.threshold = 0.25, assay = "SCT")
write.csv(Ple_cell_markers, file = "data/Pleidyi/Ple_cell_markers.csv", row.names = FALSE)
save(Ple_cell_markers, file = "data/Pleidyi/Ple_cell_markers.rdata")

# Compute average expression per cell type
Ple_average <- AverageExpression(seob, group.by = "cell_type")$SCT
Ple_average <- rownames_to_column(as.data.frame(Ple_average), "GeneID")
save(Ple_average, file = "data/Pleidyi/Ple_average.rdata")

# Save the Seurat object
save(seob, file = "data/Pleidyi/Pleidyi.scRNAseq.rdata")
