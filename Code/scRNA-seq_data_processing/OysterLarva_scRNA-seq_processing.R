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

# Get sample names (directory names under data/Cgigas/larva)
samples = basename(list.dirs("data/Cgigas/larva", recursive = F))

# Loop through samples and create individual Seurat objects
for (sample in samples) {
  scrna_data <- Read10X(
    data.dir = str_c("data/Cgigas/larva/", sample))
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
pdf("data/Cgigas/larva/Cgi.filter_plot.pdf", width = 10, height = 12)
p0 /p2
dev.off()

# Overwrite seob with filtered object
seob <- seob_filter2

# SCTransform normalization (no regression at this stage)
seob_sct_noregress <- SCTransform(seob)

# Read S phase and G2M phase gene IDs from file
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
pdf("data/Cgigas/larva/Cgi_PCA_UMAP.pdf", width = 12, height = 6)
print(p_pca + p_umap)
dev.off()

# Update object and compute nearest-neighbor graph
DefaultAssay(seob) <- "SCT"
seob <- FindNeighbors(seob, reduction = "pca", dims = 1:30)

# Explore clustering resolutions (1 to 7)
resolutions <- seq(1, 8, by=0.1)
seob <- FindClusters(seob,
                     resolution = resolutions,
                     random.seed = 1)

# Visualize cluster relationships across resolutions
clustree_plot <- clustree(seob, prefix = "SCT_snn_res.")
ggsave("data/Cgigas/larva/clustree_plot.pdf", plot = clustree_plot, width = 40, height = 20)

# Perform clustering at the optimal resolution (6.7)
seob <- FindNeighbors(seob, reduction = "pca", dims = 1:30)
seob <- FindClusters(seob, resolution = 6.7, random.seed = 1)

# Plot and save UMAP with cluster labels
cluster_plot <- DimPlot(seob, reduction = "umap", group.by = "seurat_clusters", label = T)
ggsave("data/Cgigas/larva/cluster_plot.pdf", plot = cluster_plot, width = 11, height = 9)

# Doublet removal using DoubletFinder
## 1. Split by sample
seob_list <- SplitObject(seob, split.by = "samples")
options(future.globals.maxSize = 50 * 1024^3)

for (i in 1:length(seob_list)) {
  s <- seob_list[[i]]
  sample_name <- names(seob_list)[i]
  
  ## 2. Estimate expected number of doublets
  nExp_poi <- round(0.04 * nrow(s@meta.data))
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

# Subset Seurat object based on Singlets
seob_singlets <- subset(seob, cells = Singlets)

# Overwrite seob with singlet-only object for downstream analysis
seob <- seob_singlets

# Marker genes identification
Cg_markers <- FindAllMarkers(seob, only.pos = TRUE, min.pct = 0.25, logfc.threshold = 0.25, assay = "SCT")
Cg_markers <- Cg_markers %>%
  arrange(cluster, desc(avg_log2FC))
write.csv(Cg_markers, file = "data/Cgigas/larva/Cg_cluster_markers.csv", row.names = FALSE)

# Define the list of marker genes
ordered_genes <- c(
  "LOC105319624",
  "LOC105344040",
  "LOC105340517",
  "LOC105319410",
  "LOC105344230",
  "LOC105344097",
  "LOC105341857",
  "LOC105329605",
  "LOC105324772",
  "LOC105348917",
  "LOC105345135",
  "LOC105322264",
  "LOC105325470",
  "LOC105328534",
  "LOC105320177",
  "LOC105331920",
  "LOC105338883",
  "LOC105338950",
  "LOC105317061",
  "LOC105329634",
  "LOC105348570",
  "LOC105346273",
  "LOC105336521",
  "LOC105329786",
  "LOC105335486",
  "LOC105348490")

# Extract expression data for DotPlot
marker_data <- DotPlot(seob, 
                   assay = "SCT", 
                   features = ordered_genes,
                   group.by = "seurat_clusters")$data
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
cluster_order <- rev(c(24,13,27,51,57,28,42,11,49,32,50,43,37,40,41,60,12,15,19,34,39,53,17,22,38,54,4,9,18,36,62,1,3,5,8,14,20,25,26,29,30,31,44,56,59,0,2,6,7,10,16,21,23,33,35,45,46,47,48,52,55,58,61))
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
pdf("data/Cgigas/larva/marker_dotplot.pdf", width = 12, height = 15)
print(marker_plot)
dev.off()

# Define mapping from cluster ID to cell type
cluster2type <- c(
  "24" = "Shell gland1",
  "13" = "Shell gland2",
  "27" = "Shell gland2",
  "51" = "Shell gland3",
  "57" = "Shell gland3",
  "28" = "Shell gland4",
  "42" = "Shell gland4",
  "11" = "Proliferative cells",
  "49" = "Proliferative cells",
  "32" = "Proliferative cells",
  "50" = "Proliferative cells",
  "43" = "Neurons",
  "37" = "Hemocytes",
  "40" = "Hemocytes",
  "41" = "Hemocytes",
  "60" = "Hemocytes",
  "12" = "Hemocytes",
  "15" = "Hemocytes",
  "19" = "Hemocytes",
  "34" = "Hemocytes",
  "39" = "Hemocytes",
  "53" = "Hemocytes",
  "17" = "Hemocytes",
  "22" = "Hemocytes",
  "38" = "Hemocytes",
  "54" = "Hemocytes",
  "4" = "Cilia",
  "9" = "Cilia",
  "18" = "Cilia",
  "36" = "Cilia",
  "62" = "Cilia",
  "59" = "Myocytes",
  "56" = "Myocytes",
  "44" = "Myocytes",
  "31" = "Myocytes",
  "30" = "Myocytes",
  "29" = "Myocytes",
  "26" = "Myocytes",
  "25" = "Myocytes",
  "20" = "Myocytes",
  "14" = "Myocytes",
  "8" = "Myocytes",
  "5" = "Myocytes",
  "3" = "Myocytes",
  "1" = "Myocytes",
  "0" = "Unknown",
  "2" = "Unknown",
  "6" = "Unknown",
  "7" = "Unknown",
  "10" = "Unknown",
  "16" = "Unknown",
  "21" = "Unknown",
  "23" = "Unknown",
  "33" = "Unknown",
  "35" = "Unknown",
  "45" = "Unknown",
  "46" = "Unknown",
  "47" = "Unknown",
  "48" = "Unknown",
  "52" = "Unknown",
  "55" = "Unknown",
  "58" = "Unknown",
  "61" = "Unknown")

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
pdf("data/Cgigas/larva/cluster_plot.cell_type.pdf", width = 15, height = 15)
print(cluster_plot)
dev.off()


# Generate DotPlot for cell types
marker_data <- DotPlot(seob, 
                       assay = "SCT", 
                       features = ordered_genes,
                       group.by = "cell_type")$data

cluster_order <- c("Shell gland1","Shell gland2","Shell gland3", "Shell gland4","Proliferative cells",
                       "Neurons","Hemocytes","Cilia","Myocytes","Unknown")

marker_data$id <- factor(marker_data$id, levels = as.character(cluster_order))

marker_plot <- ggplot(marker_data, aes(x =id , y = Gene_name)) +
  geom_point(
    aes(size = pct.exp, color = avg.exp.scaled)) +
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

pdf("data/Cgigas/larva/marker_dotplot.pdf", width = 12, height = 15)
print(marker_plot)
dev.off()

# Identify marker genes for each cell type
Cg_cell_markers <- FindAllMarkers(seob, only.pos = TRUE, group.by = "cell_type", min.pct = 0.25, logfc.threshold = 0.25, assay = "SCT")
write.csv(Cg_cell_markers, file = "data/Cgigas/larva/Cg_cell_markers.csv", row.names = FALSE)
save(Cg_cell_markers, file = "data/Cgigas/larva/Cg_cell_markers.rdata")

# Compute average expression per cell type
Cg_average <- AverageExpression(seob, group.by = "cell_type")$SCT
Cg_average <- rownames_to_column(as.data.frame(Cg_average), "GeneID")
save(Cg_average, file = "data/Cgigas/larva/Cg_average.rdata")
