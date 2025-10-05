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

# Assign developmental stage labels based on sample identity
seob$stages <- ifelse(seob$samples == "T", "Trochophore", "Gastrula")

# Integrate datasets using Harmony to correct batch effects
seob_harmony <- IntegrateLayers(object = seob,
                                normalization.method = "SCT",
                                orig.reduction = "pca",
                                method = HarmonyIntegration,
                                new.reduction = "integrated.harmony")

# Run UMAP on integrated data
seob_harmony <- RunUMAP(seob_harmony, assay = "SCT", dims = 1:30,
                        reduction = "integrated.harmony")

# Save sample-colored UMAP visualization
pdf("data/Cgigas/larva/Cgi_UMAP_harmony.pdf", width = 6, height = 6)
DimPlot(seob_harmony, reduction = "umap", group.by = "samples", label = F)
dev.off()

# Update object and compute nearest-neighbor graph
seob <- seob_harmony
seob <- FindNeighbors(seob, reduction = "integrated.harmony", dims = 1:30)

# Explore clustering resolutions (1 to 7)
resolutions <- seq(1, 7, by=0.5)
seob <- FindClusters(seob,
                     resolution = resolutions,
                     random.seed = 1)

# Visualize cluster relationships across resolutions
clustree_plot <- clustree(seob, prefix = "SCT_snn_res.")
ggsave("data/Cgigas/larva/clustree_plot.pdf", plot = clustree_plot, width = 40, height = 20)

# Perform clustering at the optimal resolution (2.0)
seob <- seob_harmony
seob <- FindNeighbors(seob, reduction = "integrated.harmony", dims = 1:30)
seob <- FindClusters(seob, resolution = 2.0, random.seed = 1)

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
options(future.globals.maxSize = 10 * 1024^3) 
seob <- SCTransform(seob)
seob <- RunPCA(seob, assay = "SCT")
seob <- RunUMAP(seob, reduction = "pca", dims = 1:30)

# Batch correction across developmental stages
seob$stages <- ifelse(seob$samples == "T", "Trochophore", "Gastrula")
seob <- IntegrateLayers(object = seob,
                        normalization.method = "SCT",
                        orig.reduction = "pca",
                        method = HarmonyIntegration,
                        group.by.vars = "stages",
                        new.reduction = "integrated.harmony")

# Re-run UMAP on Harmony-integrated data
seob <- RunUMAP(seob, assay = "SCT", dims = 1:30, reduction = "integrated.harmony")

# Compute neighbors and cluster cells
seob <- FindNeighbors(seob, reduction = "integrated.harmony", dims = 1:30)
seob <- FindClusters(seob, resolution = 2.0, random.seed = 1)

# Plot UMAP with custom colors and cluster labels at centroids
cluster_plot <- DimPlot(seob, reduction = "umap", group.by = "seurat_clusters", label = TRUE) +
  theme_void(base_size = 12) +
  theme(legend.position = "na") +
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
pdf("data/Cgigas/larva/cluster_plot.pdf", width = 15, height = 15)
print(cluster_plot)
dev.off()

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
  "LOC105319624", "LOC105340018", "LOC105344040", "LOC105344230", "LOC105344097",
  "LOC105340330", "LOC105317751", "LOC105340517", "LOC105341857","LOC105337940",
  "LOC105339039", "LOC105341104", "LOC105324772", "LOC105347355", "LOC105326938",
  "LOC105328755", "LOC105319677", "LOC105323511", "LOC105328903", "LOC105339769",
  "LOC105325470", "LOC105328534", "LOC105343071", "LOC105320177", "LOC105331920",
  "LOC105338883", "LOC105338950", "LOC105329634", "LOC105336521", "LOC105348570",
  "LOC105348490", "LOC105335486", "LOC105346273", "LOC105329481", "LOC105335713",
  "LOC105342762", "LOC105323195"
)

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
cluster_order <- rev(c(76, 19, 30, 36, 15, 34, 27, 73, 29, 32, 41, 68, 83, 85, 33, 40, 44, 50, 52, 53, 54, 61, 77, 8, 12, 14, 16, 18, 21, 25, 
26, 28, 31, 35, 45, 49, 56, 58, 59, 66, 70, 79, 80, 81, 60, 0, 7, 47, 48, 63, 2, 4, 6, 9, 11, 22, 38, 39, 62, 65, 69, 74, 75, 78, 86, 1, 10, 
17, 24, 72, 67, 87, 57, 20, 23, 84, 42, 51, 55, 37, 64, 3, 43, 71, 5, 13, 46, 82))
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
  "0" = "Hemocytes",
  "1" = "Myocytes",
  "2" = "Cilia",
  "3" = "Cilia + Myocytes",
  "4" = "Cilia",
  "5" = "Unknown",
  "6" = "Cilia",
  "7" = "Hemocytes",
  "8" = "Proliferative cells",
  "9" = "Cilia",
  "10" = "Myocytes",
  "11" = "Cilia",
  "12" = "Proliferative cells",
  "13" = "Unknown",
  "14" = "Proliferative cells",
  "15" = "Shell gland1",
  "16" = "Proliferative cells",
  "17" = "Myocytes",
  "18" = "Proliferative cells",
  "19" = "Shell gland1",
  "20" = "Proliferative cells + Myocytes",
  "21" = "Proliferative cells",
  "22" = "Cilia",
  "23" = "Proliferative cells + Myocytes",
  "24" = "Myocytes",
  "25" = "Proliferative cells",
  "26" = "Proliferative cells",
  "27" = "Shell gland2",
  "28" = "Proliferative cells",
  "29" = "Shell gland3",
  "30" = "Shell gland1",
  "31" = "Proliferative cells",
  "32" = "Shell gland3",
  "33" = "Shell gland4",
  "34" = "Shell gland1",
  "35" = "Proliferative cells",
  "36" = "Shell gland1",
  "37" = "Proliferative cells + Myocytes",
  "38" = "Cilia",
  "39" = "Cilia",
  "40" = "Shell gland4",
  "41" = "Shell gland3",
  "42" = "Proliferative cells + Myocytes",
  "43" = "Cilia + Myocytes",
  "44" = "Shell gland4",
  "45" = "Proliferative cells",
  "46" = "Unknown",
  "47" = "Hemocytes",
  "48" = "Hemocytes",
  "49" = "Proliferative cells",
  "50" = "Shell gland4",
  "51" = "Proliferative cells + Myocytes",
  "52" = "Shell gland4",
  "53" = "Shell gland4",
  "54" = "Shell gland4",
  "55" = "Proliferative cells + Myocytes",
  "56" = "Proliferative cells",
  "57" = "Hemocytes + Myocytes",
  "58" = "Proliferative cells",
  "59" = "Proliferative cells",
  "60" = "Neurons",
  "61" = "Shell gland4",
  "62" = "Cilia",
  "63" = "Hemocytes",
  "64" = "Neurons + Myocytes",
  "65" = "Cilia",
  "66" = "Proliferative cells",
  "67" = "Myocytes",
  "68" = "Shell gland3",
  "69" = "Cilia",
  "70" = "Proliferative cells",
  "71" = "Cilia + Myocytes",
  "72" = "Myocytes",
  "73" = "Shell gland2",
  "74" = "Cilia",
  "75" = "Cilia",
  "76" = "Shell gland1",
  "77" = "Shell gland4",
  "78" = "Cilia",
  "79" = "Proliferative cells",
  "80" = "Proliferative cells",
  "81" = "Proliferative cells",
  "82" = "Unknown",
  "83" = "Shell gland3",
  "84" = "Proliferative cells + Myocytes",
  "85" = "Shell gland3",
  "86" = "Cilia",
  "87" = "Myocytes"
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
  summarise(UMAP_1 = median(UMAP_1),
            UMAP_2 = median(UMAP_2))

# Define colors for each cell type (no cluster numbers in legend)
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

# Marker gene preparation
Idents(seob) <- "seurat_clusters"
seob <- PrepSCTFindMarkers(seob)

# Identify marker genes for each cell type
Cg_cell_markers <- FindAllMarkers(seob, only.pos = TRUE, group.by = "cell_type", min.pct = 0.10, logfc.threshold = 0.25, assay = "SCT")
write.csv(Cg_cell_markers, file = "data/Cgigas/larva/Cg_cell_markers.csv", row.names = FALSE)
save(Cg_cell_markers, file = "data/Cgigas/larva/Cg_cell_markers.rdata")

# Compute average expression per cell type
Cg_average <- AverageExpression(seob, group.by = "cell_type")$SCT
Cg_average <- rownames_to_column(as.data.frame(Cg_average), "GeneID")
save(Cg_average, file = "data/Cgigas/larva/Cg_average.rdata")

# Save the Seurat object
save(seob, file = "data/Cgigas/larva/OysterLarva.scRNAseq.rdata")

# Select only cells belonging to gastrula
cells_use <- rownames(seob@meta.data[seob@meta.data$sample %in% c("G1", "G2", "G3", "G4"), ])

# Generate UMAP plot colored by annotated cell types
cluster_plot <- DimPlot(
  seob, 
  reduction = "umap", 
  group.by = "cell_type_label",
  label = FALSE, 
  cols = my_colors, 
  cells = cells_use) +
  theme_void(base_size = 12) + 
  theme(legend.position = "right") + 
  ggtitle(NULL)

# Add arrows to indicate UMAP axes
xrange <- ggplot_build(cluster_plot)$layout$panel_params[[1]]$x.range
yrange <- ggplot_build(cluster_plot)$layout$panel_params[[1]]$y.range
arrow_base_x <- xrange[1] + 0.1 * diff(xrange)
arrow_base_y <- yrange[1] + 0.1 * diff(yrange)

# Generate Plot for gastrula cells
cluster_plot <- cluster_plot + 
  annotate("segment", x = arrow_base_x, xend = arrow_base_x + 0.1 * diff(xrange),
           y = arrow_base_y, yend = arrow_base_y, 
           arrow = arrow(length = unit(0.2, "cm")), size = 0.5) +
  annotate("text", x = arrow_base_x + 0.1 * diff(xrange), y = arrow_base_y,
           label = "UMAP1", hjust = 0, vjust = 1, size = 3.5) +
  annotate("segment", x = arrow_base_x, xend = arrow_base_x,
           y = arrow_base_y, yend = arrow_base_y + 0.1 * diff(yrange),
           arrow = arrow(length = unit(0.2, "cm")), size = 0.5) +
  annotate("text", x = arrow_base_x + 0.05 * diff(xrange), 
           y = arrow_base_y + 0.1 * diff(yrange), 
           label = "UMAP2", hjust = 1, vjust = 0, size = 3.5) +
  coord_fixed(ratio = 1)
pdf("data/Cgigas/larva/cluster_plot.gastrula.pdf", width = 15, height = 15)
print(cluster_plot)
dev.off()

# Select only cells belonging to trochophore
cells_use <- rownames(seob@meta.data[seob@meta.data$sample %in% c("T"), ])

# Generate UMAP plot colored by annotated cell types
cluster_plot <- DimPlot(
  seob, 
  reduction = "umap", 
  group.by = "cell_type_label",
  label = FALSE, 
  cols = my_colors, 
  cells = cells_use) +
  theme_void(base_size = 12) + 
  theme(legend.position = "right") + 
  ggtitle(NULL)

# Add arrows to indicate UMAP axes
xrange <- ggplot_build(cluster_plot)$layout$panel_params[[1]]$x.range
yrange <- ggplot_build(cluster_plot)$layout$panel_params[[1]]$y.range
arrow_base_x <- xrange[1] + 0.1 * diff(xrange)
arrow_base_y <- yrange[1] + 0.1 * diff(yrange)

# Generate Plot for trochophore cells
cluster_plot <- cluster_plot + 
  annotate("segment", x = arrow_base_x, xend = arrow_base_x + 0.1 * diff(xrange),
           y = arrow_base_y, yend = arrow_base_y, 
           arrow = arrow(length = unit(0.2, "cm")), size = 0.5) +
  annotate("text", x = arrow_base_x + 0.1 * diff(xrange), y = arrow_base_y,
           label = "UMAP1", hjust = 0, vjust = 1, size = 3.5) +
  annotate("segment", x = arrow_base_x, xend = arrow_base_x,
           y = arrow_base_y, yend = arrow_base_y + 0.1 * diff(yrange),
           arrow = arrow(length = unit(0.2, "cm")), size = 0.5) +
  annotate("text", x = arrow_base_x + 0.05 * diff(xrange), 
           y = arrow_base_y + 0.1 * diff(yrange), 
           label = "UMAP2", hjust = 1, vjust = 0, size = 3.5) +
  coord_fixed(ratio = 1)
pdf("data/Cgigas/larva/cluster_plot.trochophore.pdf", width = 15, height = 15)
print(cluster_plot)
dev.off()
