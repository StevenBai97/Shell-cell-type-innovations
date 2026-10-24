# Load packages
library(SeuratWrappers)
library(monocle3)
library(ggplot2)

# Load Seurat object and subset unknown cell types
load(file = "OysterMantle.scRNAseq.rdata")
seob <- subset(seob, subset = cell_type != "Unknown")
cds <- as.cell_data_set(seob,
                        assay = "SCT",
                        reductions = c("pca", "umap"),
                        default.reduction = "umap")

# Assign gene short names for plotting
rowData(cds) <- data.frame(
  row.names = rownames(cds),
  gene_short_name = rownames(cds)
)

# Keep only relevant shell-secreting epithelial cells and proliferative cells
keep_types <- c("Shell-secreting epithelial cells1","Shell-secreting epithelial cells2",
                "Shell-secreting epithelial cells3","Shell-secreting epithelial cells4",
                "Shell-secreting epithelial cells5","Shell-secreting epithelial cells6",
                "Proliferative cells")
cds <- cds[, colData(cds)$cell_type %in% keep_types]

# Cluster cells in UMAP space
cds <- cluster_cells(cds, reduction_method = "UMAP")

# Plot UMAP colored by cell type
pdf("plot_cells.pdf", width = 6, height = 6)
plot_cells(cds, label_groups_by_cluster = FALSE, color_cells_by = "cell_type")
dev.off()

# Learn trajectory graph
cds <- learn_graph(cds, use_partition = F, close_loop = F)

# Function to select root node based on a known early cell type
get_earliest_principal_node <- function(cds, cell_type){
  cell_ids <- which(colData(cds)[, "cell_type"] == cell_type)
  
  closest_vertex <-
    cds@principal_graph_aux[["UMAP"]]$pr_graph_cell_proj_closest_vertex
  closest_vertex <- as.matrix(closest_vertex[colnames(cds), ])
  root_pr_nodes <-
    igraph::V(principal_graph(cds)[["UMAP"]])$name[as.numeric(names
                                                              (which.max(table(closest_vertex[cell_ids,]))))]
  
  root_pr_nodes
}

# Order cells along pseudotime using proliferative cells as the root
cds <- order_cells(cds, root_pr_nodes = get_earliest_principal_node(cds, "Proliferative cells"))

# Plot UMAP with cell type labels and trajectory
pdf("plot_cells.cell_type.pdf", width = 6, height = 6)
plot_cells(cds,
           reduction_method = "UMAP",
           color_cells_by = "cell_type",
           cell_size = 0.8,
           label_roots = FALSE, 
           show_trajectory_graph = TRUE, 
           trajectory_graph_color = "black",
           trajectory_graph_segment_size = 1,
           label_cell_groups = TRUE,
           label_branch_points = FALSE,
           label_leaves = FALSE)
dev.off()

# Plot UMAP colored by pseudotime
pdf("plot_cells.pseudotime.pdf", width = 7.1, height = 6)
plot_cells(cds,
           reduction_method = "UMAP",
           color_cells_by = "pseudotime",
           cell_size = 0.8,
           label_roots = FALSE, 
           show_trajectory_graph = TRUE,
           trajectory_graph_color = "black",
           trajectory_graph_segment_size = 1,
           label_cell_groups = FALSE,
           label_branch_points = FALSE,
           label_leaves = FALSE)
dev.off()
