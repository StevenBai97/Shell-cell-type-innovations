library(Seurat)
library(tidyverse)
library(dplyr)
library(SeuratDisk)

obj <- readRDS("data/Smed_ACME.rds")
seob <- UpdateSeuratObject(obj)
cell_type_map <- read.table("data/Smediterranea/cell_type", header = TRUE, sep = "\t", stringsAsFactors = FALSE)
seob$cluster <- as.numeric(as.character(seob$RNA_snn_res.1.6))
seob@meta.data <- seob@meta.data %>%
  left_join(cell_type_map, by = c("cluster" = "number"))
colnames(seob@meta.data)[colnames(seob@meta.data) == "name"] <- "cell_type"
colnames(seob@meta.data)[colnames(seob@meta.data) == "group"] <- "cell_type_label"

seob_samap <- seob
seob_samap@assays$RNA@var.features <- character(0)
seob_samap@assays$RNA@scale.data <- matrix(nrow = 0, ncol = 0)
SaveH5Seurat(seob_samap, filename = "data/Smediterranea/Sme.samap.h5Seurat", overwrite = TRUE)
Convert("data/Smediterranea/Sme.samap.h5Seurat", dest = "h5ad", overwrite = TRUE)
