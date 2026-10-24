library(srtdisk)
load(file = "data/Cgigas/mantle/OysterMantle.rdata")
seob_samap <- seob
DefaultAssay(seob_samap) <- "RNA"
seob_samap[["SCT"]] <- NULL
seob_samap[["RNA3"]] <- as(object = seob_samap[["RNA"]], Class = "Assay")
DefaultAssay(seob_samap) <- "RNA3"
seob_samap[["RNA"]] <- NULL
seob_samap[["RNA"]] <- seob_samap[["RNA3"]]
DefaultAssay(seob_samap) <- "RNA"
seob_samap[["RNA3"]] <- NULL
SaveH5Seurat(seob_samap, filename = "oyster.samap.h5Seurat", overwrite = TRUE)
Convert("oyster.samap.h5Seurat", dest = "h5ad", overwrite = TRUE)

seob$cell_type_label <- gsub("Shell-secreting epithelial cells[1-6]+", "Shell-secreting epithelial cells", seob$cell_type_label)
seob$cell_type_label <- gsub("^Hemocytes.*", "Hemocytes", seob$cell_type_label)
seob_samap <- seob
DefaultAssay(seob_samap) <- "RNA"
seob_samap[["SCT"]] <- NULL
seob_samap[["RNA3"]] <- as(object = seob_samap[["RNA"]], Class = "Assay")
DefaultAssay(seob_samap) <- "RNA3"
seob_samap[["RNA"]] <- NULL
seob_samap[["RNA"]] <- seob_samap[["RNA3"]]
DefaultAssay(seob_samap) <- "RNA"
seob_samap[["RNA3"]] <- NULL
SaveH5Seurat(seob_samap, filename = "oyster.SEC.samap.h5Seurat", overwrite = TRUE)
Convert("oyster.SEC.samap.h5Seurat", dest = "h5ad", overwrite = TRUE)
