library(SeuratDisk)
load(file = "data/Pleidyi/Pleidyi.scRNAseq.rdata")
seob_samap <- seob
DefaultAssay(seob_samap) <- "RNA"
seob_samap[["SCT"]] <- NULL
seob_samap[["RNA3"]] <- as(object = seob_samap[["RNA"]], Class = "Assay")
DefaultAssay(seob_samap) <- "RNA3"
seob_samap[["RNA"]] <- NULL
seob_samap[["RNA"]] <- seob_samap[["RNA3"]]
DefaultAssay(seob_samap) <- "RNA"
seob_samap[["RNA3"]] <- NULL
SaveH5Seurat(seob_samap, filename = "Ple.samap.h5Seurat", overwrite = TRUE)
Convert("Ple.samap.h5Seurat", dest = "h5ad", overwrite = TRUE)
