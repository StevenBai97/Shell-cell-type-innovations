#!/usr/bin/env Rscript

###############################################
# parse parameter 
###############################################
library(argparser, quietly=TRUE)

p <- arg_parser("make OrgDB from emapper")
p <- add_argument(p, "emapper_anno", help="emapper annotation result", type="character")
p <- add_argument(p, "proteins", help="proteins in fasta format", type="character")

argv <- parse_args(p)

# set script dir
script_dir <- dirname(strsplit(commandArgs(trailingOnly = FALSE)[4],"=")[[1]][2])

###############################################
# test input 
###############################################
# argv <- list()
# argv$emapper_anno <- 'out.emapper.annotations'
# argv$proteins <- 'proteins.fa'
# script_dir <- 'emcp'

library(tidyverse, quietly = TRUE)
library(formattable, quietly = TRUE)
library(AnnotationForge, quietly = TRUE)
library(seqinr, quietly = TRUE)
library(clusterProfiler, quietly = TRUE)


###############################################
# parse parameter 
###############################################
emapper <- read_delim(argv$emapper_anno, 
                      "\t", escape_double = FALSE, col_names = FALSE, 
                      comment = "#", trim_ws = TRUE) %>%
  dplyr::select(GID = X1, 
                COG = X7,
                Gene_Name = X8,
                Gene_Symbol = X9,
                GO = X10,
                KO = X12,
                Pathway = X13
                )

###############################################
# make OrgDB 
###############################################
# gene name
gene_info <- dplyr::select(emapper,  GID, Gene_Name) %>%
  dplyr::filter(!is.na(Gene_Name)) %>%
  dplyr::filter(Gene_Name != '-') 
eggnog_anno = length(gene_info$GID)

# gene to gene ontology
gene2go <- dplyr::select(emapper, GID, GO) %>%
  separate_rows(GO, sep = ',', convert = F) %>%
  filter(!is.na(GO)) %>%
  filter(str_detect(GO, '^GO')) %>%
  mutate(EVIDENCE = 'IEA') 
  
# make org package
makeOrgPackage(gene_info=gene_info,
               go=gene2go,
               maintainer='zhangsan <zhangsan@genek.cn>',
               author='zhangsan',
               outputDir="./",
               tax_id=0000,
               genus='M',
               species='y',
               goTable="go",
               version="1.0")

# build package
pkgbuild::build('.//org.My.eg.db', dest_path = ".")

# install package  
dir.create('R_Library', recursive = T)
install.packages('org.My.eg.db_1.0.tar.gz',
                 repos = NULL,
                 lib = 'R_Library')
library(org.My.eg.db, lib = 'R_Library')

###############################################
# GO statistics and plot
###############################################
all_gene <- getName.list(read.fasta(file = argv$proteins, 
                                    seqtype = 'AA'))
go_anno = length(unique(gene2go$GID))
go_bp <- clusterProfiler::groupGO(gene     = all_gene,
                 OrgDb    = org.My.eg.db,
                 keyType  = "GID",
                 ont      = "BP",
                 level    = 2,
                 readable = FALSE)
  
go_bp <- as.data.frame(go_bp)
go_bp$GO_Class <- "Biological Process"
  
go_cc <- clusterProfiler::groupGO(gene     = all_gene,
                 OrgDb    = org.My.eg.db,
                 keyType  = "GID",
                 ont      = "CC",
                 level    = 2,
                 readable = FALSE)
  
go_cc <- as.data.frame(go_cc)
go_cc$GO_Class <- "Cellular Component"
  
go_mf <- clusterProfiler::groupGO(gene     = all_gene,
                 OrgDb    = org.My.eg.db,
                 keyType  = "GID",
                 ont      = "MF",
                 level    = 2,
                 readable = FALSE)
go_mf <- as.data.frame(go_mf)
go_mf$GO_Class <- "Molecular Function"
  
go_all <- rbind(go_bp, go_cc, go_mf)
write.table(go_all, "go.txt", sep = "\t", quote = F)
p <- ggplot(go_all) + 
  geom_bar(aes(x = Description, 
               y = Count,
               fill = GO_Class),
           stat = "identity") + 
  facet_wrap(~GO_Class, scales = "free_x") + 
  labs(title = "GO function classification", y = "Number of genes") +
  theme_classic() +
  theme(plot.title = element_text(hjust = 0.5),
        axis.text.x = element_text(angle = 90, hjust = 1, vjust = 0.4),
        axis.title.x = element_blank(),
        legend.position = "none")
  
ggsave("go.pdf", p, width = 20, height = 7)
  
###############################################
# Pathway statistics and plot 
###############################################
gene2pathway <- dplyr::select(emapper, GID, Pathway) %>%
  separate_rows(Pathway, sep = ',', convert = F) %>%
  filter(str_detect(Pathway, 'ko'))
  
load(file = paste(script_dir, "kegg_info.RData", sep = "/"))
  
gene2pathway <- gene2pathway %>%
  left_join(pathway2name, by = "Pathway") %>%
  dplyr::select(GID, Pathway, Pathway_Name, Pathway_Class, Pathway_Subclass) %>%
  distinct() %>%
  na.omit()
  
pathway_anno = length(unique(gene2pathway$GID))
  
pathway_stat <- dplyr::select(gene2pathway, GID, Pathway_Class, Pathway_Subclass) %>% 
  distinct() %>% 
  group_by(Pathway_Class, Pathway_Subclass) %>%
  summarise(Count = n(), Percentage = percent(n()/pathway_anno))
  
pathway_stat$Pathway_Subclass <- ordered(pathway_stat$Pathway_Subclass, levels = pathway_stat$Pathway_Subclass) 
  
p <- ggplot(pathway_stat, aes(x = Pathway_Subclass, y = Percentage)) +
  geom_bar(aes(fill = Pathway_Class), stat = 'identity') +
  geom_text(aes(label = Count), nudge_y = 0.005) +
  scale_y_continuous(labels=percent) + 
  labs(y = "Percent of genes(%)", x ="", fill = "Class") +
  coord_flip() +
  theme_classic()
  
ggsave("pathway.pdf", p, width = 20, height = 7)
write.table(gene2pathway, file = "pathway.txt", sep = "\t", quote = F, row.names = F)
write.table(pathway2name, file = 'pathway_name.txt', sep = '\t', quote = F,row.names = F)
write.table(pathway_stat, file = "pathway_stat.txt", sep = "\t", quote = F, row.names = F)
  

###############################################
# COG statistics and plot 
###############################################  
cog_funclass <- read_delim(paste(script_dir, "cog_funclass.tab", sep = "/"), 
                             "\t", escape_double = FALSE, trim_ws = TRUE)
  
insert_comma <- function(x){
  str_c(x, sep = '', collapse = ',')
}
  
gene2cog <- dplyr::select(emapper, GID, COG) %>%
  filter(!is.na(COG)) %>%
  filter(COG != '-') %>%
  mutate(COG = sapply(str_split(COG, ''), insert_comma)) %>%
  separate_rows(COG, sep = ',', convert = F) %>%
  left_join(cog_funclass, by = c('COG' = 'COG'))
  
cog_anno = length(unique(gene2cog$GID))
  
write.table(gene2cog, file = "cog.txt", sep = "\t", quote = F, row.names = F)
  
p <- ggplot(data = gene2cog) + 
  geom_bar(aes(x = COG, 
               fill = COG_Name)) +
  labs(title = "COG/KOG Function Classification ", 
       x = "",
       y = "Number of genes") +
  theme_classic() +
  theme(plot.title = element_text(hjust = 0.5),
        legend.title = element_blank(),
        legend.key.size=unit(1,"line"),
        legend.text = element_text(size = 7.5)) +
  guides(fill=guide_legend(ncol=1))
ggsave("cog.pdf", p, width = 16, height = 7)

###############################################
# number and percentage 
###############################################  
total_gene = length(all_gene)
anno_stat <- tibble(
  database = c("EggNOG", "GO", "COG/KOG", "KEGG Pathway"),
  number = comma(c(eggnog_anno, go_anno, cog_anno, pathway_anno), digits = 0),
  percentage = percent(c(eggnog_anno, go_anno, cog_anno, pathway_anno)/total_gene)
)
  
write.table(anno_stat, "anno_stat.txt", quote = F, row.names = F, sep = "\t")
