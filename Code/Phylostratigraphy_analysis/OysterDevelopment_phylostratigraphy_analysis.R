library(myTAI)
library(ggplot2)
library(dplyr)

gene_ages <- read.table("gene_ages.tsv", header = TRUE, row.names = 1, sep = "\t")
Cgi_exp <- read.delim("salmon.gene.TMM.EXPR.matrix", header = TRUE, row.names = 1)
y_order <- c("TC", "FC", "EM", "M", "B", "RM", "FS", "EG", "G", "T1", "T2", "T3", "T4", "T5",
             "ED1", "ED2", "D1", "D2", "D3", "D4", "D5", "D6", "D7", "EU1", "EU2",
             "U1", "U2", "U3", "U4", "U5", "U6", "LU1", "LU2", "P1", "P2", "S", "J")
Cgi_exp <- Cgi_exp[, y_order]
merged_data <- merge(gene_ages, Cgi_exp, by = "row.names", all = FALSE)
colnames(merged_data)[1:3] <- c("GeneID", "Other", "Phylostratum")

dataset <- merged_data %>%
  dplyr::select(Phylostratum, GeneID, all_of(y_order))

tai_values <- TAI(dataset)

pdf("TAI.dev.pdf", width = 12, height = 8)
PlotSignature(
  ExpressionSet = dataset,
  permutations = 20000,
  measure       = "TAI",
  TestStatistic = "FlatLineTest",
  xlab          = "Ontogeny",
  ylab          = "TAI",
  lwd           = 2
  )
dev.off()

partialTAI <- pMatrix(dataset)
write.csv(partialTAI,"partialTAI.dev.csv")

flatline_test <- FlatLineTest(ExpressionSet = dataset, permutations = 20000, plotHistogram = FALSE)

tai_df <- data.frame(
  Stage = factor(colnames(dataset)[3:ncol(dataset)], levels = colnames(dataset)[3:ncol(dataset)]),
  TAI = TAI(dataset),
  SD = flatline_test$std.dev
)

ggplot(tai_df, aes(x = Stage, y = TAI, group = 1)) +
  geom_ribbon(aes(ymin = TAI - SD, ymax = TAI + SD), fill = "#0571b0", alpha = 0.2) +
  geom_line(color = "#0571b0", size = 1.5) +
  theme_bw() +
  theme(
    axis.text = element_text(color = "black"),
    axis.title = element_text(color = "black"),
    axis.ticks = element_line(color = "black"),
    plot.title = element_text(color = "black", hjust = 0.5)
  ) +
  labs(x = "Stage", y = "TAI") +
  scale_y_continuous(
  expand = expansion(mult = c(0.05, 0.05)),
  breaks = seq(floor(min(tai_df$TAI - tai_df$SD)*10)/10,
               ceiling(max(tai_df$TAI + tai_df$SD)*10)/10,
               by = 0.5)
  )

ggsave("TAI.dev.ggplot.pdf", width = 12, height = 6)
