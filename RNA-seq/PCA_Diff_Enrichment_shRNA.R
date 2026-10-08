setwd('/data/wangyy/huxiao/RNA_seq/reanalysis/shRNA/04_diff')

library(DESeq2)

count <- read.csv('shRNA.csv')
rownames(count) <- sub("\\..*$", "", count$Geneid)
count_filtered <- count[rowSums(count[, 5:ncol(count), drop = FALSE] == 0 ) == 0 , ]
count_filtered <- count_filtered[, 5:ncol(count)]
count_filtered$shRNA_Me2_PBS_2 <- NULL
metadata <- data.frame(
  Interference = factor(c(rep("shMe2",5), rep("shNC",6))),
  Treatment = factor(c(rep("PBS",2), rep("Silica",3), rep("PBS",3), rep("Silica",3))))
rownames(metadata) <- colnames(count_filtered)

# PCA plot
dds <- DESeqDataSetFromMatrix(countData = count_filtered,
                              colData = metadata,
                              design = ~ Interference + Treatment)
rld <- rlog(dds, blind=TRUE)
# plotPCA(rld, intgroup="Sample")
rld_mat <- assay(rld)

library(ggbiplot)
library(ggforce)
groups <- c(rep('shMe2_PBS', 2),
            rep('shMe2_Silica', 3),
            rep('shNC_PBS', 3),
            rep('shNC_Silica', 3))
pca <- prcomp(t(rld_mat), scale = T)
coordinate <- data.frame(
  PC1 = pca$x[, 1],
  PC2 = pca$x[, 2],
  Group = groups
)
cols  <- c("#5B1A6E", "#2C9C9C", "#F1D32B", "#E64B35")
fills <- c("#C8B8D0", "#B9D8DB", "#F6EDB5", "#F7C6C2")

ggbiplot(
  pca,
  obs.scale = 1,
  var.scale = 1,
  groups = groups,
  circle = FALSE,
  var.axes = FALSE
) +
  geom_mark_ellipse(
    data = coordinate,
    aes(x = PC1, y = PC2, group = groups, fill = groups),
    alpha = 0.20,
    expand = unit(5, "mm"),
    color = NA,
    show.legend = FALSE
  ) +
  scale_shape_manual(values = rep(16, 4), guide = "none") + 
  scale_color_manual(values = cols, name = "Group") +
  scale_fill_manual(values = fills, guide = "none") +
  labs(
    title = "shRNA",
    color = "Group"
  ) +
  scale_x_continuous(
    limits = c(-170, 170),
    breaks = seq(-150, 150, 50)
  ) +
  scale_y_continuous(
    limits = c(-120, 120),
    breaks = seq(-150, 150, 50)
  ) +
  theme_classic(base_size = 14) +
  theme(
    plot.title = element_text(hjust = 0.5, face = "bold", size = 14),
    legend.title = element_text(face = "bold"),
    legend.text = element_text(size = 11),
    axis.title = element_text(size = 12),
    axis.text = element_text(size = 11),
    plot.margin = margin(10, 10, 10, 10)
  )

# Differential analysis
# shNC_PBS vs shNC_Silica
countdata1 <- count_filtered[, c(6:11)]
metadata1 <- metadata[c(6:11), ]
dds1 <- DESeqDataSetFromMatrix(countData = countdata1, colData = metadata1, design = ~ Treatment)
dds1$Treatment <- relevel(dds1$Treatment, ref = "PBS")
dds1 <- DESeq(dds1)
res1 <- lfcShrink(dds1, coef = "Treatment_Silica_vs_PBS")
summary(res1, alpha = 0.05)
res1_df <- data.frame(res1)
sig_res1 <- dplyr::filter(res1_df, padj < 0.05 & abs(log2FoldChange) > 0.263)
up_res1 <- dplyr::filter(res1_df, padj < 0.05 & log2FoldChange > 0.263)
down_res1 <- dplyr::filter(res1_df, padj < 0.05 & log2FoldChange < -0.263)

res1_df$ENSEMBL <- rownames(res1_df)
res1_df$SYMBOL <- unname(
  AnnotationDbi::mapIds(
    org.Mm.eg.db,
    keys = rownames(res1_df),
    keytype = "ENSEMBL",
    column = "SYMBOL",
    multiVals = "first"
  )
)

res1_df$GENENAME <- unname(
  AnnotationDbi::mapIds(
    org.Mm.eg.db,
    keys = rownames(res1_df),
    keytype = "ENSEMBL",
    column = "GENENAME",
    multiVals = "first"
  )
)

write.csv(res1_df, "Silica_induced_DEG_shRNA.csv")

library(ggplot2)
library(dplyr)
library(grid)

volcano1_df <- res1_df %>%
  dplyr::select(log2FoldChange, padj) %>%
  mutate(neglog10_padj = -log10(padj))

volcano1_df <- volcano1_df %>%
  mutate(group = case_when(log2FoldChange < -0.263 & neglog10_padj > 1.3 ~ "DOWN",
                           log2FoldChange >  0.263 & neglog10_padj > 1.3 ~ "UP",
                           TRUE                                          ~ "NS"))

make_grad_raster <- function(col_start, col_end, n = 800, h = 30) {
  cols <- colorRampPalette(c(col_start, col_end))(n)
  as.raster(matrix(cols, nrow = h, ncol = n, byrow = TRUE))
}

left_arrow_raster <- make_grad_raster(col_start = "#6670A0", col_end   = "#E9ECF5", n = 1000, h = 30)
right_arrow_raster <- make_grad_raster(col_start = "#F3E7E7", col_end   = "#A24F4F", n = 1000, h = 30)

volcano1 <- ggplot(volcano1_df, aes(x = log2FoldChange, y = neglog10_padj)) +
  geom_point(data = subset(volcano1_df, group == "NS"), color = "#CFCFCF", alpha = 0.55, size = 1.2) +
  geom_point(data = subset(volcano1_df, group == "DOWN"), color = "#606A97", alpha = 0.70, size = 1.2) +
  geom_point(data = subset(volcano1_df, group == "UP"), color = "#8E3737", alpha = 0.68, size = 1.2) +
  
  scale_x_continuous(limits = c(-4, 4), breaks = c(-3, -2, -1, 0, 1, 2, 3), expand = c(0, 0)) +
  scale_y_continuous(breaks = c(0, 10, 20, 30, 40), expand = c(0, 0)) +
  
  labs(x = expression(log[2]*"[FC]"), y = expression(-log[10]*"[q]")) +
  theme_classic(base_size = 16) +
  theme(axis.line = element_line(linewidth = 0.8,colour = "black"),
        axis.ticks = element_line(linewidth = 0.7, colour = "black"),
        axis.title.x = element_text(margin = margin(t = 18)),
        axis.title.y = element_text(margin = margin(r = 10)),
        legend.position = "none",
        plot.margin = margin(t = 15, r = 15, b = 40, l = 15)) +
  
  coord_cartesian(ylim = c(-1, 50), clip = "off") +
  
  annotation_raster(left_arrow_raster, xmin = -2.40, xmax = -1.15, ymin = -7.40, ymax = -6.60, interpolate = TRUE) +
  annotate("polygon", x = c(-2.35, -2.60, -2.35), y = c(-6.1, -7.0, -7.9), fill = "#6670A0", colour = NA) +
  annotation_raster(right_arrow_raster, xmin = 1.15, xmax = 2.40, ymin =-7.40, ymax = -6.60, interpolate = TRUE) +
  annotate("polygon", x = c(2.35, 2.60, 2.35), y = c(-6.1, -7.0, -7.9), fill = "#A24F4F", colour = NA)+
  
  annotate("text", x = -3.3, y = -7.0, label = "PBS", size = 5) +
  annotate("text", x = 3.3 ,y = -7.0, label = "Silica", size = 5) 

volcano1

# shNC_Silica vs shMe2_Silica
countdata2 <- count_filtered[, c(3:5, 9:11)]
metadata2 <- metadata[c(3:5, 9:11), ]
dds2 <- DESeqDataSetFromMatrix(countData = countdata2, colData = metadata2, design = ~ Interference)
dds2$Interference <- relevel(dds2$Interference, ref = "shNC")
dds2 <- DESeq(dds2)
res2 <- lfcShrink(dds2, coef = "Interference_shMe2_vs_shNC")
summary(res2, alpha = 0.05)
res2_df <- data.frame(res2)
sig_res2 <- dplyr::filter(res2_df, padj < 0.05 & abs(log2FoldChange) > 0.263)
up_res2 <- dplyr::filter(res2_df, padj < 0.05 & log2FoldChange > 0.263)
down_res2 <- dplyr::filter(res2_df, padj < 0.05 & log2FoldChange < -0.263)

res2_df$ENSEMBL <- rownames(res2_df)
res2_df$SYMBOL <- unname(
  AnnotationDbi::mapIds(
    org.Mm.eg.db,
    keys = rownames(res2_df),
    keytype = "ENSEMBL",
    column = "SYMBOL",
    multiVals = "first"
  )
)

res2_df$GENENAME <- unname(
  AnnotationDbi::mapIds(
    org.Mm.eg.db,
    keys = rownames(res2_df),
    keytype = "ENSEMBL",
    column = "GENENAME",
    multiVals = "first"
  )
)

write.csv(res2_df, "shMe2_interference_DEG_shRNA.csv")

volcano2_df <- res2_df %>%
  dplyr::select(log2FoldChange, padj) %>%
  mutate(neglog10_padj = -log10(padj))

volcano2_df <- volcano2_df %>%
  mutate(group = case_when(log2FoldChange < -0.263 & neglog10_padj > 1.3 ~ "DOWN",
                           log2FoldChange >  0.263 & neglog10_padj > 1.3 ~ "UP",
                           TRUE                                          ~ "NS"))

volcano2 <- ggplot(volcano2_df, aes(x = log2FoldChange, y = neglog10_padj)) +
  geom_point(data = subset(volcano2_df, group == "NS"), color = "#CFCFCF", alpha = 0.55, size = 1.2) +
  geom_point(data = subset(volcano2_df, group == "DOWN"), color = "#606A97", alpha = 0.70, size = 1.2) +
  geom_point(data = subset(volcano2_df, group == "UP"), color = "#8E3737", alpha = 0.68, size = 1.2) +
  
  scale_x_continuous(limits = c(-4, 4), breaks = c(-3, -2, -1, 0, 1, 2, 3), expand = c(0, 0)) +
  scale_y_continuous(breaks = c(0, 5, 10, 15, 20), expand = c(0, 0)) +
  
  labs(x = expression(log[2]*"[FC]"), y = expression(-log[10]*"[q]")) +
  theme_classic(base_size = 16) +
  theme(axis.line = element_line(linewidth = 0.8,colour = "black"),
        axis.ticks = element_line(linewidth = 0.7, colour = "black"),
        axis.title.x = element_text(margin = margin(t = 18)),
        axis.title.y = element_text(margin = margin(r = 10)),
        legend.position = "none",
        plot.margin = margin(t = 15, r = 15, b = 40, l = 15)) +
  
  coord_cartesian(ylim = c(-1, 25), clip = "off") +
  
  annotation_raster(left_arrow_raster, xmin = -2.40, xmax = -1.15, ymin = -4.50, ymax = -4.10, interpolate = TRUE) +
  annotate("polygon", x = c(-2.35, -2.60, -2.35), y = c(-3.9, -4.3, -4.8), fill = "#6670A0", colour = NA) +
  annotation_raster(right_arrow_raster, xmin = 1.15, xmax = 2.40, ymin =-4.50, ymax = -4.10, interpolate = TRUE) +
  annotate("polygon", x = c(2.35, 2.60, 2.35), y = c(-3.9, -4.3, -4.8), fill = "#A24F4F", colour = NA)+
  
  annotate("text", x = -3.3, y = -4.3, label = "shNC", size = 5) +
  annotate("text", x = 3.3 ,y = -4.3, label = "shMe2", size = 5) 

volcano2

# Enrichment analysis
library(clusterProfiler)
library(enrichplot)
library(org.Mm.eg.db)
up1_down2 <- intersect(rownames(up_res1), rownames(down_res2))
down1_up2 <- intersect(rownames(down_res1), rownames(up_res2))
gene <- union(up1_down2, down1_up2)

library(VennDiagram)
library(grid)

grid.newpage()

pushViewport(viewport(x = 0.50, y = 0.75, width  = unit(3.2, "inches"), height = unit(3.2, "inches")))

venn1 <- draw.pairwise.venn(
  area1 = 200,
  area2 = 200,
  cross.area = 50,
  category = c("", ""),
  scaled = TRUE,
  euler.d = TRUE,
  fill = c("#6e48fb", "#ffa500"),
  alpha = c(0.5, 0.5),
  col = c("black", "black"),
  lwd = 1.2,
  lty = 1,
  cex = 0,
  cat.cex = 0,
  ind = FALSE)

grid.draw(venn1)

grid.text("1084", x = unit(0.20, "npc"), y = unit(0.50, "npc"), gp = gpar(fontsize = 13))
grid.text("221", x = unit(0.50, "npc"), y = unit(0.50, "npc"), gp = gpar(fontsize = 13))
grid.text("165",x = unit(0.80, "npc"),y = unit(0.50, "npc"),gp = gpar(fontsize = 13))
popViewport()


pushViewport(viewport(x = 0.50, y = 0.25, width  = unit(3.2, "inches"), height = unit(3.2, "inches")))

venn2 <- draw.pairwise.venn(
  area1 = 200,
  area2 = 200,
  cross.area = 50,
  category = c("", ""),
  scaled = TRUE,
  euler.d = TRUE,
  fill = c("#ff7d82", "#7f7f7f"),
  alpha = c(0.5, 0.5),
  col = c("black", "black"),
  lwd = 1.2,
  lty = 1,
  cex = 0,
  cat.cex = 0,
  ind = FALSE)

grid.draw(venn2)
grid.text("797", x = unit(0.20, "npc"), y = unit(0.50, "npc"), gp = gpar(fontsize = 13))
grid.text("44", x = unit(0.50, "npc"), y = unit(0.50, "npc"), gp = gpar(fontsize = 13))
grid.text("108", x = unit(0.80, "npc"), y = unit(0.50, "npc"), gp = gpar(fontsize = 13))
popViewport()


grid.text(
  "Silica Induced\nUP (1305)", 
  x = unit(0.17, "npc"), 
  y = unit(0.75, "npc"), 
  just = "centre", 
  gp = gpar(fontsize = 12, col = "#9999EE", lineheight = 1.05))

grid.text(
  "Me2 Knockout\nDOWN (386)",
  x = unit(0.83, "npc"),
  y = unit(0.75, "npc"),
  just = "centre",
  gp = gpar(fontsize = 12, col = "#FFA500", lineheight = 1.05))

grid.text(
  "Silica Induced\nDOWN (841)",
  x = unit(0.17, "npc"),
  y = unit(0.25, "npc"),
  just = "centre",
  gp = gpar(fontsize = 12, col = "#FF6666", lineheight = 1.05))

grid.text(
  "Me2 Knockout\nUP (152)",
  x = unit(0.83, "npc"),
  y = unit(0.25, "npc"),
  just = "centre",
  gp = gpar(fontsize = 12, col = "#777777", lineheight = 1.05))

grid.text(
  "Enrichment",
  x = unit(0.50, "npc"),
  y = unit(0.50, "npc"),
  just = "centre",
  gp = gpar(fontsize = 16, col = "#000000", lineheight = 1.05))


grid.segments(
  x0 = unit(0.50, "npc"), y0 = unit(0.695, "npc"),
  x1 = unit(0.50, "npc"), y1 = unit(0.545, "npc"), 
  gp = gpar(col = "black", lwd = 2.0, lineend = "round"),
  arrow = arrow(length = unit(0.14, "inches"), type = "open", ends = "last"))

grid.segments(
  x0 = unit(0.50, "npc"), y0 = unit(0.305, "npc"),   
  x1 = unit(0.50, "npc"), y1 = unit(0.455, "npc"),
  gp = gpar(col = "black", lwd = 2.0, lineend = "round"),
  arrow = arrow(length = unit(0.14, "inches"), type = "open", ends = "last"))


GO_BP <- enrichGO(gene = gene,  
                  OrgDb = "org.Mm.eg.db",  
                  keyType = 'ENSEMBL',  
                  ont = 'BP', 
                  pAdjustMethod = 'BH',  
                  pvalueCutoff = 0.05,  
                  qvalueCutoff = 0.25,  
                  readable = FALSE)
GO_BP <- GO_BP[,]
write.csv(GO_BP, "Enrichment_shRNA.csv")
