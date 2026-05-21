# Script Name: Single-Cell RNA-Seq Analysis of Mouse Meninges Immune Cells
# Description: Full analysis pipeline including preprocessing, integration, 
#              clustering, DEG analysis, GO/KEGG enrichment, and visualization.
# Target Journal: Nature Communications

#####----- 0. Package Loading and Initialization -----#####
library(Seurat)
library(dplyr)
library(ggplot2)
library(patchwork)
library(harmony)
library(celldex)
library(SingleR)
library(scales)
library(ggrepel)
library(org.Mm.eg.db)
library(stringr)
library(clusterProfiler)
library(GOplot)
library(tidyr)
library(enrichplot)
library(circlize)
library(ComplexHeatmap)
library(RColorBrewer)




# USER CONFIGURATION: Define File Paths

DIR_YOUNG <- "./your directory"
DIR_OLD   <- "./your directory"
DIR_OLDEV <- "./your directory"

set.seed(1234)

# Pre-computed DEG results file path for GO/KEGG downstream analysis
FILE_DEG_RESULTS <- "./rawdata.csv"


#####----- 1. Data Loading and Seurat Object Creation -----#####

young_data <- Read10X(DIR_YOUNG)
old_data   <- Read10X(DIR_OLD)
oldEV_data <- Read10X(DIR_OLDEV)

young <- CreateSeuratObject(young_data, project = "Young", min.cells = 3, min.features = 200)
young$age <- "young"

old <- CreateSeuratObject(old_data, project = "Old", min.cells = 3, min.features = 200)
old$age <- "old"

oldEV <- CreateSeuratObject(oldEV_data, project = "OldEV", min.cells = 3, min.features = 200)
oldEV$age <- "oldEV"

# Merge datasets
combined <- merge(
  old, y = list(oldEV, young),
  add.cell.ids = c("old", "oldEV", "young"),
  project = "Meninges"
)


#####----- 2. Quality Control and Filtering -----#####

combined$mitoRatio <- PercentageFeatureSet(object = combined, pattern = "^MT-") / 100
combined$log10GenesPerUMI <- log10(combined$nFeature_RNA) / log10(combined$nCount_RNA)

metadata <- combined@meta.data %>%
  dplyr::mutate(seq_folder = orig.ident,
                nUMI = nCount_RNA,
                nGene = nFeature_RNA)
combined@meta.data <- metadata

Meninges <- subset(x = combined, 
                   subset = (nUMI >= 500) & 
                     (nGene >= 200) & 
                     (log10GenesPerUMI > 0.80) & 
                     (mitoRatio < 0.20))


#####----- 3. Global Normalization and Integration (Harmony) -----#####
Meninges <- NormalizeData(Meninges)
Meninges <- FindVariableFeatures(Meninges, nfeatures = 2000)
Meninges <- ScaleData(Meninges)
Meninges <- RunPCA(Meninges)

Meninges <- IntegrateLayers(
  object = Meninges, 
  method = HarmonyIntegration,
  orig.reduction = "pca", 
  new.reduction = "integration",
  verbose = FALSE,
  group.by = "orig.ident" 
)


#####----- 4. Global Clustering and Broad Annotation -----#####
Meninges <- JoinLayers(Meninges)
Meninges <- FindNeighbors(Meninges, reduction = "integration", dims = 1:20)
Meninges <- FindClusters(Meninges, resolution = 0.5) 
Meninges <- RunUMAP(Meninges, reduction = "integration", dims = 1:20)
meninges_markers <- FindAllMarkers(Meninges, only.pos = TRUE, min.pct = 0.25, logfc.threshold = 0.25)
top10_markers <- meninges_markers %>% 
  group_by(cluster) %>% 
  top_n(n = 10, wt = avg_log2FC)
write.csv(top10_markers, file = "Meninges_Cluster_Top10_Markers.csv", row.names = FALSE)


Idents(Meninges) <- "seurat_clusters"
new.cluster.ids.broad <- c(
  "0"  = "B_cells",          "1"  = "Fibroblasts",      "2"  = "Macrophages",
  "3"  = "Neural_Cells",    "4"  = "T_cells",          "5"  = "Endothelial",
  "6"  = "Neutrophils",     "7"  = "Macrophages",      "8"  = "Choroid_Plexus",
  "9"  = "Neural_Cells",    "10" = "Monocytes",        "11" = "Fibroblasts",
  "12" = "Fibroblasts",     "13" = "Proliferating",    "14" = "Fibroblasts",
  "15" = "Dendritic_Cells", "16" = "Fibroblasts",      "17" = "Neural_Cells",
  "18" = "Neural_Cells",    "19" = "Neural_Cells",     "20" = "Neural_Cells",
  "21" = "Others",          "22" = "Mural_Cells",      "23" = "Neutrophils",
  "24" = "Others",          "25" = "Neural_Cells",     "26" = "Neural_Cells",
  "27" = "Choroid_Plexus",  "28" = "Neural_Cells",     "29" = "NK_cells",
  "30" = "T_cells"
)

Meninges <- RenameIdents(Meninges, new.cluster.ids.broad)
Meninges$CellType_Major <- Idents(Meninges)

my_cols_refined <- c(
  "Macrophages"     = "#1F77B4", "Monocytes"       = "#AEC7E8", "Dendritic_Cells" = "#9467BD",
  "Neutrophils"     = "#D62728",
  "T_cells"         = "#2CA02C", "NK_cells"        = "#98DF8A", "B_cells"         = "#FF7F0E",
  "Fibroblasts"     = "#8C564B", "Endothelial"     = "#17BECF", "Mural_Cells"     = "#7F7F7F", 
  "Choroid_Plexus"  = "#BCBD22",
  "Neural_Cells"    = "#E377C2", "Proliferating"   = "#F7B6D2", "Others"          = "#E7E7E7" 
)

umap_all <- DimPlot(Meninges, reduction = "umap", group.by = "CellType_Major", 
                    cols = my_cols_refined, label = TRUE, label.size = 3.5, 
                    label.box = FALSE, repel = TRUE, raster = FALSE) + 
  ggtitle("Figure 1. Single-Cell Atlas of Mouse Meninges") +
  theme(legend.position = "right") 

print(umap_all)
# NOTE: Use ggsave("Figure_1_UMAP.pdf", width = 8, height = 6) to save the output locally.


#####----- 5. Immune Subset Extraction and Fine Annotation -----#####
target_immune <- c("Macrophages", "Monocytes", "Neutrophils", "Dendritic_Cells", "T_cells", "B_cells", "NK_cells")
Meninges_immune <- subset(Meninges, idents = target_immune)

Meninges_immune[["RNA"]] <- split(Meninges_immune[["RNA"]], f = Meninges_immune$orig.ident)
Meninges_immune <- NormalizeData(Meninges_immune)
Meninges_immune <- FindVariableFeatures(Meninges_immune)
Meninges_immune <- ScaleData(Meninges_immune)
Meninges_immune <- RunPCA(Meninges_immune)

Meninges_immune <- IntegrateLayers(
  object = Meninges_immune, method = HarmonyIntegration,
  orig.reduction = "pca", new.reduction = "integration",
  verbose = FALSE, group.by = "orig.ident"
)

Meninges_immune <- JoinLayers(Meninges_immune)
Meninges_immune <- FindNeighbors(Meninges_immune, reduction = "integration", dims = 1:20)
Meninges_immune <- FindClusters(Meninges_immune, resolution = 0.5) 
Meninges_immune <- RunUMAP(Meninges_immune, reduction = "integration", dims = 1:20)

Idents(Meninges_immune) <- "seurat_clusters"
immune_markers <- FindAllMarkers(
  Meninges_immune, 
  only.pos = TRUE,     
  min.pct = 0.25,       
  logfc.threshold = 0.25 
)
top10_immune <- immune_markers %>%
  group_by(cluster) %>%
  top_n(n = 10, wt = avg_log2FC)


Idents(Meninges_immune) <- "seurat_clusters"
cluster_names <- c(
  "0"  = "BAMs",                  "1"  = "B cells",               "2"  = "B cells",
  "3"  = "Macro_Inflammatory",   "4"  = "Active_Neutrophils",    "5"  = "CD8+ T Cells",
  "6"  = "Resident_BAMs",        "7"  = "Monocytes",             "8"  = "Mature_Neutrophils",
  "9"  = "CD4+ T Cells",         "10" = "B cells",               "11" = "B cells",
  "12" = "DCs",                  "13" = "ILC2",                  "14" = "Microglia",
  "15" = "Monocytes",            "16" = "B cells",               "17" = "NK_Cells",
  "18" = "Immature_Neutrophils", "19" = "GammaDelta T Cells",    "20" = "DCs",
  "21" = "B cells",              "22" = "Macro_Inflammatory"
)

Meninges_immune <- RenameIdents(Meninges_immune, cluster_names)
Meninges_immune$CellType_Immune_Fine <- Idents(Meninges_immune) 

umap_immune <- DimPlot(Meninges_immune, reduction = "umap", group.by = "CellType_Immune_Fine", 
                       label = TRUE, label.size = 3.5, label.box = FALSE, repel = TRUE, raster = FALSE) + 
  ggtitle("Figure. Meninges Immune Cell Subtypes") +
  theme(plot.title = element_text(hjust = 0.5, face = "bold", size = 15), legend.position = "right") 

print(umap_immune)


#####----- 6. Immune Cell Composition Bar Plot (Figure 7b) -----#####
label_order <- c("Young", "Old", "OldEV")
target_cells_fig7b <- c("BAMs", "Resident_BAMs", "Immature_Neutrophils", 
                        "Mature_Neutrophils", "Active_Neutrophils", "B cells")

cell_prop_fig7b <- Meninges_immune@meta.data %>%
  group_by(orig.ident, CellType_Immune_Fine) %>%
  summarise(n = n(), .groups = 'drop') %>%
  group_by(orig.ident) %>%
  mutate(percentage = n / sum(n) * 100) %>%
  filter(CellType_Immune_Fine %in% target_cells_fig7b)

cell_prop_fig7b$orig.ident <- factor(cell_prop_fig7b$orig.ident, levels = label_order)
cell_prop_fig7b$CellType_Immune_Fine <- factor(cell_prop_fig7b$CellType_Immune_Fine, levels = target_cells_fig7b)

p_fig7b <- ggplot(cell_prop_fig7b, aes(x = CellType_Immune_Fine, y = percentage, fill = orig.ident)) +
  geom_bar(stat = "identity", position = position_dodge(width = 0.8), width = 0.7) +
  geom_text(aes(label = sprintf("%.1f%%", percentage)), 
            position = position_dodge(width = 0.8), vjust = -0.5, size = 3) +
  scale_fill_manual(values = c("Young" = "#1F77B4", "Old" = "#D62728", "OldEV" = "#2CA02C")) +
  scale_y_continuous(limits = c(0, max(cell_prop_fig7b$percentage) + 5)) +
  labs(x = "", y = "Proportion (% within Group)", fill = "", title = "b") +
  theme_classic() +
  theme(
    axis.text.x = element_text(angle = 45, hjust = 1, face = "bold", color = "black", size = 10),
    axis.text.y = element_text(color = "black", size = 10),
    axis.title.y = element_text(size = 12, face = "bold"),
    legend.position = "top",
    plot.title = element_text(face = "bold", size = 16, hjust = -0.05)
  )

print(p_fig7b)


#####----- 7. Global Immune DEGs Volcano Plot (Old vs OldEV) -----#####
fc_threshold <- 0.25   
pval_threshold <- 0.05 

Idents(Meninges_immune) <- "orig.ident"
deg_results <- FindMarkers(Meninges_immune, ident.1 = "Old", ident.2 = "OldEV", 
                           min.pct = 0.05, logfc.threshold = 0.0, test.use = "wilcox")

deg_results$gene <- rownames(deg_results)
deg_results$log10_fdr <- -log10(deg_results$p_val_adj)

deg_results <- deg_results %>%
  mutate(Significance = case_when(
    p_val_adj < pval_threshold & avg_log2FC > fc_threshold ~ "Up",   
    p_val_adj < pval_threshold & avg_log2FC < -fc_threshold ~ "Down", 
    TRUE ~ "NS"
  ))

top10_up <- deg_results %>% filter(Significance == "Up") %>% arrange(p_val_adj) %>% head(10) %>% pull(gene)
top10_down <- deg_results %>% filter(Significance == "Down") %>% arrange(p_val_adj) %>% head(10) %>% pull(gene)
top_genes_auto <- c(top10_up, top10_down)

deg_results$delabel <- ifelse(deg_results$gene %in% top_genes_auto, deg_results$gene, NA)
total_cells <- ncol(Meninges_immune)

p_volcano <- ggplot(deg_results, aes(x = avg_log2FC, y = log10_fdr)) +
  geom_point(aes(color = Significance), alpha = 0.6, size = 1.5, stroke = 0) +
  scale_color_manual(values = c("Up" = "#B22222", "Down" = "#4169E1", "NS" = "#D3D3D3")) +
  geom_vline(xintercept = c(-fc_threshold, fc_threshold), linetype = "dashed", color = "black", linewidth = 0.8) +
  geom_hline(yintercept = -log10(pval_threshold), linetype = "dashed", color = "black", linewidth = 0.8) +
  geom_text_repel(aes(label = delabel), fontface = "italic", size = 3.5, color = "black",
                  box.padding = 0.5, point.padding = 0.3, segment.color = "grey30", 
                  segment.size = 0.4, min.segment.length = 0, max.overlaps = Inf, na.rm = TRUE) +
  labs(title = "MSC-EV effects on global immune cell",
       subtitle = paste0("Total Immune Population (n = ", total_cells, ")"),
       x = "Log2 Fold Change", y = "-Log10(FDR)") + 
  theme_classic() +
  theme(
    plot.title = element_text(face = "bold", size = 18, hjust = 0.5, color = "black"),
    plot.subtitle = element_text(size = 12, hjust = 0.5, margin = margin(b = 15), color = "black"),
    axis.title.x = element_text(face = "bold", size = 14, margin = margin(t = 10)),
    axis.title.y = element_text(face = "bold", size = 14, margin = margin(r = 10)),
    axis.text = element_text(size = 11, color = "black"),
    axis.line = element_line(linewidth = 1.2, color = "black"),  
    axis.ticks = element_line(linewidth = 1.2, color = "black"), 
    legend.position = "none" 
  )

print(p_volcano)


#####----- 8. Inflammatory Regulation Dot Plot (Figure 7d) -----#####
features_fig7d <- c("S100a8", "S100a9", "Lcn2", "Chil3", "Cxcl13", "Il1b", "Ttr", "Apoe", "Col1a1", "Fn1")
target_cells_fig7d <- c("BAMs", "Active_Neutrophils", "B cells", "Macro_Inflammatory", "Monocytes")

dotplot_order <- c(
  "BAMs - OldEV", "BAMs - Old", "BAMs - Young",
  "Active_Neutrophils - OldEV", "Active_Neutrophils - Old", "Active_Neutrophils - Young",
  "B cells - OldEV", "B cells - Old", "B cells - Young",
  "Macro_Inflammatory - OldEV", "Macro_Inflammatory - Old", "Macro_Inflammatory - Young",
  "Monocytes - OldEV", "Monocytes - Old", "Monocytes - Young"
)

Meninges_immune$CellType_Group <- paste(Meninges_immune$CellType_Immune_Fine, Meninges_immune$orig.ident, sep = " - ")
subset_obj <- subset(Meninges_immune, CellType_Immune_Fine %in% target_cells_fig7d)

dp_base <- DotPlot(subset_obj, features = features_fig7d, group.by = "CellType_Group", scale = FALSE)
plot_data_dp <- dp_base$data

plot_data_dp <- plot_data_dp %>%
  mutate(BaseCellType = sub(" - .*", "", id)) %>%
  group_by(features.plot, BaseCellType) %>%
  mutate(avg.exp.scaled = as.numeric(scale(avg.exp))) %>%
  ungroup()

plot_data_dp$avg.exp.scaled[is.na(plot_data_dp$avg.exp.scaled)] <- 0
plot_data_dp$id <- factor(plot_data_dp$id, levels = rev(dotplot_order))
plot_data_dp$features.plot <- factor(plot_data_dp$features.plot, levels = features_fig7d)

p_fig7d <- ggplot(plot_data_dp, aes(x = features.plot, y = id)) +
  geom_point(aes(size = pct.exp, fill = avg.exp.scaled), shape = 21, color = "black", stroke = 0.8) +
  scale_fill_gradientn(colors = c("#4169E1", "white", "#B22222"), limits = c(-2, 2), oob = scales::squish, breaks = c(-2, -1, 0, 1, 2)) +
  scale_size_continuous(range = c(1, 8), breaks = c(0, 25, 50, 75), labels = c("0.00", "0.25", "0.50", "0.75"), limits = c(0, 100)) +
  labs(title = "Inflammatory regulation of MSC-EVs", x = "", y = "", fill = "Expression\n(Log2FC)", size = "Percent\nExpressed") +
  theme_bw() + 
  theme(
    plot.title = element_text(face = "bold", size = 16, hjust = 0.5, margin = margin(b = 15)),
    axis.text.x = element_text(angle = 45, hjust = 1, face = "italic", color = "black", size = 11),
    axis.text.y = element_text(face = "bold", color = "black", size = 10),
    panel.grid.major = element_line(color = "grey90"), 
    legend.title = element_text(size = 10, face = "bold"),
    legend.text = element_text(size = 9)
  )

print(p_fig7d)


#####----- 9. Target Gene UMAP Feature Plots (Figure 7e) -----#####
Meninges_immune$orig.ident <- factor(Meninges_immune$orig.ident, levels = c("Young", "Old", "OldEV"))

fp_a8 <- FeaturePlot(Meninges_immune, features = "S100a8", split.by = "orig.ident", cols = c("lightgrey", "#D62728"), order = TRUE, pt.size = 0.8, keep.scale = "all")
fp_a9 <- FeaturePlot(Meninges_immune, features = "S100a9", split.by = "orig.ident", cols = c("lightgrey", "#D62728"), order = TRUE, pt.size = 0.8, keep.scale = "all")

apply_pub_theme <- function(p, row_label = "") {
  p <- p + labs(title = NULL, subtitle = NULL, x = "umap_1", y = "umap_2") +
    theme_classic() +
    theme(
      axis.title = element_text(face = "bold", size = 14),
      axis.text = element_text(size = 12, color = "black"),
      axis.line = element_line(linewidth = 1.2, color = "black"),
      axis.ticks = element_line(linewidth = 1.2, color = "black"),
      legend.position = "none" 
    )
  if (row_label != "") {
    p <- p + guides(y.sec = guide_axis(title = row_label)) +
      theme(axis.title.y.right = element_text(face = "bold", size = 16, angle = 270, vjust = 1.5))
  }
  return(p)
}

p1_fp <- apply_pub_theme(fp_a8[[1]]) + labs(title = "Young") + theme(plot.title = element_text(hjust = 0.5, face = "bold", size = 18, margin = margin(b=10)))
p2_fp <- apply_pub_theme(fp_a8[[2]]) + labs(title = "Old") + theme(plot.title = element_text(hjust = 0.5, face = "bold", size = 18, margin = margin(b=10)))
p3_fp <- apply_pub_theme(fp_a8[[3]], row_label = "S100a8") + labs(title = "Old-EV") + theme(plot.title = element_text(hjust = 0.5, face = "bold", size = 18, margin = margin(b=10)))
p4_fp <- apply_pub_theme(fp_a9[[1]])
p5_fp <- apply_pub_theme(fp_a9[[2]])
p6_fp <- apply_pub_theme(fp_a9[[3]], row_label = "S100a9")

p_fig7e <- wrap_plots(p1_fp, p2_fp, p3_fp, p4_fp, p5_fp, p6_fp, ncol = 3) +
  plot_annotation(tag_levels = list(c('e'))) & theme(plot.tag = element_text(size = 22, face = "bold"))

print(p_fig7e)


#####----- 10. GO Biological Process Enrichment Analysis (Figure 7f) -----#####
gene_entrez_go <- bitr(merged_rescue_genes, fromType = "SYMBOL", toType = "ENTREZID", OrgDb = org.Mm.eg.db)

go_merged <- enrichGO(gene = gene_entrez_go$ENTREZID,
                      OrgDb = org.Mm.eg.db,
                      ont = "BP",
                      pAdjustMethod = "BH",
                      pvalueCutoff = 0.05,
                      readable = TRUE)

go_plot_df <- go_merged@result

plot_data <- go_plot_df %>%
  mutate(
    Term_Size = as.numeric(sub("/.*", "", BgRatio)),
    RichFactor = Count / Term_Size,
    Log10Pval = -log10(p.adjust)
  )

keywords_inflammation <- c("response to biotic", "defense response", "response to bacterium", "oxidative stress", "inflammatory")
keywords_regulation   <- c("leukocyte migration", "chemotaxis", "innate immune response", "receptor signaling", "activation of innate", "pattern recognition")
keywords_homeostasis  <- c("homeostasis", "apoptotic", "secretion", "transport", "metabolic", "insulin", "Golgi")

plot_data <- plot_data %>%
  mutate(Category = case_when(
    str_detect(Description, paste(keywords_inflammation, collapse = "|")) ~ "Alleviation of\nChronic Inflammation",
    str_detect(Description, paste(keywords_regulation, collapse = "|")) ~ "Regulation of\nImmune Signaling",
    str_detect(Description, paste(keywords_homeostasis, collapse = "|")) ~ "Restoration of\nHomeostasis",
    TRUE ~ "Others" 
  )) %>%
  filter(Category != "Others" & p.adjust < 0.05) %>%
  group_by(Category) %>%
  arrange(p.adjust) %>%        
  slice_head(n = 5) %>%        
  ungroup()


cat_levels <- c("Alleviation of\nChronic Inflammation", "Regulation of\nImmune Signaling", "Restoration of\nHomeostasis")
plot_data$Category <- factor(plot_data$Category, levels = cat_levels)
plot_data$Description <- factor(plot_data$Description, levels = rev(unique(plot_data$Description)))

plot_data <- plot_data %>%
  mutate(
    Plot_Count = pmax(100, pmin(Count, 175)),
    Plot_RichFactor = pmax(0.325, pmin(RichFactor, 0.450))
  )


x_min <- min(13, min(plot_data$Log10Pval) * 0.85)
x_max <- max(32, max(plot_data$Log10Pval) * 1.05)



make_go_plot <- function(df, cat_name, strip_color, show_x_axis = FALSE) {
  sub_df <- df %>% filter(Category == cat_name)
  
  p <- ggplot(sub_df, aes(x = Log10Pval, y = Description)) +
    geom_point(aes(size = Plot_Count, color = Plot_RichFactor)) +
    facet_grid(Category ~ ., scales = "free_y", space = "free_y") +
    
    scale_color_gradientn(
      colors = c("#2F358F", "#9A1B81", "#E61434"), 
      breaks = c(0.325, 0.350, 0.375, 0.400, 0.425, 0.450),
      name = "Rich Factor",
      guide = guide_colorbar(order = 1, barheight = unit(3, "cm"), ticks.colour = "white")
    ) +
    
    scale_size_continuous(
      range = c(3, 9), 
      breaks = c(100, 125, 150, 175),
      name = "Gene Count",
      guide = guide_legend(order = 2, override.aes = list(color = "black"))
    ) +
    
    scale_x_continuous(breaks = c(15, 20, 25, 30)) +
    coord_cartesian(xlim = c(x_min, x_max)) + 
    
    theme_bw() +
    theme(
      axis.title.y = element_blank(),
      axis.text.y = element_text(color = "black", size = 11, face = "bold"),
      strip.text.y = element_text(angle = 270, face = "bold", size = 11, margin = margin(l = 10, r = 10)),
      strip.background = element_rect(fill = strip_color, color = "black", linewidth = 1.2),
      panel.border = element_rect(color = "black", linewidth = 1.2),
      panel.grid.minor = element_blank()
    )
  
  if (!show_x_axis) {
    p <- p + theme(axis.title.x = element_blank(), axis.text.x = element_blank(), axis.ticks.x = element_blank())
  } else {
    p <- p + labs(x = bquote(-log[10]("Adjusted P-value"))) +
      theme(axis.title.x = element_text(face = "bold", size = 12, margin = margin(t = 10)),
            axis.text.x = element_text(color = "black", size = 10))
  }
  return(p)
}

p1 <- make_go_plot(plot_data, "Alleviation of\nChronic Inflammation", "#92C5DE", show_x_axis = FALSE)
p2 <- make_go_plot(plot_data, "Regulation of\nImmune Signaling", "#A1D99B", show_x_axis = FALSE)     
p3 <- make_go_plot(plot_data, "Restoration of\nHomeostasis", "#C994C7", show_x_axis = TRUE)          

p_fig7f <- (p1 / p2 / p3) + 
  plot_layout(guides = 'collect') + 
  plot_annotation(
    title = "GO : Biological Process of MSC-EVs treatment",
    theme = theme(plot.title = element_text(hjust = 0.5, face = "bold", size = 18))
  ) & 
  theme(
    legend.position = "right",
    legend.box = "vertical", 
    legend.title = element_text(size = 11, face = "bold"),
    legend.text = element_text(size = 10)
  )

print(p_fig7f)



#####----- 11. KEGG Pathway Chord Diagram (Figure 7g) -----#####
gene_entrez_kegg <- bitr(merged_rescue_genes, fromType = "SYMBOL", toType = "ENTREZID", OrgDb = org.Mm.eg.db)

kegg_merged <- enrichKEGG(gene = gene_entrez_kegg$ENTREZID, organism = 'mmu', pvalueCutoff = 0.05)
kegg_readable <- setReadable(kegg_merged, OrgDb = org.Mm.eg.db, keyType = "ENTREZID")
kegg_data <- kegg_readable@result

target_pathways <- c(
  "Autophagy - animal", "TNF signaling pathway", "Chemokine signaling pathway",
  "Mitophagy - animal", "Cellular senescence", "NF-kappa B signaling pathway"
)

link_data_full <- kegg_data %>% filter(Description %in% target_pathways) %>%
  select(Term = Description, Genes = geneID) %>%
  mutate(Genes = str_replace_all(Genes, "/", ", ")) %>%
  separate_rows(Genes, sep = ", ") %>% select(Gene = Genes, Term = Term)

gene_fc_info <- deg_data %>% filter(Gene %in% link_data_full$Gene & Comparison == "Old_vs_OldEV") %>%
  group_by(Gene) %>% summarise(logFC = avg_log2FC[which.max(abs(avg_log2FC))])

fc_vec <- gene_fc_info$logFC
names(fc_vec) <- gene_fc_info$Gene

top_genes_list <- link_data_full %>% inner_join(gene_fc_info, by = "Gene") %>%
  group_by(Term) %>% arrange(desc(abs(logFC))) %>% slice_head(n = 8) %>% pull(Gene) %>% unique()

link_data <- link_data_full %>% filter(Gene %in% top_genes_list)

desired_pathway_order <- c(
  "Autophagy - animal", "TNF signaling pathway", "Chemokine signaling pathway",
  "Mitophagy - animal", "Cellular senescence", "NF-kappa B signaling pathway"
)

link_data$Term <- factor(link_data$Term, levels = desired_pathway_order)
link_data <- link_data %>% arrange(Term)

pathway_cols <- c(
  "NF-kappa B signaling pathway" = "#E41A1C", "TNF signaling pathway"        = "#FF7F00",
  "Chemokine signaling pathway"  = "#984EA3", "Cellular senescence"          = "#FFFF33",
  "Autophagy - animal"           = "#377EB8", "Mitophagy - animal"           = "#4DAF4A" 
)

existing_pathways <- intersect(desired_pathway_order, unique(link_data$Term))
term_colors_vec <- pathway_cols[existing_pathways]

col_fun_gene <- colorRamp2(c(-2, 0, 2), c("#313695", "white", "#A50026"))
valid_genes <- unique(link_data$Gene)
gene_colors_vec <- col_fun_gene(fc_vec[valid_genes])

grid_col <- c(gene_colors_vec, term_colors_vec)
link_colors <- term_colors_vec[link_data$Term]

circos.clear()
n_genes <- length(valid_genes)
n_terms <- length(existing_pathways)

gene_gaps <- if(n_genes > 1) rep(1, n_genes - 1) else c()
term_gaps <- if(n_terms > 1) rep(1, n_terms - 1) else c()
circos.par(gap.after = c(gene_gaps, 15, term_gaps, 15))

chordDiagram(link_data, grid.col = grid_col, col = link_colors, transparency = 0.4, annotationTrack = "grid", preAllocateTracks = 1)

circos.trackPlotRegion(track.index = 1, panel.fun = function(x, y) {
  xlim = get.cell.meta.data("xlim"); ylim = get.cell.meta.data("ylim"); sector.name = get.cell.meta.data("sector.index")
  if(sector.name %in% existing_pathways) {
    circos.text(mean(xlim), ylim[1] + 0.5, sector.name, facing = "inside", niceFacing = TRUE, adj = c(0.5, 0), cex = 0.9, font = 2)
  } else {
    circos.text(mean(xlim), ylim[1] + 0.3, sector.name, facing = "clockwise", niceFacing = TRUE, adj = c(0, 0.5), cex = 0.8, font = 3)
  }
}, bg.border = NA)

title("Mechanisms of Rejuvenation by MSC-EV", cex.main = 1.5)
lgd_gene = Legend(title = "Log2FC", col_fun = col_fun_gene, at = c(-2, 0, 2))
draw(lgd_gene, x = unit(0.9, "npc"), y = unit(0.2, "npc"), just = c("center", "center"))