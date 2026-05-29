#####----- 0. Package Loading and Initialization -----#####
library(readxl)
library(dplyr)
library(tidyr)
library(clusterProfiler)
library(org.Hs.eg.db)
library(ggplot2)
library(circlize)
library(tidyverse)
library(ggalluvial)

#####----- 1. Data Loading -----#####
file_path <- "C:/Users/User/Desktop/EV proteomic/EV proteomic/EV proteomic Hek vs MSC.xlsx"

raw_data <- read_excel(file_path, sheet = 1)

#####----- 2. Reading data -----#####
data_unique <- raw_data %>% distinct(Gene, .keep_all = TRUE)
hek_cols <- grep("HEK", names(data_unique), value = TRUE)
msc_cols <- grep("MSC", names(data_unique), value = TRUE)

data_unique$HEK_mean <- rowMeans(data_unique[, hek_cols], na.rm = TRUE)
data_unique$MSC_mean <- rowMeans(data_unique[, msc_cols], na.rm = TRUE)

#####----- 2. Data Separation (Common / MSC Unique / HEK Unique) -----#####
common_proteins <- data_unique %>% 
  filter(!is.nan(HEK_mean) & !is.nan(MSC_mean))
msc_unique_proteins <- data_unique %>% 
  filter(is.nan(HEK_mean) & !is.nan(MSC_mean) & MSC_mean > 0)
hek_unique_proteins <- data_unique %>% 
  filter(!is.nan(HEK_mean) & is.nan(MSC_mean))


cat("- Common Protein:", nrow(common_proteins), "protein\n")
cat("- MSC Unique Protein:", nrow(msc_unique_proteins), "protein\n")
cat("- HEK Unique Protein Protein:", nrow(hek_unique_proteins), "protein\n")

#####----- 3. MSC Profile analysis -----##### 

msc_top_1000 <- data_unique %>%
  filter(!is.nan(MSC_mean) & MSC_mean > 0) %>%
  arrange(desc(MSC_mean)) %>%
  head(1000)

ids_total <- bitr(msc_top_1000$Gene, fromType="SYMBOL", toType="ENTREZID", OrgDb=org.Hs.eg.db)
kegg_msc_total <- enrichKEGG(gene = ids_total$ENTREZID, organism = 'hsa', pvalueCutoff = 0.5)

  # --- Visualization 1: Dotplot (Blue-Red) ---
  print(">> Figure 1-1: MSC Total Dotplot")
  
  df_total <- as.data.frame(kegg_msc_total)
  
  target_keywords <- c("Endocytosis", "Cytoskeleton", "Adhesion", "Migration", 
                       "Junction", "PI3K", "Rap1", "HIF-1", "Lysosome", "Metabolism")
  
  df_filtered <- df_total %>% 
    filter(grepl(paste(target_keywords, collapse = "|"), Description, ignore.case = TRUE)) %>%
    mutate(
      Pathway_Total = as.numeric(sub("/.*", "", BgRatio)),
      RichFactor = Count / Pathway_Total,
      logP = -log10(p.adjust)
    ) %>%
    arrange(desc(Count))
  
  top_plot <- head(df_filtered, 15) %>% arrange(logP)
  top_plot$Description <- factor(top_plot$Description, levels = top_plot$Description)
  
  p_dot <- ggplot(top_plot, aes(x = logP, y = Description)) +
    geom_point(aes(size = Count, color = RichFactor)) + 
    scale_size_continuous(range = c(5, 8)) +
    scale_color_gradient(low = "blue", high = "red") + 
    theme_bw() +
    labs(title = "Functional Profile of MSC Exosomes (Total)",
         subtitle = "Sorted by Statistical Significance (-log10 p-value)",
         x = "-log10 (Adjusted P-value)", y = "") +
    theme(axis.text.y = element_text(size = 12, face = "bold", color = "black"))
  
  print(p_dot)

  print(">> Figure 1-2: MSC Total Chord Diagram...")
  
  msc_readable <- setReadable(kegg_msc_total, OrgDb = org.Hs.eg.db, keyType = "ENTREZID")
  
  chord_keywords <- c("Endocytosis", "Regulation of actin cytoskeleton", 
                      "Focal adhesion", "Tight junction", "Rap1 signaling", "PI3K-Akt")
  
  chord_df_raw <- as.data.frame(msc_readable) %>%
    dplyr::filter(grepl(paste(chord_keywords, collapse = "|"), Description, ignore.case = TRUE)) %>%
    dplyr::select(Description, geneID) %>%
    tidyr::separate_rows(geneID, sep = "/") %>%
    dplyr::rename(Pathway = Description, Gene = geneID) %>%
    dplyr::select(Gene, Pathway)
  
  gene_counts <- table(chord_df_raw$Gene)
  top_genes <- names(sort(gene_counts, decreasing = TRUE))[1:min(30, length(gene_counts))]
  chord_df <- chord_df_raw %>% dplyr::filter(Gene %in% top_genes)
  pathway_names <- unique(chord_df$Pathway)
  pathway_colors <- setNames(rainbow(length(pathway_names), s = 0.8, v = 0.9), pathway_names)
  gene_names <- unique(chord_df$Gene)
  gene_colors <- setNames(rainbow(length(gene_names), s = 0.5, v = 0.9), gene_names)
  grid_col <- c(pathway_colors, gene_colors)
  link_colors <- pathway_colors[chord_df$Pathway]
  
  order_list <- c(pathway_names, gene_names)
  gap_size <- c(rep(2, length(pathway_names) - 1), 10, rep(2, length(gene_names) - 1), 10)
  
  circos.clear()
  circos.par(start.degree = 90, gap.degree = gap_size, track.margin = c(-0.01, 0.01), points.overflow.warning = FALSE)
  
  chordDiagram(chord_df, grid.col = grid_col, col = link_colors, transparency = 0.4,
               order = order_list, annotationTrack = "grid", preAllocateTracks = 1,
               directional = -1, direction.type = c("arrows"), link.arr.type = "big.arrow", link.arr.length = 0.1)
  
  circos.trackPlotRegion(track.index = 1, panel.fun = function(x, y) {
    xlim = get.cell.meta.data("xlim")
    ylim = get.cell.meta.data("ylim")
    sector.name = get.cell.meta.data("sector.index")
    if(sector.name %in% pathway_names) {
      circos.text(mean(xlim), ylim[1], sector.name, facing = "inside", niceFacing = TRUE, adj = c(0.5, 0), cex = 0.8, font = 2)
    } else {
      circos.text(mean(xlim), 0, sector.name, facing = "clockwise", niceFacing = TRUE, adj = c(0, 0.5), cex = 0.5)
    }
  }, bg.border = NA)


#####----- 4-1. Ligand-Receptor-Pathway axis anylsis -----##### 

common_proteins <- common_proteins %>% mutate(logFC = MSC_mean - HEK_mean)
msc_up_regulated <- common_proteins %>% filter(logFC > 0.26)

target_data <- msc_up_regulated


gene_info <- bitr(target_data$Gene, fromType = "SYMBOL", toType = "GENENAME", OrgDb = org.Hs.eg.db)
gene_info <- gene_info[!duplicated(gene_info$SYMBOL), ]
target_data_annotated <- left_join(target_data, gene_info, by = c("Gene" = "SYMBOL"))
target_data_annotated$GENENAME[is.na(target_data_annotated$GENENAME)] <- target_data_annotated$Gene[is.na(target_data_annotated$GENENAME)]

keywords <- c(
  "BINDING", "LIGAND", "FACTOR", "ADHESION", "MATRIX", 
  "COLLAGEN", "LAMININ", "FIBRONECTIN", "THROMBOSPONDIN", "INTEGRIN",
  "CHEMOKINE", "CYTOKINE",
  "ANNEXIN", "LECTIN", "GALECTIN", "TETRASPANIN", "SCAVENGER", "CLATHRIN", "CAVEOLIN",
  "EPHRIN", "SEMAPHORIN", "RAB", "GTPASE",
  "METALLO", "MMP", "ADAM", "PROTEASE",
  "CADHERIN", "CATENIN", "CD44", "CAM", "SELECTIN"
)

ligand_candidates <- target_data_annotated %>%
  dplyr::filter(grepl(paste(keywords, collapse = "|"), toupper(GENENAME))) 

cat("- Ligand profile:", nrow(ligand_candidates), "protein\n")
  
  ligand_list_view <- ligand_candidates %>%
    dplyr::select(Gene, GENENAME, any_of(c("logFC", "MSC_mean", "HEK_mean"))) %>%
    arrange(desc(MSC_mean)) 

  cat(">> [Top 20] Ligand Candidates (Sorted by Expression):\n")
  print(head(ligand_list_view, 20))

  View(ligand_list_view)
  
  write.csv(ligand_list_view, "MSC_Exosome_Ligand_Candidates_87.csv", row.names = FALSE)
  
  keyword_summary <- data.frame(
    Keyword = c("COLLAGEN|FIBRONECTIN|LAMININ", "INTEGRIN", "VEGF|FGF|PDGF|HGF", "ANNEXIN", "MMP|ADAM", "CD44"),
    Category = c("ECM (Matrix)", "Adhesion (Integrin)", "Growth Factor", "Uptake (Annexin)", "Enzyme (MMP)", "Hyaluronan-R")
  )
  
  for(i in 1:nrow(keyword_summary)) {
    count <- sum(grepl(keyword_summary$Keyword[i], toupper(ligand_list_view$GENENAME)))
    cat(paste0("- ", keyword_summary$Category[i], ": ", count, "개\n"))
  }
  

#####----- 4-2. Ligand-Receptor-Pathway axis anylsis lymphatic vessel targeting -----##### 
lec_receptors_expanded <- data.frame(
  Receptor_Gene = c("FLT4", "KDR", "LYVE1", "PDPN", "PROX1", "FGFR3", "NRP2", "MET", "EGFR", "PDGFRA",
                    "STAB1", "STAB2", "MRC1", "CD209", "CLEC4M", "SCARA3", "CD36",
                    "ICAM1", "VCAM1", "SELE", "SELP", "PECAM1", "CDH5", "MADCAM1",
                    "CCR7", "CXCR4", "ACKR1", "ACKR2", "ACKR3", "ACKR4",
                    "ITGA1", "ITGA4", "ITGA5", "ITGA9", "ITGB1", "ITGB3", "ITGAL", "ITGAM"),
  Category = c(rep("Growth/Marker", 10), rep("Scavenger/Uptake", 7), rep("Cell Adhesion", 7),
               rep("Chemokine/Homing", 6), rep("Integrin Binding", 8))
)

predicted_interactions <- data.frame()

if(nrow(ligand_candidates) > 0) 
  {
  for (i in 1:nrow(ligand_candidates)) {
    gene <- ligand_candidates$Gene[i]
    desc <- toupper(ligand_candidates$GENENAME[i])
    if (grepl("COLLAGEN", desc) || grepl("FIBRONECTIN", desc) || grepl("LAMININ", desc)) {
      targets <- lec_receptors_expanded %>% dplyr::filter(Category == "Integrin Binding")
      if(nrow(targets) > 0) predicted_interactions <- rbind(predicted_interactions, data.frame(Ligand = gene, Receptor = targets$Receptor_Gene))
    }
    if (grepl("CHEMOKINE", desc) || grepl("CXC", desc) || grepl("CCL", desc)) {
      targets <- lec_receptors_expanded %>% dplyr::filter(Category == "Chemokine/Homing")
      if(nrow(targets) > 0) predicted_interactions <- rbind(predicted_interactions, data.frame(Ligand = gene, Receptor = targets$Receptor_Gene))
    }
    if (grepl("VEGF", gene) || grepl("FGF", gene) || grepl("PDGF", gene) || grepl("HGF", gene)) {
      if (grepl("VEGFC", gene)) predicted_interactions <- rbind(predicted_interactions, data.frame(Ligand = gene, Receptor = "FLT4"))
      else if (grepl("VEGFA", gene)) predicted_interactions <- rbind(predicted_interactions, data.frame(Ligand = gene, Receptor = "KDR"))
      else predicted_interactions <- rbind(predicted_interactions, data.frame(Ligand = gene, Receptor = "FGFR3")) 
    }
    if (grepl("INTEGRIN", desc) || grepl("CADHERIN", desc) || grepl("CAM", desc)) {
      targets <- lec_receptors_expanded %>% dplyr::filter(Category == "Cell Adhesion")
      if(nrow(targets) > 0) predicted_interactions <- rbind(predicted_interactions, data.frame(Ligand = gene, Receptor = targets$Receptor_Gene))
    }
    if (grepl("ANNEXIN", desc) || grepl("LECTIN", desc) || grepl("GALECTIN", desc)) {
      targets <- lec_receptors_expanded %>% dplyr::filter(Category == "Scavenger/Uptake")
      if(nrow(targets) > 0) predicted_interactions <- rbind(predicted_interactions, data.frame(Ligand = gene, Receptor = targets$Receptor_Gene))
    }
    if (grepl("SEMAPHORIN", desc)) predicted_interactions <- rbind(predicted_interactions, data.frame(Ligand = gene, Receptor = "NRP2"))
    if (gene == "FN1") predicted_interactions <- rbind(predicted_interactions, data.frame(Ligand = gene, Receptor = "ITGA5"))
    if (gene == "SPP1") predicted_interactions <- rbind(predicted_interactions, data.frame(Ligand = gene, Receptor = "CD44"))
    if (gene == "CD44") predicted_interactions <- rbind(predicted_interactions, data.frame(Ligand = gene, Receptor = "LYVE1"))
  }
}
predicted_interactions <- distinct(predicted_interactions)


#####----- 4-3. Kegg mapping -----##### 
target_receptors <- unique(predicted_interactions$Receptor)
rec_ids <- bitr(target_receptors, fromType="SYMBOL", toType="ENTREZID", OrgDb=org.Hs.eg.db)

kk_rec <- enrichKEGG(gene = rec_ids$ENTREZID, organism = 'hsa', pvalueCutoff = 1, qvalueCutoff = 1)
rec_pathway_clean <- data.frame()

interesting_pathways <- c("Focal adhesion", "ECM-receptor interaction", "Endocytosis", "Phagosome",
                            "PI3K-Akt signaling pathway", "MAPK signaling pathway", "Rap1 signaling pathway",
                            "Chemokine signaling pathway", "Regulation of actin cytoskeleton",
                            "Leukocyte transendothelial migration", "Gap junction", "Tight junction",
                            "VEGF signaling pathway")
  
  rec_pathway_clean <- as.data.frame(kk_rec) %>%
    dplyr::filter(Description %in% interesting_pathways) %>%
    dplyr::select(Description, geneID) %>%
    tidyr::separate_rows(geneID, sep = "/") %>%
    left_join(rec_ids, by = c("geneID" = "ENTREZID")) %>%
    dplyr::rename(Pathway = Description, Receptor = SYMBOL) %>%
    dplyr::select(Receptor, Pathway) %>%
    na.omit()


sankey_3step_up <- predicted_interactions %>%
  inner_join(rec_pathway_clean, by = "Receptor") %>%
  dplyr::select(Ligand, Receptor, Pathway) %>%
  distinct()

#####----- 4-4. Axis Grouping -----##### 

rec_category_map <- data.frame(
  Receptor = c("FLT4", "KDR", "FGFR3", "MET", "NRP2", "EGFR", "PDGFRA",
               "ITGA1", "ITGA4", "ITGA5", "ITGA9", "ITGB1", "ITGB3", "ITGAL",
               "CD44", "LYVE1", 
               "ICAM1", "VCAM1", "SELE", "SELP", "CDH5",
               "CXCR4", "CCR7", "ACKR1", 
               "STAB1", "STAB2", "MRC1", "CD36"),
  Category = c(rep("Growth Factor R", 7),
               rep("Integrins", 7),
               rep("Hyaluronan Axis", 2),
               rep("Adhesion Molecules", 5),
               rep("Chemokine Receptors", 3),
               rep("Scavenger Receptors", 4))
)

sankey_mapped <- sankey_3step_up %>%
  left_join(rec_category_map, by = "Receptor") %>%
  mutate(Receptor_Display = ifelse(is.na(Category), "Other Receptors", Category))


sankey_grouped <- sankey_mapped %>%
  mutate(Ligand_Group = case_when(
    grepl("^COL", Ligand) ~ "Collagens",
    grepl("^LAM", Ligand) ~ "Laminins",
    grepl("^ITG", Ligand) ~ "Integrins (Exo)",
    grepl("FN1", Ligand) ~ "Fibronectin",
    grepl("THBS", Ligand) ~ "Thrombospondins",
    grepl("VEGF|FGF|PDGF|HGF|EGF", Ligand) ~ "Growth Factors",
    grepl("MMP|ADAM", Ligand) ~ "MMPs/ADAMs",
    grepl("ANXA", Ligand) ~ "Annexins",
    grepl("CD44", Ligand) ~ "Hyaluronan-R", 
    TRUE ~ "Other Ligands" 
  )) %>%
  filter(Ligand_Group != "Other Ligands") 


plot_data_balanced <- sankey_grouped %>%
  group_by(Ligand_Group) %>%
  distinct(Receptor_Display, Pathway, .keep_all = TRUE) %>%
  slice_head(n = 5) %>% 
  ungroup() %>%
  mutate(freq = 1) 

plot_data_balanced$Pathway <- stringr::str_wrap(plot_data_balanced$Pathway, width = 20)



#####----- 4-5. Visulalization  -----##### 
unique_groups <- sort(unique(plot_data_balanced$Ligand_Group))
n_groups <- length(unique_groups)

mycol.vector <- c('#9e0142','#d53e4f','#f46d43','#fdae61','#fee08b',
                  '#e6f598','#abdda4','#66c2a5','#3288bd','#5e4fa2')
final_colors <- rep(mycol.vector, times = ceiling(n_groups / length(mycol.vector)))[1:n_groups]

p_sankey_final <- ggplot(plot_data_balanced,
                         aes(y = freq, axis1 = Ligand_Group, axis2 = Receptor_Display, axis3 = Pathway)) +
  geom_alluvium(aes(fill = Ligand_Group), width = 1/20, alpha = 0.7) +
  geom_stratum(width = 1/20, fill = "white", color = "grey30", size = 0.3) +
  geom_label(stat = "stratum", aes(label = after_stat(stratum)), 
             size = 3, label.size = 0, fill = NA, fontface = "bold") +
  scale_fill_manual(values = final_colors) +
  scale_x_discrete(limits = c("Ligand Family\n(MSC Upregulated)", "Receptor Group\n(Lymphatic)", "Target Pathway\n(Intracellular)"), 
                   expand = c(.05, .05)) +
  labs(title = "Mechanism of MSC Exosome Uptake (Balanced View)",
       subtitle = "Key Ligand Families -> Functional Receptors -> Signaling Pathways",
       y = "Representative Interactions (Top 3 per Group)") +
  theme_void() +
  theme(
    plot.title = element_text(hjust = 0.5, face = "bold", size = 16),
    plot.subtitle = element_text(hjust = 0.5, size = 12),
    axis.text.x = element_text(size = 11, face = "bold", vjust = 1, color = "black"), 
    legend.position = "none"
  )

print(p_sankey_final)

library(dplyr)


#####----- 4-6. Final result Ligand list  -----##### 
  ligand_summary <- sankey_grouped %>%
    dplyr::select(Ligand, Ligand_Group, Receptor_Display) %>%
    distinct() %>%
    arrange(Ligand_Group, Ligand) %>% 
    group_by(Ligand, Ligand_Group) %>%
    summarise(Target_Receptors = paste(unique(Receptor_Display), collapse = ", "), .groups = "drop")
  

  cat("\n>> [Total Count] Found ligand number:", nrow(ligand_summary), "number\n")
  print(as.data.frame(ligand_summary), row.names = FALSE)
  print(table(ligand_summary$Ligand_Group))







