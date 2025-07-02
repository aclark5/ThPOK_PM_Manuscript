# Thpok mutaion scRNA data analysis
# Author: Xinxin chen
# Data: 21/05/2025
# Update Date: 02/07/2025


library(Seurat) # version: 5.0.1
library("ggplot2")
library(reshape2)
library(ggpubr)
library(scales)


# separately analysis Neg and stimulated data, but the pipeline is similar
samples <- list(
  list("path" = "/Singlecellrna/8_samples/cellranger_count/HC4-ThN-D7-Neg-1_out/outs/filtered_feature_bc_matrix", 
       "project" = "HC4_Neg"),
  list("path" = "/Singlecellrna/8_samples/cellranger_count/HC5-ThN-D7-Neg-1_out/outs/filtered_feature_bc_matrix", 
       "project" = "HC5_Neg"),
  list("path" = "/Singlecellrna/8_samples/cellranger_count/LCKE-ThN-D7-Neg-1_out/outs/filtered_feature_bc_matrix", 
       "project" = "LCKE_Neg"),
  list("path" = "/Singlecellrna/8_samples/cellranger_count/LKKS-ThN-D7-Neg-1_out/outs/filtered_feature_bc_matrix", 
       "project" = "LKKS_Neg")
)
# separately analysis Neg and stimulated data, but the pipeline is similar
samples <- list(
  list("path" = "/Singlecellrna/8_samples/cellranger_count/HC4-ThN-D7-iono-1_out/outs/filtered_feature_bc_matrix", "project" = "HC4_iono"),
  list("path" = "/Singlecellrna/8_samples/cellranger_count/HC5-ThN-D7-iono-1_out/outs/filtered_feature_bc_matrix", "project" = "HC5_iono"),
  list("path" = "/Singlecellrna/8_samples/cellranger_count/LCKE-ThN-D7-iono-1_out/outs/filtered_feature_bc_matrix", "project" = "LCKE_iono"),
  list("path" = "/Singlecellrna/8_samples/cellranger_count/LKKS-ThN-D7-iono-1_out/outs/filtered_feature_bc_matrix", "project" = "LKKS_iono")
)


# Initialize a list to store Seurat objects
seurat_objects <- list()

# Process each sample
for (sample in samples) {
  input_matrix <- sample$path
  print(input_matrix)
  project_name <- sample$project
  print(project_name)
  
  # Read and process the data
  seurat_object <- Seurat::Read10X(input_matrix)
  seurat_object <- Seurat::CreateSeuratObject(counts = seurat_object, project = project_name, min.cells = 3, min.features = 200)
  seurat_object[["percent.mt"]] <- Seurat::PercentageFeatureSet(seurat_object, pattern = "^MT-")
  seurat_object[["percent.rb"]] <- Seurat::PercentageFeatureSet(seurat_object, pattern = "^RP[SL]")
  seurat_object[["percent.hb"]] <- Seurat::PercentageFeatureSet(seurat_object, pattern = "^HB[A-Z0-9]*")
  
  # Run external Python script for doublet detection (assumed to be run manually or separately)
  system(paste("python /Singlecellrna/8_samples/cellranger_count/my_scrublet_run.py", input_matrix))
  
  # Load doublet results
  doublets <- read.table(file.path(input_matrix, "doublet.txt"), header = TRUE, sep = ",", row.names = 1)
  colnames(doublets) <- c("Doublet_score", "Is_doublet")
  seurat_object <- Seurat::AddMetaData(seurat_object, doublets)
  
  # Print head of metadata
  print(head(seurat_object[[]]))
  
  # Filter the Seurat object 200, 9000, 20% mt
  seurat_object <- subset(seurat_object, subset = nCount_RNA > 1000 & nFeature_RNA > 200 & nFeature_RNA < 7500 & percent.mt < 10 & Doublet_score < 0.22 & percent.hb < 1)
  
  
  # Store the object in the list with a named key
  seurat_objects[[project_name]] <- seurat_object
}

# scale
# Normalize and find variable features for each Seurat object
for (name in names(seurat_objects)) {
  seurat_objects[[name]] <- NormalizeData(seurat_objects[[name]], verbose = FALSE)
  seurat_objects[[name]] <- FindVariableFeatures(seurat_objects[[name]], selection.method = "vst", nfeatures = 2000)
}

T_anchors <- Seurat::FindIntegrationAnchors(object.list = seurat_objects, dims = 1:30)
seurat_integrated  <- Seurat::IntegrateData(anchorset = T_anchors, dims = 1:30)

saveRDS(seurat_integrated, file = "/Singlecellrna/8_samples/seurat_out/T_seurat_integrated.rds")


# read data
seurat_integrated <- readRDS("/Singlecellrna/8_samples/seurat_out/T_seurat_integrated.rds")
all.genes <- rownames(seurat_integrated)
seurat_integrated <- Seurat::ScaleData(seurat_integrated, features = all.genes, verbose = F)
seurat_integrated <- Seurat::RunPCA(seurat_integrated, npcs = 30, features = Seurat::VariableFeatures(object = seurat_integrated), verbose = F)


saveRDS(seurat_integrated, file = "/Singlecellrna/8_samples/seurat_out/T_seurat_integrated_pca30.rds")


seurat_integrated <- readRDS("/Singlecellrna/8_samples/seurat_out/T_seurat_integrated_pca30.rds")
seurat_integrated <- Seurat::FindNeighbors(seurat_integrated, dims = 1:20)
# seurat_integrated <- Seurat::FindClusters(seurat_integrated, resolution = seq(0.1, 1.5, by = 0.1)) find the best resolution
seurat_integrated <- Seurat::FindClusters(seurat_integrated, resolution = 1)
seurat_integrated <- Seurat::RunUMAP(seurat_integrated, reduction = "pca", dims = 1:20, verbose = F)
saveRDS(seurat_integrated, file = "/Singlecellrna/8_samples/seurat_out/T_seurat_umap_d20_mt10_4_samples.rds")


# annotation, Neg and iono are separately analysis
T_seurat_4_samples <- readRDS("./T_seurat_umap_d20_mt10_4_samples.rds")

# feature plot
DefaultAssay(T_seurat_4_samples) <- "RNA"
marker <- c("CD3E","CD3D","CD4","CD8A","NCAM1","FCGR3A")
p <- Seurat::FeaturePlot(T_seurat_4_samples, cols = c("lightgrey", "red"),
                         features = marker, combine = TRUE, raster = FALSE)  #max.cutoff = 3
p|Seurat::DimPlot(T_seurat_4_samples, label = TRUE)


# VlnPlot visulization
library(paletteer)
VlnPlot(T_seurat_4_samples, features = marker, cols = c(paletteer_d("awtools::bpalette"),
                                                                 paletteer_d("awtools::a_palette"),
                                                                 paletteer_d("awtools::mpalette"),
                                                                 paletteer_d("awtools::spalette")),combine = T, raster=F)

Seurat::FeaturePlot(T_seurat_4_samples, features = c("CD4","CD8A"), split.by = "orig.ident")

# Neg groups
cells_to_remove <- Seurat::WhichCells(T_seurat_4_samples, idents = c(7,19,21,23,24))
T_seurat_4_samples <- subset(T_seurat_4_samples, cells = setdiff(Cells(T_seurat_4_samples), cells_to_remove))
# iono groups
cells_to_remove <- Seurat::WhichCells(T_seurat_4_samples, idents = c(14,17,20,21)) # 17,19,20,21,22,23
T_seurat_4_samples <- subset(T_seurat_4_samples, cells = setdiff(Cells(T_seurat_4_samples), cells_to_remove))

marker <- c("CD3E","CD4","CD8A")
p <- Seurat::FeaturePlot(T_seurat_4_samples, cols = c("lightgrey", "red"),
                         features = marker, combine = TRUE, raster = FALSE, ncol = 3)
ggsave(p,
       filename = "./Tcell_marker_plot_4samples_Neg.pdf",
       width = 10, height = 2.3)

# Neg, CD4+ & CD8+ plot
p|Seurat::DimPlot(T_seurat_4_samples, label = TRUE)
old_levels <- levels(T_seurat_4_samples)
cd4_clusters <- c(0,3,4,9,12,14,16,18)  # CD4+ clusters
new.cluster.ids <- setNames(rep("CD4+", length(old_levels)), old_levels)
cd8_clusters <- c(1,2,5,6,7,8,10,11,13,15)  # CD8+ clusters
cd8_clusters <- as.character(cd8_clusters)
new.cluster.ids[cd8_clusters] <- "CD8+"
T_seurat_4_samples_rename <-  Seurat::RenameIdents(T_seurat_4_samples, new.cluster.ids)
p <- Seurat::DimPlot(T_seurat_4_samples_rename, label = F, split.by = "orig.ident", cols = c("#FFA500","#4B0082"))
ggsave(p,
       filename = "./Tcell_two_groups_plot_4samples_Neg.pdf",
       width = 6.3, height = 3)

# T cell proportion, Neg sample
bar_plot <- plot_celltype_proportions(T_seurat_4_samples_rename)
ggsave(bar_plot,
       filename = "./bar_plot_Tcell_4samples_Neg.pdf",
       width = 3.2, height = 4.8)

# Neg, CD4+ 
T_seurat_4_samples_cd4 <- subset(T_seurat_4_samples, idents = c(0,3,4,9,12,14,16,18))
Seurat::DimPlot(T_seurat_4_samples_cd4, label = TRUE)
# iono, CD4+
T_seurat_4_samples_cd4 <- subset(T_seurat_4_samples, idents = c(0,3,5,9,10,15,16,19))
Seurat::DimPlot(T_seurat_4_samples_cd4, label = TRUE)

# Neg, CD8+
T_seurat_4_samples_cd8 <- subset(T_seurat_4_samples, idents = c(1,2,5,6,7,8,10,11,13,15))
Seurat::DimPlot(T_seurat_4_samples_cd8, label = TRUE)
# iono, CD8+
T_seurat_4_samples_cd8 <- subset(T_seurat_4_samples, idents = c(1,2,4,6,7,8,11,12,13,18))
Seurat::DimPlot(T_seurat_4_samples_cd8, label = TRUE)

# Neg, recluster
DefaultAssay(T_seurat_4_samples_cd4) <- "integrated"
T_seurat_4_samples_cd4_re <- RunPCA(T_seurat_4_samples_cd4, npc = 20)
T_seurat_4_samples_cd4_re <- FindNeighbors(T_seurat_4_samples_cd4_re, dim = 1:20)
T_seurat_4_samples_cd4_re <- FindClusters(T_seurat_4_samples_cd4_re,resolution = 0.8)  # change 1 to 1.5
T_seurat_4_samples_cd4_re <- RunUMAP(T_seurat_4_samples_cd4_re, reduction = "pca", dim = 1:20)
Seurat::DimPlot(T_seurat_4_samples_cd4_re, label = TRUE)
T_seurat_4_samples_cd4_re <- subset(T_seurat_4_samples_cd4_re, idents = c(0,1,2,3,4,5,6,7,8,9,12,13,14,15))
Seurat::DimPlot(T_seurat_4_samples_cd4_re, label = TRUE)

# Neg, marker gene heatmap --------------
Seurat::FeaturePlot(T_seurat_4_samples_cd4_re, features = c("CD4","CD8A","CCR7","CXCR6","FOXP3","TYMS"), cols = c("lightgrey", "red"),
                    ncol = 3, max.cutoff = 2)

marker <- c("CCR7","IL7R","TCF7","SELL","CXCR5","TOX2","SLC9A9","ZC3H12D","PRDM1","CXCR6","TBX21","RUNX3","ID2","IFNG","PRF1","CCL4",
            "GZMA","GZMB","GZMK","FOXP3","IKZF2","IL1R1","TIGIT","TYMS","MKI67","CLSPN")
scRNAtoolVis::averageHeatmap(object = T_seurat_4_samples_cd4_re,htCol = c("steelblue3","white","red3"),border=F,
                             markerGene = marker, fontsize = 10, annoCol = T,
                             gene.order = c("CCR7","IL7R","TCF7","SELL","CXCR5","TOX2","SLC9A9","ZC3H12D","PRDM1","CXCR6","TBX21","RUNX3","ID2","IFNG","PRF1","CCL4",
                                            "GZMA","GZMB","GZMK","FOXP3","IKZF2","IL1R1","TIGIT","TYMS","MKI67","CLSPN"),width = 3, height = 8)

scRNAtoolVis::averageHeatmap(object = T_seurat_4_samples_cd4_re_rename,htCol = c("steelblue3","white","red3"),border=F,
                             markerGene = marker, fontsize = 10, annoCol = T, 
                             cluster.order = c("Tcm","Tfh","Th1/cytotoxic","Treg","Cycling"),
                             gene.order = c("CCR7","IL7R","TCF7","SELL","CXCR5","TOX2","SLC9A9","ZC3H12D","PRDM1","CXCR6","TBX21","RUNX3","ID2","IFNG","PRF1","CCL4",
                                            "GZMA","GZMB","GZMK","FOXP3","IKZF2","IL1R1","TIGIT","TYMS","MKI67","CLSPN"),width = 3, height = 8,
                             myanCol = c("#FFA500", "#FF0000","#00FF00","#0000FF", "gray"))


# DEG (clusters)
T_seurat_4_samples_cd4_re <- SeuratObject::JoinLayers(T_seurat_4_samples_cd4_re)
T_seurat_4_samples_cd4_re_findallmaker_df <- FindAllMarkers(T_seurat_4_samples_cd4_re, only.pos = T, min.pct = 0.25, logfc.threshold = 0.5)
write.csv(T_seurat_4_samples_cd4_re_findallmaker_df, file = "./12_samples_20240722/only_8_samples/T_seurat_4_samples_neg_cd4_findallmakers_df.csv", row.names = FALSE)

# Neg rename
T_seurat_4_samples_cd4_re@meta.data$orig.ident[T_seurat_4_samples_cd4_re@meta.data$orig.ident == "HC4_Neg"] <- "Controls"
T_seurat_4_samples_cd4_re@meta.data$orig.ident[T_seurat_4_samples_cd4_re@meta.data$orig.ident == "HC5_Neg"] <- "Controls"
T_seurat_4_samples_cd4_re@meta.data$orig.ident[T_seurat_4_samples_cd4_re@meta.data$orig.ident == "LCKE_Neg"] <- "Patients"
T_seurat_4_samples_cd4_re@meta.data$orig.ident[T_seurat_4_samples_cd4_re@meta.data$orig.ident == "LKKS_Neg"] <- "Patients"

new.cluster.ids <- c("Tcm","Tcm","Tcm","Tcm","Cycling","Th1/cytotoxic","Tcm","Tfh","Cycling","Tcm",
                     "Tcm","Treg","Tcm","Tcm")
names(new.cluster.ids) <- levels(T_seurat_4_samples_cd4_re)
T_seurat_4_samples_cd4_re_rename <- Seurat::RenameIdents(T_seurat_4_samples_cd4_re, new.cluster.ids)
p <- Seurat::DimPlot(T_seurat_4_samples_cd4_re_rename, label = F, split.by = "orig.ident", cols = c("#FFA500","gray","#00FF00","#FF0000","#0000FF" ,"#FF00FF"))
ggsave(p,
       filename = "./CD4_Tcell_two_groups_plot_4samples_Neg.pdf",
       width = 6.3, height = 3)


# iono rename
new.cluster.ids <- c("Tcm-stimulated","Tcm-stimulated","Tcm-stimulated","Tcm-stimulated","Cycling-stimulated","Tcm-stimulated","Th1/cytotoxic-stimulated","Cycling-stimulated","Tfh-stimulated","Cycling-stimulated",
                     "Cycling-stimulated","Treg-stimulated","Tcm-stimulated")
names(new.cluster.ids) <- levels(T_seurat_4_samples_cd4_re)
T_seurat_4_samples_cd4_re_rename <- Seurat::RenameIdents(T_seurat_4_samples_cd4_re, new.cluster.ids)
p <- Seurat::DimPlot(T_seurat_4_samples_cd4_re_rename, label = F, split.by = "orig.ident", cols = c("#BDD7E7","#00FF00","#FF0000","#0000FF" ,"#FF00FF"))
ggsave(p,
       filename = "./CD4_Tcell_two_groups_plot_4samples_iono.pdf",
       width = 6.3, height = 3)

# DEG (Patients vs controls)
T_seurat_4_samples_cd4_re <- JoinLayers(T_seurat_4_samples_cd4_re, assay = "RNA", layer = "data")
T_seurat_4_samples_cd4_re_markers <- Seurat::FindMarkers(T_seurat_4_samples_cd4_re, assay = "RNA", layer = "data",ident.1 = "Patients", 
                                                         group.by = 'orig.ident', min.pct = 0.25, logfc.threshold = 0.5)
T_seurat_4_samples_cd4_re_markers_df <- T_seurat_4_samples_cd4_re_markers %>% as.data.frame(.) %>% rownames_to_column(., "genes")
write.csv(T_seurat_4_samples_cd4_re_markers_df, file = "./12_samples_20240722/only_8_samples/T_seurat_4_samples_cd4_re_rename_DEG_neg_df.csv", row.names = FALSE)

# DEG (clusters)
T_seurat_4_samples_cd4_re_rename <- JoinLayers(T_seurat_4_samples_cd4_re_rename, assay = "RNA", layer = "data")
T_seurat_4_samples_cd4_re_rename_findallmaker_df <- FindAllMarkers(T_seurat_4_samples_cd4_re_rename, only.pos = T, min.pct = 0.25, logfc.threshold = 0.5)
write.csv(T_seurat_4_samples_cd4_re_rename_findallmaker_df, file = "./12_samples_20240722/only_8_samples/T_seurat_4_samples_cd4_re_rename_findallmaker_df.csv", row.names = FALSE)

# DEG (Tfh and Tcm)
T_seurat_4_samples_cd4_re_rename_Tfh_marker <- Seurat::FindMarkers(T_seurat_4_samples_cd4_re_rename, assay = "RNA", layer = "data",ident.1 = "Patients",
                                                                   group.by = 'orig.ident', min.pct = 0.01, logfc.threshold = 0.1, subset.ident = "Tfh")
T_seurat_4_samples_cd4_re_rename_Tcm_marker <- Seurat::FindMarkers(T_seurat_4_samples_cd4_re_rename, assay = "RNA", layer = "data",ident.1 = "Patients",
                                                                   group.by = 'orig.ident', min.pct = 0.01, logfc.threshold = 0.1, subset.ident = "Tcm")
T_seurat_4_samples_cd4_re_rename_Tfh_marker_df <- T_seurat_4_samples_cd4_re_rename_Tfh_marker %>% as.data.frame(.) %>% rownames_to_column(., "genes")
T_seurat_4_samples_cd4_re_rename_Tcm_marker_df <- T_seurat_4_samples_cd4_re_rename_Tcm_marker %>% as.data.frame(.) %>% rownames_to_column(., "genes")
write.csv(T_seurat_4_samples_cd4_re_rename_Tfh_marker_df, file = "./12_samples_20240722/only_8_samples/T_seurat_4_samples_cd4_re_rename_Tfh_DEG_neg_df.csv", row.names = FALSE)
write.csv(T_seurat_4_samples_cd4_re_rename_Tcm_marker_df, file = "./12_samples_20240722/only_8_samples/T_seurat_4_samples_cd4_re_rename_Tcm_DEG_neg_df.csv", row.names = FALSE)

# Feature plot
DefaultAssay(T_seurat_4_samples_cd4_re) <- "RNA"
Seurat::FeaturePlot(T_seurat_4_samples_cd4_re, features = c("CD4","CD8A"), cols = c("lightgrey", "red"), split.by = "orig.ident")
features_plot <- Seurat::FeaturePlot(T_seurat_4_samples_cd4_re_rename, features = c("CD4","CCR7","CXCR5","FOXP3","IFNG","MKI67"), cols = c("lightgrey", "red"), ncol = 3)
ggsave(features_plot, file = "./features_plot_4samples_neg.pdf",
       width = 8.5, height = 5)

# violin_plot
library(dplyr)
library(tidyr)
gene_of_interest <- c("ZBTB7B")
gene_of_interest <- c("RUNX3")
gene_of_interest <- c("IFNG","NKG7","SLAMF7","CRTAM","PRF1","KLRC1","CCL4","CCL5","GZMA","GZMB","GZMH","GZMM")
df <- FetchData(T_seurat_4_samples_cd4_re, vars = c(gene_of_interest, "orig.ident")) %>%
  pivot_longer(cols = all_of(gene_of_interest), 
               names_to = "gene", 
               values_to = "gene_expression") %>%
  rename(Group = orig.ident) %>%
  mutate(Group = factor(Group, levels = c("Controls", "Patients")))
violin_plot_thpok <- ggplot(df, aes(x = Group, y = gene_expression, fill = Group)) +
  geom_jitter(aes(color = Group), width = 0.15, alpha = 0.2, size = 0.5) +  # Jitter for individual points
  geom_violin(alpha = 0.6, trim = FALSE, width = 0.6) +  # Violin plot for distribution
  geom_boxplot(width = 0.15, outlier.shape = NA, alpha = 0.3) +  # Boxplot inside violin
  facet_wrap(~gene, scales = "free_y") +  # Facet by gene, allowing independent y scales
  theme_classic() +
  scale_fill_manual(values = c("Controls" = "#1f77b4", "Patients" = "#ff7f0e")) + 
  scale_color_manual(values = c("Controls" = "#1f77b4", "Patients" = "#ff7f0e")) +
  labs(title = "Gene Expression Across Cell Types",
       x = "Group",
       y = "Gene Expression") +
  theme(legend.position = "top",
        plot.title = element_text(size = 12),        # Title size
        axis.title.x = element_text(size = 12),                     # X-axis title size
        axis.title.y = element_text(size = 12),                     # Y-axis title size
        axis.text = element_text(size = 10),                        # Axis tick labels size
        strip.text = element_text(size = 10),        # Facet strip text size
        legend.text = element_text(size = 10),                      # Legend text size
        legend.title = element_text(size = 12))

ggsave(violin_plot_thpok, file = "./violin_plot_cytotoxic_4samples_neg.pdf",
       width = 10, height = 5.5)
ggsave(violin_plot_thpok, file = "./violin_plot_cellmarkers_4samples_neg.pdf",
       width = 4, height = 3.5)

# volcano plot for Tfh
T_seurat_4_samples_cd4_re_rename_Tfh_marker_df$diffexpressed <- "NO"
T_seurat_4_samples_cd4_re_rename_Tfh_marker_df$diffexpressed[T_seurat_4_samples_cd4_re_rename_Tfh_marker_df$avg_log2FC > 0.5 & T_seurat_4_samples_cd4_re_rename_Tfh_marker_df$p_val_adj < 0.05] <- "UP"
T_seurat_4_samples_cd4_re_rename_Tfh_marker_df$diffexpressed[T_seurat_4_samples_cd4_re_rename_Tfh_marker_df$avg_log2FC < -0.5 & T_seurat_4_samples_cd4_re_rename_Tfh_marker_df$p_val_adj < 0.05] <- "DOWN"

mycolors <- c("grey", "red", "black")
names(mycolors) <- c("DOWN", "UP", "NO")
T_seurat_4_samples_cd4_re_rename_Tfh_marker_df$label <- NA

T_seurat_4_samples_cd4_re_rename_Tfh_marker_df$label[T_seurat_4_samples_cd4_re_rename_Tfh_marker_df$genes %in% go_genes_up$gene] <- 
  T_seurat_4_samples_cd4_re_rename_Tfh_marker_df$genes[T_seurat_4_samples_cd4_re_rename_Tfh_marker_df$genes %in%go_genes_up$gene] 


# library(ggrepel)
volcano_DP <- ggplot(data=T_seurat_4_samples_cd4_re_rename_Tfh_marker_df, 
                     aes(x=avg_log2FC, 
                         y=-log10(p_val), 
                         col=diffexpressed, 
                         label=label)) +
  geom_point(size = 0.5) + 
  theme_classic() +
  geom_text_repel(size = 4,
                  min.segment.length = 0) +
  scale_color_manual(values=mycolors) +
  xlab("Average Log2 Fold Change") +
  ylab("-log10 p-value") +
  xlim(c(-2.8, 2.8)) +
  theme(axis.text=element_text(size=12),
        axis.title=element_text(size=12)) +
  NoLegend()

ggsave(volcano_DP,
       filename = "./volcano_Tfh.pdf",
       width = 4, height = 4)

# volcano plot for Tcm
T_seurat_4_samples_cd4_re_rename_Tcm_marker_df$diffexpressed <- "NO"
T_seurat_4_samples_cd4_re_rename_Tcm_marker_df$diffexpressed[T_seurat_4_samples_cd4_re_rename_Tcm_marker_df$avg_log2FC > 0.5 & T_seurat_4_samples_cd4_re_rename_Tcm_marker_df$p_val_adj < 0.05] <- "UP"
T_seurat_4_samples_cd4_re_rename_Tcm_marker_df$diffexpressed[T_seurat_4_samples_cd4_re_rename_Tcm_marker_df$avg_log2FC < -0.5 & T_seurat_4_samples_cd4_re_rename_Tcm_marker_df$p_val_adj < 0.05] <- "DOWN"

mycolors <- c("grey", "red", "black")
names(mycolors) <- c("DOWN", "UP", "NO")
T_seurat_4_samples_cd4_re_rename_Tcm_marker_df$label <- NA

T_seurat_4_samples_cd4_re_rename_Tcm_marker_df$label[T_seurat_4_samples_cd4_re_rename_Tcm_marker_df$genes %in% go_genes_up$gene] <- 
  T_seurat_4_samples_cd4_re_rename_Tcm_marker_df$genes[T_seurat_4_samples_cd4_re_rename_Tcm_marker_df$genes %in%go_genes_up$gene] 


# library(ggrepel)
volcano_DP <- ggplot(data=T_seurat_4_samples_cd4_re_rename_Tcm_marker_df, 
                     aes(x=avg_log2FC, 
                         y=-log10(p_val), 
                         col=diffexpressed, 
                         label=label)) +
  geom_point(size = 0.5) + 
  theme_classic() +
  geom_text_repel(size = 4,
                  min.segment.length = 0) +
  scale_color_manual(values=mycolors) +
  xlab("Average Log2 Fold Change") +
  ylab("-log10 p-value") +
  xlim(c(-2.8, 2.8)) +
  theme(axis.text=element_text(size=12),
        axis.title=element_text(size=12)) +
  NoLegend()

ggsave(volcano_DP,
       filename = "./volcano_Tcm.pdf",
       width = 4, height = 4)


# CD4+ proportion, Neg
T_seurat_4_samples_cd4_re_rename <- subset(T_seurat_4_samples_cd4_re_rename, idents = c("Tcm","Tfh","Th1/cytotoxic","Treg"))
bar_plot <- plot_celltype_proportions(T_seurat_4_samples_cd4_re_rename)
ggsave(bar_plot,
       filename = "./bar_plot_CD4_Tcell_4samples_Neg.pdf",
       width = 3.2, height = 4.8)

# proportion function
library(Seurat)
library(ggplot2)
library(dplyr)
library(ggpubr)
library(scales)
plot_celltype_proportions <- function(seurat_obj) {
  seurat_obj@meta.data$orig.ident[seurat_obj@meta.data$orig.ident == "HC4_Neg"]    <- "Controls"
  seurat_obj@meta.data$orig.ident[seurat_obj@meta.data$orig.ident == "HC5_Neg"]    <- "Controls"
  seurat_obj@meta.data$orig.ident[seurat_obj@meta.data$orig.ident == "LCKE_Neg"]   <- "Patients"
  seurat_obj@meta.data$orig.ident[seurat_obj@meta.data$orig.ident == "LKKS_Neg"]   <- "Patients"
  metadata <- seurat_obj@meta.data
  cell_type_sample_table <- table(Idents(seurat_obj), metadata$orig.ident)
  cell_type_sample_df <- as.data.frame(cell_type_sample_table)
  colnames(cell_type_sample_df) <- c("Cell_Type", "Sample", "Count")
  # Calculate total cells per sample
  total_cells_per_sample <- aggregate(Count ~ Sample, data = cell_type_sample_df, sum)
  # Merge to compute proportions
  cell_type_sample_df <- merge(cell_type_sample_df, total_cells_per_sample, by = "Sample", suffixes = c("", "_Total"))
  cell_type_sample_df$Proportion <- cell_type_sample_df$Count / cell_type_sample_df$Count_Total
  # Define custom colors
  custom_colors <- c("#FFA500", "#00FF00", "#FF0000", "#0000FF", "#8c564b", "#8a8cd8", "gray")
  # Plot
  p <- ggplot(cell_type_sample_df, aes(x = Sample, y = Proportion, fill = Cell_Type)) +
    geom_bar(stat = "identity", position = "fill") +
    theme_bw() +
    labs(
      title = "RNA Cluster Proportions",
      x = NULL,
      y = NULL,
      fill = "RNA Cluster",
      size = 2
    ) +
    theme(axis.text.x = element_text(angle = 45, hjust = 1, size = 10)) +
    scale_y_continuous(labels = scales::percent) +
    scale_fill_manual(values = custom_colors)
  return(p)
}






