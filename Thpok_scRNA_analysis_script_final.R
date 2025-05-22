# Thpok mutaion scRNA data analysis
# Author: Xinxin chen 
# Date: 21/05/2025


library(Seurat) # version: 5.0.1
library("ggplot2")
library(reshape2)
library(ggpubr)
library(scales)



samples <- list(
  list("path" = "/Singlecellrna/8_samples/cellranger_count/HC4-ThN-D7-iono-1_out/outs/filtered_feature_bc_matrix", "project" = "HC4_iono"),
  list("path" = "/Singlecellrna/8_samples/cellranger_count/HC5-ThN-D7-iono-1_out/outs/filtered_feature_bc_matrix", "project" = "HC5_iono"),
  list("path" = "/Singlecellrna/8_samples/cellranger_count/LCKE-ThN-D7-iono-1_out/outs/filtered_feature_bc_matrix", "project" = "LCKE_iono"),
  list("path" = "/Singlecellrna/8_samples/cellranger_count/LKKS-ThN-D7-iono-1_out/outs/filtered_feature_bc_matrix", "project" = "LKKS_iono"),
  list("path" = "/Singlecellrna/8_samples/cellranger_count/HC4-ThN-D7-Neg-1_out/outs/filtered_feature_bc_matrix", 
       "project" = "HC4_Neg"),
  list("path" = "/Singlecellrna/8_samples/cellranger_count/HC5-ThN-D7-Neg-1_out/outs/filtered_feature_bc_matrix", 
       "project" = "HC5_Neg"),
  list("path" = "/Singlecellrna/8_samples/cellranger_count/LCKE-ThN-D7-Neg-1_out/outs/filtered_feature_bc_matrix", 
       "project" = "LCKE_Neg"),
  list("path" = "/Singlecellrna/8_samples/cellranger_count/LKKS-ThN-D7-Neg-1_out/outs/filtered_feature_bc_matrix", 
       "project" = "LKKS_Neg")
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
  seurat_object <- subset(seurat_object, subset = nCount_RNA > 1000 & nFeature_RNA > 200 & nFeature_RNA < 7500 & percent.mt < 15 & Doublet_score < 0.22 & percent.hb < 1)
  
  
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
seurat_integrated <- Seurat::FindClusters(seurat_integrated, resolution = 1)
seurat_integrated <- Seurat::RunUMAP(seurat_integrated, reduction = "pca", dims = 1:20, verbose = F)
saveRDS(seurat_integrated, file = "/Singlecellrna/8_samples/seurat_out/T_seurat_umap_d20_mt15_8_samples.rds")


# annotation
T_seurat_d20_mt15_8_samples <- readRDS("./T_seurat_umap_d20_mt15_8_samples.rds")

# feature plot
DefaultAssay(T_seurat_d20_mt15_8_samples) <- "RNA"
marker <- c("CD3E","CD3D","CD4","CD8A","NCAM1","FCGR3A")
p <- Seurat::FeaturePlot(T_seurat_d20_mt15_8_samples, cols = c("lightgrey", "red"),
                         features = marker, combine = TRUE, raster = FALSE)  #max.cutoff = 3
p|Seurat::DimPlot(T_seurat_d20_mt15_8_samples, label = TRUE)


# VlnPlot visulization
library(paletteer)
VlnPlot(T_seurat_d20_mt15_8_samples, features = marker, cols = c(paletteer_d("awtools::bpalette"),
                                                                 paletteer_d("awtools::a_palette"),
                                                                 paletteer_d("awtools::mpalette"),
                                                                 paletteer_d("awtools::spalette")),combine = T, raster=F)

Seurat::FeaturePlot(T_seurat_d20_mt15_8_samples, features = c("CD4","CD8A"), split.by = "orig.ident")

# findsubcluster
T_seurat2 <- Seurat::FindSubCluster(T_seurat_d20_mt15_8_samples, cluster = "3", graph.name = "integrated_snn", 
                                    subcluster.name = "integrated_snn_res.1_c1_sub", resolution = 0.2)
sel.clust = "integrated_snn_res.1_c1_sub"
sce1 <- SetIdent(T_seurat2, value = sel.clust)
Seurat::DimPlot(sce1, label = T)

sce2 <- Seurat::FindSubCluster(sce1, cluster = "20", graph.name = "integrated_snn", 
                               subcluster.name = "integrated_snn_res.1_c1_sub", resolution = 0.5)
sel.clust = "integrated_snn_res.1_c1_sub"
sce2 <- SetIdent(sce2, value = sel.clust)
Seurat::DimPlot(sce2, label = T)

sce3 <- Seurat::FindSubCluster(sce2, cluster = "19", graph.name = "integrated_snn", 
                               subcluster.name = "integrated_snn_res.1_c1_sub", resolution = 0.2)
sel.clust = "integrated_snn_res.1_c1_sub"
sce3 <- SetIdent(sce3, value = sel.clust)
Seurat::DimPlot(sce3, label = T)

cells_to_remove <- Seurat::WhichCells(sce3, idents = c(24,25,26,21))
sce3_T <- subset(sce3, cells = setdiff(Cells(sce3), cells_to_remove))

new.cluster.ids <- rep("CD4+", length(levels(sce3_T)))  # Default to CD4+
cd4_clusters <- c(0,2,"3_0",6,7,8,18,"19_0","20_1","20_2",22,23)  # CD4+ clusters
cd8_clusters <- c(1,"3_1","3_2",4,5,9,10,11,12,13,14,15,16,17,"19_1","20_0","20_3",27)  # CD8+ clusters
new.cluster.ids[which(levels(sce3_T) %in% cd8_clusters)] <- "CD8+"
names(new.cluster.ids) <- levels(sce3_T)
sce3_T_rename <-  Seurat::RenameIdents(sce3_T, new.cluster.ids)
Seurat::DimPlot(sce3_T_rename, label = TRUE)


# Neg and iono groups
selected_samples <- c("HC4_Neg","HC5_Neg", "LCKE_Neg","LKKS_Neg")
selected_samples <- c("HC4_iono", "HC5_iono", "LCKE_iono","LKKS_iono")
sce3_T_cd4_re_neg <- subset(sce3_T_rename, subset = orig.ident %in% selected_samples)
sce3_T_cd4_re_iono <- subset(sce3_T_rename, subset = orig.ident %in% selected_samples)
Seurat::DimPlot(sce3_T_rename_neg, label = F, split.by = "orig.ident", cols = c("#FFA500","#4B0082"))
Seurat::DimPlot(sce3_T_rename_iono, label = TRUE, split.by = "orig.ident", cols = c("#FFA500","#4B0082"))


# T cell proportion
sce3_T_rename_p <- plot_celltype_proportions(sce3_T_rename)


# CD4 
sce3_T_cd4 <- subset(sce3_T, idents = c(0,2,"3_0",6,7,8,18,"19_0","20_1","20_2",22,23))
Seurat::DimPlot(sce3_T_cd4, label = TRUE)
# CD8
sce3_T_cd8 <- subset(sce3_T, idents = c(1,"3_1","3_2",4,5,9,10,11,12,13,14,1,16,17,"19_1","20_0","20_3",27))

# recluster
DefaultAssay(sce3_T_cd4) <- "integrated"
sce3_T_cd4_re <- RunPCA(sce3_T_cd4, npc = 20)
sce3_T_cd4_re <- FindNeighbors(sce3_T_cd4_re, dim = 1:20)
sce3_T_cd4_re <- FindClusters(sce3_T_cd4_re,resolution = 0.8)
sce3_T_cd4_re <- RunUMAP(sce3_T_cd4_re, reduction = "pca", dim = 1:20)
Seurat::DimPlot(sce3_T_cd4_re, label = TRUE)
Seurat::ElbowPlot(sce3_T_cd4_cd4)


# marker gene heatmap --------------
Seurat::FeaturePlot(sce3_T_cd4_re, features = c("CD4","CD8A","CCR7","CXCR6","FOXP3","TYMS"), cols = c("lightgrey", "red"),
                    ncol = 3, max.cutoff = 2)

marker <- c("CCR7","IL7R","TCF7","CXCR5","TOX2","BCL6","PDCD1","TOX","ASCL2","PRDM1","CXCR6","TBX21","RUNX3","ID2","IFNG","PRF1","CCL4",
            "GZMB","GZMK","FOXP3","TYMS","MKI67","CLSPN","CD55")

scRNAtoolVis::averageHeatmap(object = sce3_T_cd4_re_rename,htCol = c("steelblue3","white","red3"), border=F,
                             markerGene = marker, fontsize = 10, annoColType = "light",
                             cluster.order = c("Tcm","Tfh","Th1/cytotoxic","Treg","Cycling","Stimulated"))


# DEG and plot
sce3_T_cd4_re_neg_marks <- Seurat::FindAllMarkers(sce3_T_cd4_re_neg, only.pos = TRUE, min.pct = 0.01, logfc.threshold = 0.1)
sce3_T_cd4_re_neg_marks_df <- sce3_T_cd4_re_neg_marks %>% as.data.frame(.) %>% rownames_to_column(., "genes")
write.csv(sce3_T_cd4_re_neg_marks_df, file = "./sce3_T_cd4_re_neg_marks_df.csv", row.names = FALSE)

gene_of_interest <- c("ZBTB7B")
gene_of_interest <- c("RUNX3")
df <- FetchData(sce3_T_cd4_re_neg, vars = c(gene_of_interest, "ident", "orig.ident")) %>%
  rename(gene_expression = !!gene_of_interest) %>%  # Rename gene column
  mutate(orig.ident = factor(orig.ident, levels = c("Controls", "Patients")))
colnames(df) <- c("gene_expression","ident","Group")
ggplot(df, aes(x = Group, y = gene_expression, fill = Group)) +
  geom_jitter(aes(color = Group), width = 0.2, alpha = 0.5, size = 1.5) +  # Jitter for individual points
  geom_boxplot(outlier.shape = NA, alpha = 0.6) +
  facet_wrap(~ident) +  # Separate plots for each cell type
  theme_classic() +
  scale_fill_manual(values = c("Controls" = "#1f77b4", "Patients" = "#ff7f0e")) + 
  scale_color_manual(values = c("Controls" = "#1f77b4", "Patients" = "#ff7f0e")) +
  labs(title = "Gene Expression Across Cell Types",
       x = "Group",
       y = "Gene Expression") +
  theme(legend.position = "top") +
  scale_y_continuous(limits = c(0, 2.9))

# for multiple genes plot
library(dplyr)
library(tidyr)
gene_of_interest <- c("IFNG","NKG7","SLAMF7","CRTAM","PRF1","KLRC1","CCL4","CCL5","GZMA","GZMB","GZMH","GZMM")
df <- FetchData(sce3_T_cd4_re_neg, vars = c(gene_of_interest, "orig.ident")) %>%
  pivot_longer(cols = all_of(gene_of_interest), 
               names_to = "gene", 
               values_to = "gene_expression") %>%
  rename(Group = orig.ident) %>%
  mutate(Group = factor(Group, levels = c("Controls", "Patients")))
ggplot(df, aes(x = Group, y = gene_expression, fill = Group)) +
  geom_jitter(aes(color = Group), width = 0.2, alpha = 0.3, size = 1.5) +  # Jitter for individual points
  geom_violin(alpha = 0.6, trim = FALSE) +  # Violin plot for distribution
  geom_boxplot(width = 0.2, outlier.shape = NA, alpha = 0.6) +  # Boxplot inside violin
  facet_wrap(~gene, scales = "free_y") +  # Facet by gene, allowing independent y scales
  theme_classic() +
  scale_fill_manual(values = c("Controls" = "#1f77b4", "Patients" = "#FFA500")) + 
  scale_color_manual(values = c("Controls" = "#1f77b4", "Patients" = "#FFA500")) +
  labs(title = "Gene Expression Across Cell Types",
       x = "Group",
       y = "Gene Expression") +
  theme(legend.position = "top")

# volcano plot for Tfh
sce3_T_cd4_re_neg_Tfh_marker <- Seurat::FindMarkers(sce3_T_cd4_re_neg, assay = "RNA",ident.1 = "Patients",
                                                    group.by = 'orig.ident', min.pct = 0.01, logfc.threshold = 0.1, subset.ident = "Tfh")
sce3_T_cd4_re_neg_Tcm_marker <- Seurat::FindMarkers(sce3_T_cd4_re_neg, assay = "RNA",ident.1 = "Patients",
                                                    group.by = 'orig.ident', min.pct = 0.01, logfc.threshold = 0.1, subset.ident = "Tcm")
sce3_T_cd4_re_neg_Tfh_marker_df <- sce3_T_cd4_re_neg_Tfh_marker %>% as.data.frame(.) %>% rownames_to_column(., "genes")
sce3_T_cd4_re_neg_Tcm_marker_df <- sce3_T_cd4_re_neg_Tcm_marker %>% as.data.frame(.) %>% rownames_to_column(., "genes")

sce3_T_cd4_re_neg_Tfh_marker_df$diffexpressed <- "NO"
sce3_T_cd4_re_neg_Tfh_marker_df$diffexpressed[sce3_T_cd4_re_neg_Tfh_marker_df$avg_log2FC > 0.5 & sce3_T_cd4_re_neg_Tfh_marker_df$p_val_adj < 0.05] <- "UP"
sce3_T_cd4_re_neg_Tfh_marker_df$diffexpressed[sce3_T_cd4_re_neg_Tfh_marker_df$avg_log2FC < -0.5 & sce3_T_cd4_re_neg_Tfh_marker_df$p_val_adj < 0.05] <- "DOWN"

mycolors <- c("grey", "red", "black")
names(mycolors) <- c("DOWN", "UP", "NO")
sce3_T_cd4_re_neg_Tfh_marker_df$label <- NA
sce3_T_cd4_re_neg_Tfh_marker_df$label[sce3_T_cd4_re_neg_Tfh_marker_df$genes  %in% c("GZMB", "HLA-B", "NKG7", "GZMM", "KLRB1", "KLRK1")] <-
  sce3_T_cd4_re_neg_Tfh_marker_df$genes[sce3_T_cd4_re_neg_Tfh_marker_df$genes %in% c("GZMB", "HLA-B", "NKG7", "GZMM", "KLRB1", "KLRK1")]


library(ggrepel)
volcano_Tfh <- ggplot(data=sce3_T_cd4_re_neg_Tfh_marker_df, 
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
  xlim(c(-2.8, 4)) +
  theme(axis.text=element_text(size=12),
        axis.title=element_text(size=12)) +
  NoLegend()

ggsave(volcano_Tfh,
       filename = "./volcano_Tfh.pdf",
       width = 4, height = 4)

write.csv(sce3_T_cd4_re_neg_Tfh_marker_df, file = "./sce3_T_cd4_re_neg_Tfh_marker_df.csv", row.names = FALSE)
# volcano plot for Tcm
sce3_T_cd4_re_neg_Tcm_marker_df$diffexpressed <- "NO"
sce3_T_cd4_re_neg_Tcm_marker_df$diffexpressed[sce3_T_cd4_re_neg_Tcm_marker_df$avg_log2FC > 0.5 & sce3_T_cd4_re_neg_Tcm_marker_df$p_val_adj < 0.05] <- "UP"
sce3_T_cd4_re_neg_Tcm_marker_df$diffexpressed[sce3_T_cd4_re_neg_Tcm_marker_df$avg_log2FC < -0.5 & sce3_T_cd4_re_neg_Tcm_marker_df$p_val_adj < 0.05] <- "DOWN"

mycolors <- c("grey", "red", "black")
names(mycolors) <- c("DOWN", "UP", "NO")
sce3_T_cd4_re_neg_Tcm_marker_df$label <- NA
sce3_T_cd4_re_neg_Tcm_marker_df$label[sce3_T_cd4_re_neg_Tcm_marker_df$genes  %in% c("GZMB", "HLA-B", "NKG7", "GZMM", "KLRK1")] <-
  sce3_T_cd4_re_neg_Tcm_marker_df$genes[sce3_T_cd4_re_neg_Tcm_marker_df$genes %in% c("GZMB", "HLA-B", "NKG7", "GZMM", "KLRK1")]

volcano_Tcm <- ggplot(data=sce3_T_cd4_re_neg_Tcm_marker_df, 
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
  xlim(c(-2.8, 4)) +
  theme(axis.text=element_text(size=12),
        axis.title=element_text(size=12)) +
  NoLegend()

ggsave(volcano_Tcm,
       filename = "./volcano_Tcm.pdf",
       width = 4, height = 4)


# CD4+ proportion
ce3_T_cd4_re_neg <- subset(sce3_T_cd4_re_neg, idents = c("Tcm","Tfh","Th1/cytotoxic","Treg"))
sce3_T_cd4_re_iono <- subset(sce3_T_cd4_re_iono, idents = c("Tcm","Tfh","Th1/cytotoxic","Treg"))
ce3_T_cd4_re_neg_p <- plot_celltype_proportions(ce3_T_cd4_re_neg)
sce3_T_cd4_re_iono_p <- plot_celltype_proportions(sce3_T_cd4_re_iono)


# proportion function
library(Seurat)
library(ggplot2)
library(dplyr)

plot_celltype_proportions <- function(seurat_obj) {
  seurat_obj@meta.data$orig.ident[seurat_obj@meta.data$orig.ident == "HC4_Neg"]    <- "Controls"
  seurat_obj@meta.data$orig.ident[seurat_obj@meta.data$orig.ident == "HC5_Neg"]    <- "Controls"
  seurat_obj@meta.data$orig.ident[seurat_obj@meta.data$orig.ident == "LCKE_Neg"]   <- "Patients"
  seurat_obj@meta.data$orig.ident[seurat_obj@meta.data$orig.ident == "LKKS_Neg"]   <- "Patients"
  seurat_obj@meta.data$orig.ident[seurat_obj@meta.data$orig.ident == "HC4_iono"]   <- "Controls"
  seurat_obj@meta.data$orig.ident[seurat_obj@meta.data$orig.ident == "HC5_iono"]   <- "Controls"
  seurat_obj@meta.data$orig.ident[seurat_obj@meta.data$orig.ident == "LCKE_iono"]  <- "Patients"
  seurat_obj@meta.data$orig.ident[seurat_obj@meta.data$orig.ident == "LKKS_iono"]  <- "Patients"
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






