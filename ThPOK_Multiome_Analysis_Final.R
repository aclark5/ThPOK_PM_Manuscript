# ---- Single Cell Multiome Analysis of ThPOK point mutant bone marrow chimeras ---- #
# ---- Analysis is of Tetramer Positive CD4+ T cells a day 8pi with LCMV Arm ---- #
# ---- Author of Script: Alex Clark, Date: 03/05/2025---- #

# Load the necessary libraries to perform the analysis
libs <- c(libs <- c(
  "Seurat", "Signac", "ggplot2", "ggrepel", "scDblFinder", "patchwork",
  "dplyr", "tibble", "clustree", "tidyr", "readxl", "purrr", "magrittr",
  "stringr", "readr", "EnsDb.Mmusculus.v79", "BSgenome.Mmusculus.UCSC.mm10",
  "RColorBrewer", "GenomeInfoDb", "viridis", "Polychrome", "harmony", "cowplot",
  "TFBSTools", "motifmatchr", "DBI"
))
lapply(libs, require, character.only = TRUE)

#Ensure we are working with the Seurat V5 object type
options(Seurat.object.assay.version = "v5")

# options(future.globals.maxSize = 1000 * 1024^2)

#Set a seed for reproducibility
set.seed(10)

# ---- Generate custom functtion ---- #

#Function for plotting gene heatmaps
pheatmap.gene.per.cluster <- function(obj.seurat,
                                      gene.to.use,
                                      cluster.used,
                                      color = colorRampPalette( rev(brewer.pal(5, "RdBu")) )(255),
                                      border_color,
                                      fontsize_row = 12,
                                      fontsize_col = 12,
                                      cellwidth = 16,
                                      cellheight = 16,
                                      show_colnames = T,
                                      treeheight_col = 0,
                                      cluster_rows = F,
                                      cluster_cols = F,
                                      order = NULL,
                                      gaps_row = NULL, 
                                      gaps_col = NULL
){
  
  gene.to.use <- intersect(gene.to.use, rownames(obj.seurat))
  df.select.gene <- data.frame(cluster = cluster.used,
                               t(as.matrix(obj.seurat@assays$RNA@layers$data)[gene.to.use,]),
                               check.names = FALSE)
  g <- split(df.select.gene, df.select.gene$cluster)
  expr.gene.per.cluster <- sapply(g, function(x) colMeans(x[, -1]))
  
  if (!is.null(order)){
    m <- na.omit(match(order, colnames(expr.gene.per.cluster)))
    expr.gene.per.cluster <- expr.gene.per.cluster[,m]
  }
  
  pheatmap::pheatmap(expr.gene.per.cluster, 
                     cluster_rows = cluster_rows,
                     cluster_cols = cluster_cols,
                     treeheight_col = treeheight_col,
                     scale = "row",
                     border_color = border_color,
                     color = color,
                     fontsize_row = fontsize_row,
                     fontsize_col = fontsize_col,
                     show_colnames = show_colnames,
                     cellwidth = cellwidth,
                     cellheight = cellheight,
                     name = "Scaled Expression",#,
                     gaps_row = gaps_row, 
                     gaps_col = gaps_col
  )
  
}


#Create a function to convert human to mouse gene names
to_mouse <- function(x) {
  substr(x, 2, nchar(x)) <- tolower(substr(x, 2, nchar(x)))
  x
}

# ----- Set up directory structure to save figures and processed data ----- #

#Create a directory to store the results
project.path <- ("/data/Bosselut_Lab/Alex/ThPOK_PM_LCMV_Multiome/R_jobs/")
setwd(project.path)

ResultsDir <- file.path(project.path,"Results",Sys.Date(), "/")
if (!dir.exists(ResultsDir))
  dir.create(ResultsDir, recursive = T)

# Create a directory to store figures
figDir <- file.path(ResultsDir,"Figures/")
if (!dir.exists(figDir))
  dir.create(figDir, recursive = T)

# Create a subdirectory to store ATAC QC plots
ATAC_figDir <- paste0(figDir,"ATAC_QC_plots/")
if (!dir.exists(ATAC_figDir))
  dir.create(ATAC_figDir, recursive = T)

# Create a subdirectory to store RNA QC plots
RNA_figDir <- paste0(figDir,"RNA_QC_plots/")
if (!dir.exists(RNA_figDir))
  dir.create(RNA_figDir, recursive = T)

# Create a directory to store processed data objects
DataDir <- file.path(ResultsDir,"processed_data/")
if (!dir.exists(DataDir))
  dir.create(DataDir, recursive = T)


#Write the system information and package version to file (helps with future reproducibility)
writeLines(capture.output(sessionInfo()), paste0(ResultsDir,"/sessionInfo_",Sys.Date(),".txt"))

# ----- Read in the RNA and ATAC Datasets ----- #

#Here we have eight datasets (two biological replicates from each genotype (WT/WT, PM/WT, PM/PM, PD)):
#Note RNA data has been processed using cellbender version 0.3.0

RNA_counts <- list()

for (replicate in c("R1", "R2")) {
  
  for (geno in c("wt_wt", "pm_wt", "pm_pm", "pd")) {
    
    RNA_counts[[replicate]][[geno]] <- Read10X_h5(paste0(project.path,
                                                         "counts/", 
                                                         replicate,
                                                         "/multiome_counts_", 
                                                         geno, 
                                                         "/outs/raw_feature_bc_matrix_post_cellbender_filtered_seurat.h5"))
    
  }
}

fragpaths <- list()

for (replicate in c("R1", "R2")) {
  
  for (geno in c("wt_wt", "pm_wt", "pm_pm", "pd")) {
    
    fragpaths[[replicate]][[geno]] <- paste0(project.path,
                                             "counts/", 
                                             replicate,
                                             "/multiome_counts_", 
                                             geno, 
                                             "/outs/atac_fragments.tsv.gz")
    
  }
}

#Get the annotation data in UCSC format
annotation <- GetGRangesFromEnsDb(ensdb = EnsDb.Mmusculus.v79)
seqlevelsStyle(annotation) <- "UCSC"

#Create a Seurat object
obj_all <- list()

for (replicate in c("R1", "R2")) {
  
  for (geno in c("wt_wt", "pm_wt", "pm_pm", "pd")) {
    
    obj_all[[replicate]][[geno]] <- CreateSeuratObject(counts = RNA_counts[[replicate]][[geno]]$`Gene Expression`,
                                                       assay = "RNA")
    
    obj_all[[replicate]][[geno]][["ATAC"]] <- CreateChromatinAssay(counts = RNA_counts[[replicate]][[geno]]$Peaks,
                                                                   sep = c(":", "-"),
                                                                   fragments = fragpaths[[replicate]][[geno]],
                                                                   annotation = annotation)
    
    
    atac_meta <- read_csv(file = paste0(project.path,"counts/", replicate,"/multiome_counts_", geno, "/outs/per_barcode_metrics.csv")) %>%
      dplyr::select("barcode", "atac_fragments", "atac_peak_region_fragments")
    
    obj_all[[replicate]][[geno]]@meta.data <- obj_all[[replicate]][[geno]]@meta.data %>% rownames_to_column("barcode") %>%
      left_join(atac_meta) %>%
      column_to_rownames("barcode")
    
    rm(atac_meta)
  }
}


# ----- Filter to remove empty droplets ----- #

for (replicate in c("R1", "R2")) {
  
  for (geno in c("wt_wt", "pm_wt", "pm_pm", "pd")) {
    
    obj_all[[replicate]][[geno]] <- subset(obj_all[[replicate]][[geno]],
                                           subset = nFeature_RNA > 500)
    
    print(paste0(geno," cells post filter = " , ncol(obj_all[[replicate]][[geno]])))
    
  }
}

# ----- Run scDblFinder to identify doublets ----- #

sce_list <- list()

for (replicate in c("R1", "R2")) {
  
  for (geno in c("wt_wt", "pm_wt", "pm_pm", "pd")) {
    
    sce_list[[replicate]][[geno]] <- SingleCellExperiment(assays = list(counts = obj_all[[replicate]][[geno]]@assays$RNA$counts))
    
    sce_list[[replicate]][[geno]] <- scDblFinder(sce_list[[replicate]][[geno]])
    
    #Add doublet info to object metadata
    obj_all[[replicate]][[geno]]$scDblFinder.class <- sce_list[[replicate]][[geno]]$scDblFinder.class
    
    
  }
}

rm(sce_list)

#Remove the doublets

for (replicate in c("R1", "R2")) {
  
  for (geno in c("wt_wt", "pm_wt", "pm_pm", "pd")) {
    
    obj_all[[replicate]][[geno]] <- subset(obj_all[[replicate]][[geno]],
                                           subset = scDblFinder.class == "singlet")
    
    
    print(paste0(geno," cells post filter = " , ncol(obj_all[[replicate]][[geno]])))
    
  }
}

# -----Perform quality control on the ATAC data and save plots ----- #


for (replicate in c("R1", "R2")) {
  
  for (geno in c("wt_wt", "pm_wt", "pm_pm", "pd")) {
    
    DefaultAssay(obj_all[[replicate]][[geno]]) <- "ATAC"
    
    obj_all[[replicate]][[geno]] <- NucleosomeSignal(obj_all[[replicate]][[geno]])
    obj_all[[replicate]][[geno]] <- TSSEnrichment(obj_all[[replicate]][[geno]])
    obj_all[[replicate]][[geno]]$pct_reads_in_peaks <- (obj_all[[replicate]][[geno]]$atac_peak_region_fragments / obj_all[[replicate]][[geno]]$atac_fragments) * 100
    
    atac_qc <- VlnPlot(object = obj_all[[replicate]][[geno]],
                       features = c("nCount_RNA",
                                    "nCount_ATAC",
                                    "TSS.enrichment",
                                    "nucleosome_signal"),
                       ncol = 4,
                       pt.size = 0)
    
    ggsave(filename = file.path(paste0(ATAC_figDir,replicate, "_", geno,"_RNA_ATAC_QC_metrics.png")),
           height = 6,
           width = 12,
           units = "in",
           dpi=300)
    
    rm(atac_qc)
    
  }
}

# ----- Perform quality control on the RNA data ----- #
for (replicate in c("R1", "R2")) {
  
  for (geno in c("wt_wt", "pm_wt", "pm_pm", "pd")) {
    
    DefaultAssay(obj_all[[replicate]][[geno]]) <- "RNA"
    
    obj_all[[replicate]][[geno]][["percent.mt"]] <- PercentageFeatureSet(obj_all[[replicate]][[geno]],
                                                                         pattern = "^mt-",
                                                                         assay = "RNA")
    
    obj_all[[replicate]][[geno]][["percent.rps"]] <- PercentageFeatureSet(obj_all[[replicate]][[geno]],
                                                                          pattern = "^Rps",
                                                                          assay = "RNA")
    
    obj_all[[replicate]][[geno]][["percent.rpl"]] <- PercentageFeatureSet(obj_all[[replicate]][[geno]],
                                                                          pattern = "^Rpl",
                                                                          assay = "RNA")
    
    rna_qc_1 <- VlnPlot(obj_all[[replicate]][[geno]],
                        features = c("nFeature_RNA",
                                     "nCount_RNA",
                                     "percent.mt",
                                     "percent.rps"),
                        ncol = 4)
    
    ggsave(filename = file.path(paste0(RNA_figDir,replicate, "_",geno,"_RNA_QC_metrics.png")),
           height = 6,
           width = 12,
           units = "in",
           dpi=300)
    
    rna_qc_2 <- CombinePlots(plots = list(FeatureScatter(obj_all[[replicate]][[geno]],
                                                         feature1 = "nCount_RNA",
                                                         feature2 = "percent.mt"),
                                          FeatureScatter(obj_all[[replicate]][[geno]],
                                                         feature1 = "nCount_RNA",
                                                         feature2 = "nFeature_RNA",
                                                         pt.size = 0.1)))
    
    ggsave(filename = file.path(paste0(RNA_figDir,replicate, "_",geno,"_RNA_QC_metrics_comparison.png")),
           height = 6,
           width = 12,
           units = "in",
           dpi=300)
    
    
    png(file = file.path(paste0(RNA_figDir,replicate, "_",geno,"_RNA_count_histo.png")),
        height = 6, width = 9, units = "in", res = 300)
    
    ncount_hist <-  hist(obj_all[[replicate]][[geno]]$nCount_RNA,
                         breaks = 100,
                         xlim = c(0,10000))
    
    print(ncount_hist)
    
    dev.off()
    
    
    png(file = file.path(paste0(RNA_figDir,replicate, "_",geno,"_RNA_feature_histo.png")),
        height = 6, width = 9, units = "in", res = 300)
    
    nfeat_hist <-  hist(obj_all[[replicate]][[geno]]$nFeature_RNA,
                        breaks = 100,
                        xlim = c(0,5000))
    
    print(nfeat_hist)
    
    dev.off()
    
    rm(rna_qc_1,
       rna_qc_2,
       ncount_hist,
       nfeat_hist
    )
    
  }
  
}

# ----- Filter the datasets to remove background or low quality cells ----- #

for (replicate in c("R1", "R2")) {
  
  for (geno in c("wt_wt", "pm_wt", "pm_pm", "pd")) {
    
    cat("\n---", replicate, geno, "---\n")
    cat("Number of cells prior to filtering: ", obj_all[[replicate]][[geno]], "\n")
    
    obj_all[[replicate]][[geno]] <- subset(obj_all[[replicate]][[geno]],
                                           subset = percent.mt < 10 &
                                             nCount_ATAC > 1000 &
                                             pct_reads_in_peaks > 50 &
                                             nucleosome_signal < 2 &
                                             nCount_RNA < 6000 &
                                             TSS.enrichment > 3)
    
    cat("Number of cells after filtering:  ", ncol(obj_all[[replicate]][[geno]]), "\n")
    
  }
}


# ----- Merge the data and normalize and scale ----- #

for (replicate in c("R1", "R2")) {
  
  for (geno in c("wt_wt", "pm_wt", "pm_pm", "pd")) {
    
    obj_all[[replicate]][[geno]]$sample_id <- paste0(replicate,"_",geno)
    
  }
}

obj_all_merged <- merge(obj_all[[1]][[1]], y = c(obj_all[[1]][[2]],
                                                 obj_all[[1]][[3]],
                                                 obj_all[[1]][[4]],
                                                 obj_all[[2]][[1]],
                                                 obj_all[[2]][[2]],
                                                 obj_all[[2]][[3]],
                                                 obj_all[[2]][[4]]))

rm(obj_all, RNA_counts)
gc()

# Log Normalize
obj_all_merged <- NormalizeData(obj_all_merged,
                                normalization.method = "LogNormalize",
                                scale.factor = 10000)

# Identify the top 2000 variable genes
obj_all_merged <- FindVariableFeatures(obj_all_merged,
                                       selection.method = "vst",
                                       nfeatures = 2000)

all.genes <- rownames(obj_all_merged)

# Scale the data and regress out variability associated with mitochondrial % and RNA count per cell
obj_all_merged <- ScaleData(obj_all_merged,
                            vars.to.regress = c("percent.mt","nCount_RNA"),
                            features = all.genes)

rm(obj_all)
gc()


# Perform Principal component analysis
obj_all_merged <- RunPCA(obj_all_merged,
                         features = VariableFeatures(object = obj_all_merged))

# Create new variables in the metadata to identify the replicate and genotype
obj_all_merged$replicate <- sapply(strsplit(obj_all_merged$sample_id, "_"), function(x) x[1])
obj_all_merged$genotype <- substr(obj_all_merged$sample_id, 4, nchar(obj_all_merged$sample_id))

# Reorder the factor levels (useful for plotting etc.)
obj_all_merged$replicate <- factor(obj_all_merged$replicate,
                                   levels = c("R1", "R2"))

obj_all_merged$genotype <- factor(obj_all_merged$genotype,
                                  levels = c("wt_wt", "pm_wt", "pm_pm", "pd"))


# ----- Perform clustering on the Gene Expression data ----- #

# Cluster across a range of clustering resolutions and generate a cluster tree 
# to inspect cluster assignment

obj_all_merged <- FindClusters(obj_all_merged, resolution = c(0.1,0.2,0.3,0.4,0.5,0.6,
                                                              0.7,0.8,0.9, 1, 1.5))

png(filename = paste0(figDir,"/cluster_tree_merged.png"), height = 10, width = 7,
    units = "in", res = 300)
clustree(obj_all_merged)
dev.off()

# Select a clustering resolution of 0.4
obj_all_merged <- FindClusters(obj_all_merged, resolution = c(0.4))

# Identify cluster markers 
# Only include genes expressed in >25% of one cluster and set a logFC threshold of 0.5
obj_all_merged <- JoinLayers(obj_all_merged)
cluster_markers_merged_04 <- FindAllMarkers(obj_all_merged,
                                            only.pos = T,
                                            min.pct = 0.25,
                                            logfc.threshold = 0.5)

# Save the cluster markers as a CSV file
write_csv(cluster_markers_merged_04, file = paste0(ResultsDir, "cluster_markers_res_04.csv"))

# Set the cluster assignmend as the dataset ID
Idents(obj_all_merged) <- obj_all_merged$RNA_snn_res.0.4
gc()

# Include only the first 9 clusters
# Smaller clusters determined to be contaminants (e.g. NK cells)
obj_all_merged <- subset(obj_all_merged, subset = seurat_clusters %in% c(0:8))

# Following the final filtering step - get the QC data for each genotype
write_csv(obj_all_merged@meta.data %>% 
            dplyr::select(replicate, genotype, nCount_RNA, nFeature_RNA, percent.mt,
                          nCount_ATAC,  pct_reads_in_peaks, nucleosome_percentile,
                          TSS.enrichment, TSS.percentile, nFeature_ATAC, nFeature_peaks
            ),
          file = paste0(ResultsDir,"qc_metrics.csv"))


# Add Cell annotations to the dataset based on canonical marker genes
obj_all_merged <- RenameIdents(object = obj_all_merged,
                               '0' = "Th1 1",
                               '1' = "Cytotoxic",
                               '2' = "Tcmp",
                               '3' = "Th1 2",
                               '4' = "Tfh",
                               '5' = "Cycling 1",
                               '6' = "Cycling 2",
                               '7' = "Cycling 3",
                               '8' = "ISC")

obj_all_merged$cell_type <- Idents(obj_all_merged)

# Save the object to file
saveRDS(obj_all_merged, file = paste0(DataDir, "/obj_all_merged_1.rds"))


# ----- Perform Uniform Manifold Approximation and Projection for Dimension Reduction ----- #

DefaultAssay(obj_all_merged) <- "RNA"


obj_all_merged <- FindNeighbors(obj_all_merged, dims = 1:30)

obj_all_merged <- RunUMAP(obj_all_merged,
                          reduction = "pca",
                          dims = 1:50)

# Though we generate a variable 'cell_type' for the cluster annotations, we will use
# a simplified R1:R9 annotation for the publication

Idents(obj_all_merged) <- obj_all_merged$cell_type
obj_all_merged <- RenameIdents(object = obj_all_merged,
                               "Th1 1" = "R1",
                               "Cytotoxic" = "R3",
                               "Tcmp" = "R4",
                               "Th1 2" = "R2",
                               "Tfh" = "R5",
                               "Cycling 1" = "R6",
                               "Cycling 2" = "R7",
                               "Cycling 3" = "R8",
                               "ISC" = "R9")


obj_all_merged$RNA_anno <- Idents(obj_all_merged)

# Plotting for figure 1E

RNA_UMAP <- DimPlot(obj_all_merged,
                    reduction = "umap",
                    group.by = "RNA_anno",
                    split.by = "genotype",
                    label = FALSE,
                    order = F,
                    cols = c("#edac13",
                             "#0bdb27", 
                             "#e314ce", #Tcmp
                             "#4272f5", #Tcmp
                             "#f55347", #Tfh
                             "#2f803a",
                             "#cc0652", 
                             "#1fd8db", 
                             "#9e9899"),
                    alpha = 1,
                    pt.size = 0.01,
                    ncol = 2
) +
  xlab("UMAP-1") +
  ylab("UMAP-2") +
  NoLegend() +
  theme(axis.text=element_text(size=16),
        axis.title=element_text(size=18),
        strip.text.x = element_blank()) +
  ggtitle("")

ggsave(RNA_UMAP,filename = file.path(figDir, "/Fig_1E.pdf"),
       width = 6.5,height = 6.5)


# Perfoem a 'pseudo-bulk' differential expression analysis to identify any 
# differentially expressed genes between the pm/wt and wt/wt cells
Idents(obj_all_merged) <- obj_all_merged$genotype

diff_markers_all_pm_wt <- FindMarkers(obj_all_merged, 
                                      ident.1 = "pm_wt", 
                                      ident.2 = "wt_wt",
                                      min.pct = 0.1,
                                      logfc.threshold = 0
) %>% 
  rownames_to_column("gene_name") %>% 
  mutate(absolute_fc = abs(avg_log2FC)) %>% 
  #Gm47283 is only detected in mutant cells - possibly alignment of mutant ThPOK?
  dplyr::filter(! gene_name %in% c("Gm47283")) %>% 
  dplyr::filter(p_val_adj < 0.05)

# write these results to a csv file
write_csv(diff_markers_all_pm_wt, file = paste0(ResultsDir, "diff_markers_all_pm_wt.csv"))


# Look at expression of key markers by violin plot (all cells per genotype)
# A subset of these will be used for Figure S4A

for (gene in c("Zbtb7b", "Prdm1", "Runx3", "Maf", "Bcl6", "Tbx21",
               "Socs1", "Socs2", "Socs3", "Socs4", "Cish",
               "Ifng", "Nkg7", "Gzmb", "Klrc1", "Klrc3", "Slamf7",
               "Crtam", "Klrb1", "Ccl4", "Ccl5")) {
  
  gene_violin <- VlnPlot(obj_all_merged,
                         features = c(gene),
                         idents = c("wt_wt",
                                    "pm_wt",
                                    "pm_pm"),
                         pt.size = 0,
                         alpha = 0.1,
                         cols = c("grey", "red", "#DD98A1"),
                         group.by = "genotype") +
    NoLegend() +
    xlab("")
  
  ggsave(gene_violin, file = file.path(paste0(figDir, "/",gene,"_FigS4A.pdf")),
         width = 3,height = 3)
  
}

# ----- Plot a heatmap of marker genes split by cluster and genotype ----- #

# Generate a new variable concatenating the cell genotype and cluster
obj_all_merged$geno_rna_cluster <- paste0(obj_all_merged$genotype, 
                                          "_", 
                                          obj_all_merged$RNA_anno)

Idents(obj_all_merged) <- obj_all_merged$geno_rna_cluster

# Set the order for plotting
order.cluster <- c("wt_wt_R1", "wt_wt_R2", "wt_wt_R3", "wt_wt_R4", "wt_wt_R5",
                   "pm_wt_R1", "pm_wt_R2", "pm_wt_R3",  "pm_wt_R4", "pm_wt_R5",
                   "pm_pm_R1", "pm_pm_R2", "pm_pm_R3", "pm_pm_R4", "pm_pm_R5")

# Define which genes to plot in the heatmap
gene.to.use <- c("Prdm1", "Cxcr6", "Tbx21", "Runx3", "Id2", "Ifng",
                 "Ctla2a", "Klrd1", "Ly6c2", "Prf1", "Ccl4", "Gzmb", "Gzmk",
                 "Lef1", "Bach2", "Tcf7", "Ccr7", "Bcl2", "Slamf6",
                 "Cxcr5", "Bcl6", "Pdcd1", "Tox", "Tox2", "Ascl2", "Angptl2")


DefaultAssay(obj_all_merged) <- "RNA"

Idents(obj_all_merged) <- obj_all_merged$cell_type
dimnames(obj_all_merged@assays$RNA@layers$data) <- dimnames(obj_all_merged@assays$RNA)

# Plot the heatmap
pdf(file = paste0(figDir, "/Fig_S4B.pdf"), height =14, width = 8)
pheatmap.gene.per.cluster(obj.seurat = obj_all_merged,
                          gene.to.use = gene.to.use,
                          cluster.used = Idents(obj_all_merged),
                          order = order.cluster,
                          color = colorRampPalette( rev(brewer.pal(5, "RdBu")) )(255),
                          border_color = "white",
                          fontsize_row = 12,
                          fontsize_col = 12,
                          cellwidth = 16,
                          cellheight = 16,
                          show_colnames = F,
                          treeheight_col = 0,
                          gaps_col = c(5, 10),
                          cluster_rows = F,
                          cluster_cols = F)

dev.off()


# ----- Plot the cluster proportions across genotypes ----- #

cluster_counts <- data.frame(cluster = obj_all_merged$RNA_anno,
                             sample_id = obj_all_merged$genotype) %>% 
  dplyr::filter(sample_id %in% c("wt_wt", "pm_wt", "pm_pm", "pd")) %>% 
  # Here will only plot proportions of the 5 largest clusters (Th1 (1&2), Tcmp, Tfh, Cytotoxic)
  dplyr::filter(cluster %in% c("R1", "R2", "R3", "R4", "R5"))

cluster_counts <- table(cluster_counts$cluster,cluster_counts$sample_id) %>% as.data.frame %>% 
  set_colnames(., c("cluster", "genotype","cell_number")) 



cluster_proportions <- cluster_counts %>%
  group_by(genotype) %>%
  mutate(total_count = sum(cell_number)) %>%
  group_by(cluster, add=T) %>%
  mutate(prop=round(cell_number/total_count,5))


cluster_proportions$cluster <- factor(cluster_proportions$cluster,
                                      levels = c("R5", "R4", "R3", "R2", "R1"))

cluster_proportions$genotype <- factor(cluster_proportions$genotype,
                                       levels = c("wt_wt", "pm_wt", "pm_pm", "pd"))

prop_plot <- ggplot(drop_na(cluster_proportions), aes(x=(genotype), y=prop, fill=cluster)) +
  geom_col(position = "fill",
           colour="black") +
  theme_bw() +
  xlab("") +
  ylab("Porportion") +
  theme(legend.position = "bottom",
        legend.key.size = unit(0.5, 'cm'),
        legend.title = element_blank(),
        legend.text = element_text(size=22)) +
  theme(axis.text.x = element_text(angle = 45, vjust = 1, hjust=1, size = 26),
        axis.text.y = element_text(size = 26),
        axis.title = element_text(size = 32)) +
  scale_fill_manual(values = c("#f55347", #Tfh
                               "#4272f5", #Tcmp
                               "#e314ce", #Cytotoxic
                               "#0bdb27", #Th1 (2)
                               "#edac13" # Th1 (1)
  )) +
  NoLegend()

ggsave(prop_plot, file = file.path(paste0(figDir, "/Fig_2E.pdf")),
       width = 4,height = 7)


# ----- Generate feature plots of some key marker genes ----- #

for (gene in c("Pclaf", "Mki67", "Ifit1", "Isg15", "Bcl6", "Prdm1", 
               "Ccr7", "Cxcr5", "Cxcr6", "Cd4", "Cd8a", "Pdcd1", "Tcf7", "Prf1")) {
  
  feature_plot <- FeaturePlot(obj_all_merged, 
                              features = gene,
                              reduction = "umap",
                              order = T,
                              pt.size = 0.1,
                              ncol = 1) &
    scale_color_gradientn(colors = plasma(n = 10, direction = -1)) &
    theme_classic() &
    theme(text = element_text(size=20),
          axis.title = element_blank())
  
  ggsave(feature_plot, filename = paste0(figDir,"/Feature_plot_",gene,"_Fig_S4C.png"),
         height = 3.5, width = 4)
  
  rm(feature_plot)
  
}


# ----- Differential expression between clusters and genotypes ----- #

obj_all_merged$geno_cluster_id <- paste0(obj_all_merged$genotype, 
                                         "_", 
                                         obj_all_merged$RNA_anno)

Idents(obj_all_merged) <- obj_all_merged$geno_cluster_id


# Define genotypes and clusters
genotypes <- c("pm_pm", "pm_wt", "wt_wt")
clusters <- c("R1", "R2", "R3", "R4", "R5")

# Create all genotype pair combinations (excluding identical comparisons)
geno_comparisons <- t(combn(genotypes, 2)) %>% as.data.frame()
colnames(geno_comparisons) <- c("geno1", "geno2")

# Define genotype rank: more 'pm' means higher rank
genotype_rank <- c("pm_pm" = 3, "pm_wt" = 2, "wt_wt" = 1)

# Generate all pairwise genotype comparisons within each cluster
comparison_df <- expand.grid(cluster = clusters, geno1 = genotypes, geno2 = genotypes) %>%
  dplyr::filter(geno1 != geno2) %>%
  rowwise() %>%
  mutate(
    # Order so genotype_1 has more 'pm' than genotype_2
    geno_pair = list(if (genotype_rank[geno1] >= genotype_rank[geno2]) c(geno1, geno2) else c(geno2, geno1)),
    genotype_1 = geno_pair[[1]],
    genotype_2 = geno_pair[[2]],
    ident.1 = paste(genotype_1, cluster, sep = "_"),
    ident.2 = paste(genotype_2, cluster, sep = "_"),
    contrast_id = paste(ident.1, "vs", ident.2, sep = "_")
  ) %>%
  ungroup() %>%
  distinct(contrast_id, .keep_all = TRUE) %>%   # Remove symmetrical duplicates 
  dplyr::select(cluster, genotype_1, genotype_2, ident.1, ident.2, contrast_id)

# Run differential expression for each comparison
all_DE_results <- pmap_dfr(comparison_df, function(cluster, genotype_1, genotype_2, ident.1, ident.2, contrast_id) {
  message("Comparing: ", ident.1, " vs ", ident.2)
  
  res <- tryCatch({
    FindMarkers(obj_all_merged,
                ident.1 = ident.1,
                ident.2 = ident.2,
                min.pct = 0.1,
                logfc.threshold = 0)
  }, error = function(e) {
    warning(paste("Failed comparison:", ident.1, "vs", ident.2, "-", e$message))
    return(NULL)
  })
  
  if (!is.null(res)) {
    res %>%
      rownames_to_column("gene_name") %>%
      mutate(
        cluster = cluster,
        genotype_1 = genotype_1,
        genotype_2 = genotype_2,
        absolute_fc = abs(avg_log2FC)
      ) %>%
      dplyr::filter(!gene_name %in% c("Gm47283"))
  } else {
    NULL
  }
})

write_csv(all_DE_results, file = paste0(ResultsDir,"all_DE_results.csv"))

# ----- Now generate select volcano plots to show differential expression ----- #
# pm/wt vs wt/wt for Tfh
# This will generate the plots for Figure 2F

diff_markers_Tfh_pm_wt <- all_DE_results %>% 
                          dplyr::filter(genotype_1 == "pm_wt" & genotype_2 == "wt_wt" & cluster == "R5")

# add a column of NAs
diff_markers_Tfh_pm_wt$diffexpressed <- "NO"
# if log2Foldchange > 0.6 and pvalue < 0.05, set as "UP" 
diff_markers_Tfh_pm_wt$diffexpressed[diff_markers_Tfh_pm_wt$avg_log2FC > 0.5 & diff_markers_Tfh_pm_wt$p_val_adj < 0.05] <- "UP"
# if log2Foldchange < -0.6 and pvalue < 0.05, set as "DOWN"
diff_markers_Tfh_pm_wt$diffexpressed[diff_markers_Tfh_pm_wt$avg_log2FC < -0.5 & diff_markers_Tfh_pm_wt$p_val_adj < 0.05] <- "DOWN"

#Assign colors to up and down regulated genes
mycolors <- c("grey", "red", "black")
names(mycolors) <- c("DOWN", "UP", "NO")


diff_markers_Tfh_pm_wt$label <- NA

diff_markers_Tfh_pm_wt$label[diff_markers_Tfh_pm_wt$gene_name  %in%
                               c("Crtam", "Slamf7", "Ctla2a", "Ctsw", "Hspa1a", "Hspa1b",
                                 "Mpp7", "Id2", "Vps37b", "Bhlhe40")] <-
  diff_markers_Tfh_pm_wt$gene_name[diff_markers_Tfh_pm_wt$gene_name
                                   %in% c("Crtam", "Slamf7", "Ctla2a", "Ctsw", "Hspa1a", "Hspa1b",
                                          "Mpp7", "Id2", "Vps37b", "Bhlhe40")]

#Generate the volcano plot
volcano_Tfh_pm_wt <- ggplot(data=diff_markers_Tfh_pm_wt, 
                            aes(x=avg_log2FC, 
                                y=-log10(p_val), 
                                col=diffexpressed, 
                                label=label)) +
  geom_point(size = 0.5) + 
  theme_classic() +
  geom_text_repel(size = 6,
                  min.segment.length = 0) +
  scale_color_manual(values=mycolors) +
  xlab("Average Log2 Fold Change") +
  ylab("-log10 p-value") +
  xlim(c(-2.8, 2.8)) +
  theme(axis.text=element_text(size=18),
        axis.title=element_text(size=18)) +
  NoLegend()

ggsave(volcano_Tfh_pm_wt,
       filename = paste0(figDir,"/volcano_plot_Tfh_pm_wt.pdf"),
       width = 4, height = 4)


# pm/wt vs wt/wt for Tcmp

diff_markers_Tcmp_pm_wt <- all_DE_results %>% 
  dplyr::filter(genotype_1 == "pm_wt" & genotype_2 == "wt_wt" & cluster == "R4")

# add a column of NAs
diff_markers_Tcmp_pm_wt$diffexpressed <- "NO"
# if log2Foldchange > 0.6 and pvalue < 0.05, set as "UP" 
diff_markers_Tcmp_pm_wt$diffexpressed[diff_markers_Tcmp_pm_wt$avg_log2FC > 0.5 & diff_markers_Tcmp_pm_wt$p_val_adj < 0.05] <- "UP"
# if log2Foldchange < -0.6 and pvalue < 0.05, set as "DOWN"
diff_markers_Tcmp_pm_wt$diffexpressed[diff_markers_Tcmp_pm_wt$avg_log2FC < -0.5 & diff_markers_Tcmp_pm_wt$p_val_adj < 0.05] <- "DOWN"

#Assign colors to up and down regulated genes
mycolors <- c("grey", "red", "black")
names(mycolors) <- c("DOWN", "UP", "NO")


diff_markers_Tcmp_pm_wt$label <- NA

diff_markers_Tcmp_pm_wt$label[diff_markers_Tcmp_pm_wt$gene_name  %in%
                                c("Hspa1a", "Hspa1b", "Ctla2a", "Slamf7", "Ctsw", 
                                  "Tnfsf11", "Mpp7", "Ntrk3", "Nr4a3")] <-
  diff_markers_Tcmp_pm_wt$gene_name[diff_markers_Tcmp_pm_wt$gene_name
                                    %in% c("Hspa1a", "Hspa1b", "Ctla2a", "Slamf7", "Ctsw", 
                                           "Tnfsf11", "Mpp7", "Ntrk3", "Nr4a3")]
#Generate the volcano plot
volcano_Tcmp_pm_wt <- ggplot(data=diff_markers_Tcmp_pm_wt, 
                             aes(x=avg_log2FC, 
                                 y=-log10(p_val), 
                                 col=diffexpressed, 
                                 label=label)) +
  geom_point(size = 0.5) + 
  theme_classic() +
  geom_text_repel(size = 6,
                  min.segment.length = 0) +
  scale_color_manual(values=mycolors) +
  xlab("Average Log2 Fold Change") +
  ylab("-log10 p-value") +
  xlim(c(-2.8, 2.8)) +
  theme(axis.text=element_text(size=18),
        axis.title=element_text(size=18)) +
  NoLegend()

ggsave(volcano_Tcmp_pm_wt,
       filename = paste0(figDir,"/volcano_plot_Tcmp_pm_wt.pdf"),
       width = 4, height = 4)



# ----- Do the same for pm/pm vs wt/wt comparisons ----- #
# Figures S5A

# pm/pm vs wt/wt for Tfh
diff_markers_Tfh_pm_pm <- all_DE_results %>% 
  dplyr::filter(genotype_1 == "pm_pm" & genotype_2 == "wt_wt" & cluster == "R5")

# add a column of NAs
diff_markers_Tfh_pm_pm$diffexpressed <- "NO"
# if log2Foldchange > 0.6 and pvalue < 0.05, set as "UP" 
diff_markers_Tfh_pm_pm$diffexpressed[diff_markers_Tfh_pm_pm$avg_log2FC > 0.5 & diff_markers_Tfh_pm_pm$p_val_adj < 0.05] <- "UP"
# if log2Foldchange < -0.6 and pvalue < 0.05, set as "DOWN"
diff_markers_Tfh_pm_pm$diffexpressed[diff_markers_Tfh_pm_pm$avg_log2FC < -0.5 & diff_markers_Tfh_pm_pm$p_val_adj < 0.05] <- "DOWN"

#Assign colors to up and down regulated genes
mycolors <- c("grey", "#DD98A1", "black")
names(mycolors) <- c("DOWN", "UP", "NO")


diff_markers_Tfh_pm_pm$label <- NA

diff_markers_Tfh_pm_pm$label[diff_markers_Tfh_pm_pm$gene_name  %in%
                               c("Ccl5", "Nkg7", "Ly6c2", "Ctla2a", "Plac8", "Ctsw", 
                                 "Slamf7", "Klrc1", "Runx3", "Lag3", "Nav2", "Pkp3", "Ddit4")] <-
  diff_markers_Tfh_pm_pm$gene_name[diff_markers_Tfh_pm_pm$gene_name
                                   %in% c("Ccl5", "Nkg7", "Ly6c2", "Ctla2a", "Plac8", "Ctsw", 
                                          "Slamf7", "Klrc1", "Runx3", "Lag3", "Nav2", "Pkp3", "Ddit4")]


#Generate the volcano plot
volcano_Tfh_pm_pm <- ggplot(data=diff_markers_Tfh_pm_pm, 
                            aes(x=avg_log2FC, 
                                y=-log10(p_val), 
                                col=diffexpressed, 
                                label=label)) +
  geom_point(size = 0.5) + 
  theme_classic() +
  geom_text_repel(size = 6,
                  min.segment.length = 0) +
  scale_color_manual(values=mycolors) +
  xlab("Average Log2 Fold Change") +
  ylab("-log10 p-value") +
  xlim(c(-2.8, 2.8)) +
  theme(axis.text=element_text(size=18),
        axis.title=element_text(size=18)) +
  NoLegend()

ggsave(volcano_Tfh_pm_pm,
       filename = paste0(figDir,"/volcano_plot_Tfh_pm_pm.pdf"),
       width = 4, height = 4)


# pm/pm vs wt/wt for Tcmp
diff_markers_Tcmp_pm_pm <- all_DE_results %>% 
  dplyr::filter(genotype_1 == "pm_pm" & genotype_2 == "wt_wt" & cluster == "R4")

# add a column of NAs
diff_markers_Tcmp_pm_pm$diffexpressed <- "NO"
# if log2Foldchange > 0.6 and pvalue < 0.05, set as "UP" 
diff_markers_Tcmp_pm_pm$diffexpressed[diff_markers_Tcmp_pm_pm$avg_log2FC > 0.5 & diff_markers_Tcmp_pm_pm$p_val_adj < 0.05] <- "UP"
# if log2Foldchange < -0.6 and pvalue < 0.05, set as "DOWN"
diff_markers_Tcmp_pm_pm$diffexpressed[diff_markers_Tcmp_pm_pm$avg_log2FC < -0.5 & diff_markers_Tcmp_pm_pm$p_val_adj < 0.05] <- "DOWN"

#Assign colors to up and down regulated genes
mycolors <- c("grey", "#DD98A1", "black")
names(mycolors) <- c("DOWN", "UP", "NO")


diff_markers_Tcmp_pm_pm$label <- NA

diff_markers_Tcmp_pm_pm$label[diff_markers_Tcmp_pm_pm$gene_name  %in%
                                c("Ccl5", "Nkg7", "Ly6c2", "Ctla2a", "Lgals3", "Ifng", "Ctsw",
                                  "Plac8", "Slamf7", "Klrc1", "Klrd1", "Nav2", "Il6ra", "Tle5", "Nfia", "Tmad7", "Ntrk3")] <-
  diff_markers_Tcmp_pm_pm$gene_name[diff_markers_Tcmp_pm_pm$gene_name
                                    %in% c("Ccl5", "Nkg7", "Ly6c2", "Ctla2a", "Lgals3", "Ifng", "Ctsw",
                                           "Plac8", "Slamf7", "Klrc1", "Klrd1", "Nav2", "Il6ra", "Tle5", "Nfia", "Tmad7", "Ntrk3")]

#Generate the volcano plot
volcano_Tcmp_pm_pm <- ggplot(data=diff_markers_Tcmp_pm_pm, 
                             aes(x=avg_log2FC, 
                                 y=-log10(p_val), 
                                 col=diffexpressed, 
                                 label=label)) +
  geom_point(size = 0.5) + 
  theme_classic() +
  geom_text_repel(size = 6,
                  min.segment.length = 0, 
                  max.overlaps = 50) +
  scale_color_manual(values=mycolors) +
  xlab("Average Log2 Fold Change") +
  ylab("-log10 p-value") +
  xlim(c(-2.8, 2.8)) +
  theme(axis.text=element_text(size=18),
        axis.title=element_text(size=18)) +
  NoLegend()

ggsave(volcano_Tcmp_pm_pm,
       filename = paste0(figDir,"/volcano_plot_Tcmp_pm_pm.pdf"),
       width = 4, height = 4)



# ----- Make heatmap of cytotoxic markers split genotype (pseudobulk) ----- #
# For Figure S5B

Idents(obj_all_merged) <- obj_all_merged$geno_cluster_id

dimnames(obj_all_merged@assays$RNA@layers$data) <- dimnames(obj_all_merged@assays$RNA)

gene.to.use <- c("Gzmb", "Gzmk", "Ctla2a", "Ctla2b", "Klrd1", "Ccl4", "Ccl5", "Prf1",
                 "Crtam", "Slamf7")

DefaultAssay(obj_all_merged) <- "RNA"

YellowPurple<-c("purple", "yellow")

# Define clusters
clusters <- c("R1", "R2", "R3", "R4", "R5")

# Loop over each cluster to generate individual heatmaps
for (clust in clusters) {
  
  # Define ordering for this cluster
  order.cluster <- c("wt_wt", "pm_wt", "pm_pm") %>% paste0("_", clust)
  
  # Filter to those identities present in the object
  order.cluster <- order.cluster[order.cluster %in% unique(obj_all_merged$geno_cluster_id)]
  
  # Skip cluster if none of those identities are present
  if (length(order.cluster) == 0) next
  
  message("Generating heatmap for cluster: ", clust)
  
  pdf(file = paste0(figDir, "/heatmap_cyto_sig_", clust, ".pdf"), height = 5, width = 4)
  
  pheatmap.gene.per.cluster(
    obj.seurat = obj_all_merged,
    gene.to.use = gene.to.use,
    cluster.used = Idents(obj_all_merged),
    order = order.cluster,
    color = colorRampPalette(YellowPurple)(255),
    border_color = "white",
    fontsize_row = 12,
    fontsize_col = 12,
    cellwidth = 16,
    cellheight = 12,
    show_colnames = TRUE,
    treeheight_col = 0,
    cluster_rows = FALSE,
    cluster_cols = FALSE
  )
  
  dev.off()
}




# ----- Now perform analysis on the ATAC data ----- #

PeakDir <- file.path(ResultsDir,"peaks/")
if (!dir.exists(PeakDir))
  dir.create(PeakDir, recursive = T)

peaks <- CallPeaks(obj_all_merged,
                   assay = "ATAC",
                   effective.genome.size = 2652783500,
                   macs2.path = "/usr/local/apps/macs/2.2.7.1/bin/macs2",
                   outdir = PeakDir,
                   broad = F)


peaks <- keepStandardChromosomes(peaks, pruning.mode = "coarse")
peaks <- subsetByOverlaps(x = peaks, ranges = blacklist_mm10, invert = TRUE)


macs2_counts <- FeatureMatrix(fragments = Fragments(obj_all_merged@assays$ATAC),
                              features = peaks,
                              cells = colnames(obj_all_merged))

# create a new assay using the MACS2 peak set and add it to the Seurat object
obj_all_merged[["peaks"]] <- CreateChromatinAssay(
  counts = macs2_counts,
  fragments = Fragments(obj_all_merged@assays$ATAC),
  annotation = annotation
)


DefaultAssay(obj_all_merged) <- "peaks"
obj_all_merged <- FindTopFeatures(obj_all_merged, min.cutoff = 5)
obj_all_merged <- RunTFIDF(obj_all_merged)
obj_all_merged <- RunSVD(obj_all_merged)

# first compute the GC content for each peak
obj_all_merged <- RegionStats(obj_all_merged, genome = BSgenome.Mmusculus.UCSC.mm10)

# Save the object to file
saveRDS(obj_all_merged, file = paste0(DataDir, "/obj_all_merged_2.rds"))


# ----- Perform Cluster analysis on the ATAC data ----- #

obj_all_merged <- FindNeighbors(object = obj_all_merged, 
                                reduction = 'lsi', 
                                dims = 2:30)

obj_all_merged <- FindClusters(object = obj_all_merged, 
                               verbose = FALSE, 
                               algorithm = 3, 
                               resolution = c(0.1,0.2,0.3,0.4,0.5,0.6,0.7,0.8,0.9, 1, 1.5))

#Generate a cluster tree
png(filename = paste0(figDir,"/cluster_tree_peaks_merged.png"), height = 10, width = 7,
    units = "in", res = 300)
clustree(obj_all_merged)
dev.off()


#check correspondance between ATAC and RNA clusters
table(obj_all_merged$RNA_anno, obj_all_merged$peaks_snn_res.0.5)

Idents(obj_all_merged) <- "peaks_snn_res.0.5"

obj_all_merged <- RenameIdents(object = obj_all_merged,
                               '0' = "A1", #R1 Th1 1
                               '1' = "A5", #R4 Tcmp
                               '2' = "A2", #R1/R3 Th1/Cytotoxic
                               '3' = "A4", #R3 Cytotoxic
                               '4' = "A3", #R2 Th1 2
                               '5' = "A6", #R5 Tfh
                               '6' = "A7", #R1/R3 Th1/Cytotoxic
                               '7' = "A8", #R3 Cytotoxic
                               '8' = "A9", #Cycling
                               '9' = "A10", #mix
                               '10' = "A11" #mostly Th1
)


obj_all_merged$ATAC_anno <- Idents(obj_all_merged)



# Make UMAP plots showing the ATAC and RNA clusters split by genotype
# Figure S4D

obj_all_merged <- RunUMAP(obj_all_merged, 
                          assay = "peaks",
                          reduction = "lsi", 
                          dims = 2:30, 
                          reduction.name = "umap.peak", 
                          reduction.key = "atacPEAK_")

p1 <- DimPlot(obj_all_merged,
              reduction = "umap.peak",
              group.by = "ATAC_anno",
              split.by = "genotype",
              label = FALSE,
              order = F,
              cols = c("#edac13",
                       "#f593c1",
                       "#0bdb27", 
                       "#e314ce", #Tcmp
                       "#4272f5", #Tcmp
                       "#f55347", #Tfh
                       "#b87ef2",
                       "#b0ad12",
                       "#2f803a",
                       "#9c5a40",
                       "black"),
              alpha = 1,
              pt.size = 0.1,
              ncol = 4
) +
  xlab("UMAP-1") +
  ylab("UMAP-2") +
  NoLegend() +
  theme(axis.text=element_text(size=20),
        axis.title=element_text(size=22),
        strip.text.x = element_blank()) +
  ggtitle("")


ggsave(p1,filename = file.path(figDir, "/UMAP_by_ATAC_cluster_genotype.png"),
       width = 14,height = 4)


p2 <- DimPlot(obj_all_merged,
              reduction = "umap.peak",
              group.by = "RNA_anno",
              split.by = "genotype",
              label = FALSE,
              order = F,
              cols = c("#edac13",
                       "#0bdb27", 
                       "#e314ce", 
                       "#4272f5", #Tcmp
                       "#f55347", #Tfh
                       "#2f803a",
                       "#cc0652", 
                       "#1fd8db", 
                       "#9e9899"),
              alpha = 1,
              pt.size = 0.01,
              ncol = 4
) +
  xlab("UMAP-1") +
  ylab("UMAP-2") +
  NoLegend() +
  theme(axis.text=element_text(size=20),
        axis.title=element_text(size=22),
        strip.text.x = element_blank()) +
  ggtitle("")


ggsave(p2,filename = file.path(figDir, "/UMAP_by_RNA_cluster_genotype.png"),
       width = 14,height = 4)


#Generate a heatmap to look at correspondence between ATAC and RNA clusters
# Figure S4E

correspondence <- table(obj_all_merged$RNA_anno, obj_all_merged$ATAC_anno)
correspondence <- correspondence/rowSums(correspondence)  # normalize for number of cells in each cell type
correspondence <- as.data.frame(correspondence) 

correspondence$Var1 <- factor(correspondence$Var1,
                              levels = c(paste0("R",1:9)))

correspondence$Var2 <- factor(correspondence$Var2,
                              levels = c(paste0("A",1:11)))


p3 <- ggplot(correspondence, aes(Var1, Var2, fill = Freq)) +
  geom_tile(color="black", height=0.9, width=0.9) +
  scale_fill_gradientn(name = "Fraction \n of cells",
                       colours = c("yellow", "orange", "red")) +
  xlab("Cell type annotation (RNA)") +
  ylab("Cell type annotation (ATAC)") +
  theme_cowplot() +
  theme(axis.text.x = element_text(angle = 45, 
                                   vjust = 1, 
                                   hjust=1))
ggsave(p3, file = file.path(paste0(figDir, "/cell_type_cluster_correspondence.png")),
       width = 5.5,height = 5)


# ----- Differential accessibility between clusters and genotypes ----- #

obj_all_merged$geno_ATAC_cluster_id <- paste0(obj_all_merged$genotype, 
                                         "_", 
                                         obj_all_merged$ATAC_anno)

Idents(obj_all_merged) <- obj_all_merged$geno_ATAC_cluster_id


# Define genotypes and clusters
genotypes <- c("pm_pm", "pm_wt", "wt_wt")
clusters <- c("A1", "A2", "A3", "A4", "A5", "A6")

# Create all genotype pair combinations (excluding identical comparisons)
geno_comparisons <- t(combn(genotypes, 2)) %>% as.data.frame()
colnames(geno_comparisons) <- c("geno1", "geno2")

# Define genotype rank: more 'pm' means higher rank
genotype_rank <- c("pm_pm" = 3, "pm_wt" = 2, "wt_wt" = 1)

# Generate all pairwise genotype comparisons within each cluster
comparison_df_atac <- expand.grid(cluster = clusters, geno1 = genotypes, geno2 = genotypes) %>%
  dplyr::filter(geno1 != geno2) %>%
  rowwise() %>%
  mutate(
    # Order so genotype_1 has more 'pm' than genotype_2
    geno_pair = list(if (genotype_rank[geno1] >= genotype_rank[geno2]) c(geno1, geno2) else c(geno2, geno1)),
    genotype_1 = geno_pair[[1]],
    genotype_2 = geno_pair[[2]],
    ident.1 = paste(genotype_1, cluster, sep = "_"),
    ident.2 = paste(genotype_2, cluster, sep = "_"),
    contrast_id = paste(ident.1, "vs", ident.2, sep = "_")
  ) %>%
  ungroup() %>%
  distinct(contrast_id, .keep_all = TRUE) %>%   # Remove symmetrical duplicates 
  dplyr::select(cluster, genotype_1, genotype_2, ident.1, ident.2, contrast_id)

# Run differential expression for each comparison
all_DA_results <- pmap_dfr(comparison_df_atac, function(cluster, genotype_1, genotype_2, ident.1, ident.2, contrast_id) {
  message("Comparing: ", ident.1, " vs ", ident.2)
  
  res <- tryCatch({
    FindMarkers(obj_all_merged,
                test.use = "LR",
                ident.1 = ident.1,
                ident.2 = ident.2,
                min.pct = 0.1,
                only.pos = F)
  }, error = function(e) {
    warning(paste("Failed comparison:", ident.1, "vs", ident.2, "-", e$message))
    return(NULL)
  })
  
  if (!is.null(res)) {
    res %>%
      rownames_to_column("peak") %>%
      mutate(
        cluster = cluster,
        genotype_1 = genotype_1,
        genotype_2 = genotype_2,
        absolute_fc = abs(avg_log2FC)
      ) 
  } else {
    NULL
  }
})


closest_gene <- ClosestFeature(obj_all_merged,
                               regions = all_DA_results$gene_name) %>% 
  dplyr::select("query_region", "gene_name") %>% 
  set_colnames(., c("peak", "gene_name"))

all_DA_results <- left_join(all_DA_results,
                            closest_gene)

write_csv(all_DA_results, file = paste0(ResultsDir,"all_DA_results.csv"))

# Now generate volcano plots of the differentiall accessible peaks
# labelled with tho closest gene name
# Figure S5C

# pm/pm vs wt/wt for Tfh
diff_peaks_Tfh_pm_wt <- all_DA_results %>% 
  dplyr::filter(genotype_1 == "pm_wt" & genotype_2 == "wt_wt" & cluster == "A6")

# add a column of NAs
diff_peaks_Tfh_pm_wt$diffexpressed <- "NO"
# if log2Foldchange > 0.6 and pvalue < 0.05, set as "UP" 
diff_peaks_Tfh_pm_wt$diffexpressed[diff_peaks_Tfh_pm_wt$avg_log2FC > 0.5 & diff_peaks_Tfh_pm_wt$p_val_adj < 0.05] <- "UP"
# if log2Foldchange < -0.6 and pvalue < 0.05, set as "DOWN"
diff_peaks_Tfh_pm_wt$diffexpressed[diff_peaks_Tfh_pm_wt$avg_log2FC < -0.5 & diff_peaks_Tfh_pm_wt$p_val_adj < 0.05] <- "DOWN"

#Assign colors to up and down regulated genes
mycolors <- c("grey", "red", "black")
names(mycolors) <- c("DOWN", "UP", "NO")


diff_peaks_Tfh_pm_wt$label <- NA

diff_peaks_Tfh_pm_wt$label[diff_peaks_Tfh_pm_wt$gene_name  %in%
                                        c("Crtam")] <-
  diff_peaks_Tfh_pm_wt$gene_name[diff_peaks_Tfh_pm_wt$gene_name
                                            %in% c("Crtam")]

volcano_Tfh_peaks_pm_wt <- ggplot(data=diff_peaks_Tfh_pm_wt, 
                                  aes(x=avg_log2FC, 
                                      y=-log10(p_val), 
                                      col=diffexpressed, 
                                      label=label)
) +
  geom_point(size = 0.5) + 
  theme_classic() +
  geom_text_repel(size = 6,
                  min.segment.length = 0,
                  max.overlaps = 40) +
  scale_color_manual(values=mycolors) +
  xlab("Average Log2 Fold Change") +
  ylab("-log10 p-value") +
  xlim(c(-2.8, 2.8)) +
  theme(axis.text=element_text(size=18),
        axis.title=element_text(size=18)) +
  NoLegend()

ggsave(volcano_Tfh_peaks_pm_wt,
       filename = paste0(figDir,"/volcano_plot_Tfh_peaks_pm_wt.pdf"),
       width = 4, height = 4)


# pm/pm vs wt/wt for Tcmp
diff_peaks_Tcmp_pm_wt <- all_DA_results %>% 
  dplyr::filter(genotype_1 == "pm_wt" & genotype_2 == "wt_wt" & cluster == "A5")

# add a column of NAs
diff_peaks_Tcmp_pm_wt$diffexpressed <- "NO"
# if log2Foldchange > 0.6 and pvalue < 0.05, set as "UP" 
diff_peaks_Tcmp_pm_wt$diffexpressed[diff_peaks_Tcmp_pm_wt$avg_log2FC > 0.5 & diff_peaks_Tcmp_pm_wt$p_val_adj < 0.05] <- "UP"
# if log2Foldchange < -0.6 and pvalue < 0.05, set as "DOWN"
diff_peaks_Tcmp_pm_wt$diffexpressed[diff_peaks_Tcmp_pm_wt$avg_log2FC < -0.5 & diff_peaks_Tcmp_pm_wt$p_val_adj < 0.05] <- "DOWN"

#Assign colors to up and down regulated genes
mycolors <- c("grey", "red", "black")
names(mycolors) <- c("DOWN", "UP", "NO")


diff_peaks_Tcmp_pm_wt$label <- NA

diff_peaks_Tcmp_pm_wt$label[diff_peaks_Tcmp_pm_wt$gene_name  %in%
                             c("Ctla2a", "Crtam", "Prf1",  "Inpp4b", "Tbx21",
                               "Ikzf3", "Snd1", "Mpp7", "Itgav")] <-
  diff_peaks_Tcmp_pm_wt$gene_name[diff_peaks_Tcmp_pm_wt$gene_name
                                 %in% c("Ctla2a", "Crtam", "Prf1",  "Inpp4b", "Tbx21",
                                        "Ikzf3", "Snd1", "Mpp7", "Itgav")]

volcano_Tcmp_peaks_pm_wt <- ggplot(data=diff_peaks_Tcmp_pm_wt, 
                                  aes(x=avg_log2FC, 
                                      y=-log10(p_val), 
                                      col=diffexpressed, 
                                      label=label)
) +
  geom_point(size = 0.5) + 
  theme_classic() +
  geom_text_repel(size = 6,
                  min.segment.length = 0,
                  max.overlaps = 40) +
  scale_color_manual(values=mycolors) +
  xlab("Average Log2 Fold Change") +
  ylab("-log10 p-value") +
  xlim(c(-2.8, 2.8)) +
  theme(axis.text=element_text(size=18),
        axis.title=element_text(size=18)) +
  NoLegend()

ggsave(volcano_Tcmp_peaks_pm_wt,
       filename = paste0(figDir,"/volcano_plot_Tcmp_peaks_pm_wt.pdf"),
       width = 4, height = 4)



# Perform a pseudo-bulk analysis of ThPOK deficient vs wt peaks
# This will be exported and used to produce figure 2G in prism
Idents(obj_all_merged) <- obj_all_merged$genotype

DefaultAssay(obj_all_merged) <- "peaks"

cluster_markers_peaks_PD <- FindMarkers(obj_all_merged,
                                        test.use = "LR",
                                        ident.1 = "pd",
                                        ident.2 = "wt_wt",
                                        only.pos = T,
                                        min.pct = 0.1) %>% 
  dplyr::filter(p_val_adj < 0.05) %>% 
  dplyr::filter(avg_log2FC > 0.5) %>% 
  rownames_to_column("peak")

closest_gene_PD <- ClosestFeature(obj_all_merged,
                               regions = cluster_markers_peaks_PD$peak) %>% 
  dplyr::select("query_region", "gene_name") %>% 
  set_colnames(., c("peak", "gene_name"))

cluster_markers_peaks_PD <- left_join(cluster_markers_peaks_PD,
                                      closest_gene_PD)


# Generate a heatmap of peak accessibility scores across clusters (use top 500 peaks per cluster)

peak.to.use <- cluster_markers_peaks_PD %>% 
  slice_head(n=500)

write_csv(peak.to.use, file = paste0(ResultsDir,"ThPOK_dependent_peaks.csv"))

peak.to.use <- intersect(peak.to.use$peak, rownames(obj_all_merged@assays$peaks))
df.select.peak <- data.frame(t(as.matrix(obj_all_merged@assays$peaks$data)[peak.to.use,]))

cell_genotypes <- as.data.frame(obj_all_merged@meta.data) %>% 
  rownames_to_column("barcode") %>% 
  dplyr::select(c("barcode", "genotype"))

df.select.peak <- df.select.peak %>% 
  rownames_to_column("barcode") %>% 
  left_join(cell_genotypes) %>% 
  column_to_rownames("barcode")

# #Convert 0 to NA
# df.select.gene[df.select.gene == 0] <- NA

df.select.peak <- df.select.peak %>% 
                  mutate(avg_accessibility = rowMeans(df.select.gene[,1:500], na.rm=T))


rm(cell_genotypes)

write.csv(df.select.peak[,c("genotype", "avg_accessibility")], file = paste0(ResultsDir,"avg_PD_accessibility.csv"))


# Generate coverage plots of ATAC peaks across gene regions

obj_all_merged$geno_ATAC_cluster_id <- factor(obj_all_merged$geno_ATAC_cluster_id,
                                              levels = c("wt_wt_A6",
                                                         "pm_wt_A6",
                                                         "pm_pm_A6",
                                                         "wt_wt_A5",
                                                         "pm_wt_A5",
                                                         "pm_pm_A5"))

Idents(obj_all_merged) <- obj_all_merged$geno_ATAC_cluster_id

# Crtam in Tfh (Fig 1H)
crtam_tfh <- CoveragePlot(
  object = obj_all_merged,
  region = c("chr9-40985000-41007000"),
  features = "Crtam",
  group.by = "geno_ATAC_cluster_id",
  idents = c("wt_wt_A6",
             "pm_wt_A6",
             "pm_pm_A6"),
  expression.assay = "RNA",
) &
  scale_fill_manual(values = c("grey", "#DD98A1", "red")) &
  theme(panel.grid.major.y = element_line(),
        panel.grid.minor.y = element_line())

ggsave(crtam_tfh, file = file.path(figDir, "region_plot_Crtam_Tfh.png"),
       width = 6.5,
       height = 5)


# Ctla2a Tcmp (Figure S5D)
ctla2a_tcmp <- CoveragePlot(
  object = obj_all_merged,
  region = c("Ctla2a"),
  features = "Ctla2a",
  extend.upstream = 1000,
  extend.downstream = 1000,
  group.by = "geno_ATAC_cluster_id",
  idents = c("wt_wt_A5",
             "pm_wt_A5",
             "pm_pm_A5"),
  expression.assay = "RNA",
) &
  scale_fill_manual(values = c("grey", "#DD98A1", "red")) &
  theme(panel.grid.major.y = element_line(),
        panel.grid.minor.y = element_line())

ggsave(ctla2a_tcmp, file = file.path(figDir, "region_plot_Ctla2a_Tcmp.png"),
       width = 6.5,
       height = 5)



# ----- Perform analysis of "Gene Activity" based on peak accessibility across gene region ----- #

#Calculate gene activities for each gene
DefaultAssay(obj_all_merged) <- "peaks"
gene.activities <- GeneActivity(obj_all_merged)

# add the gene activity matrix to the Seurat object as a new assay and normalize it
obj_all_merged[['Activity']] <- CreateAssayObject(counts = gene.activities)

obj_all_merged <- NormalizeData(
  object = obj_all_merged,
  assay = 'Activity',
  normalization.method = 'LogNormalize',
  scale.factor = median(obj_all_merged$nCount_Activity)
)

DefaultAssay(obj_all_merged) <- "Activity"

# Calculate gene activity scores for Tfh and Tcmp marker genes
# These are exported to plot Fig S5F with prism

# Read in the cluster marker genes
cluster_markers_res_04 <- read.csv(paste0(ResultsDir,"/cluster_markers_res_04.csv"))

cluster_markers_res_04 <- cluster_markers_res_04 %>%
  dplyr::filter(cluster %in% c(0:8)) %>% 
  mutate(cell_markers = "",
         cell_markers = ifelse(cluster == 0, "R1", cell_markers),
         cell_markers = ifelse(cluster == 1, "R3", cell_markers),
         cell_markers = ifelse(cluster == 2, "R4", cell_markers),
         cell_markers = ifelse(cluster == 3, "R2", cell_markers),
         cell_markers = ifelse(cluster == 4, "R5", cell_markers),
         cell_markers = ifelse(cluster == 5, "R6", cell_markers),
         cell_markers = ifelse(cluster == 6, "R7", cell_markers),
         cell_markers = ifelse(cluster == 7, "R8", cell_markers),
         cell_markers = ifelse(cluster == 8, "R9", cell_markers))

# Select the Tfh genes form the data frame
Tfh_genes <-  cluster_markers_res_04 %>% 
  dplyr::filter(cell_markers == "R5") %>% 
  .[,"gene"]


obj_all_merged <- AddModuleScore(obj_all_merged,
                                 features = list(Tfh_genes),
                                 assay = "Activity",
                                 name = "Tfh_gene_activity_score"
)


Tfh_activities <- as.data.frame(obj_all_merged@meta.data) %>% 
  dplyr::select("geno_cluster_id", "geno_ATAC_cluster_id", "Tfh_gene_activity_score1")

write_csv(Tfh_activities, file = paste0(ResultsDir,"Tfh_gene_activities.csv"))


# Do the same for the Tcmp signature genes

Tcmp_genes <-  cluster_markers_res_04 %>% 
  dplyr::filter(cell_markers == "R4") %>% 
  .[,"gene"]


obj_all_merged <- AddModuleScore(obj_all_merged,
                                 features = list(Tcmp_genes),
                                 assay = "Activity",
                                 name = "Tcmp_gene_activity_score"
)


Tcmp_activities <- as.data.frame(obj_all_merged@meta.data) %>% 
  dplyr::select("geno_cluster_id", "geno_ATAC_cluster_id", "Tcmp_gene_activity_score1")

write_csv(Tcmp_activities, file = paste0(ResultsDir,"Tcmp_gene_activities.csv"))





# -----Run Chromvar to determine TF "activity"----- #


# First need to get a list of motif position frequency matrices from the JASPAR database
# Here I have downloaded the JASPAR database locally

JASPAR2024_DB <- dbConnect(drv=RSQLite::SQLite(), 
                           dbname="/Users/clarkald/OneDriveNIH/Multiome_Data_All/TF_motif_datasets/JASPAR2024.sqlite")


motifs24 <- TFBSTools::getMatrixSet(JASPAR2024_DB, 
                                    opts = list(collection = "CORE", 
                                                tax_group = 'vertebrates', 
                                                all_versions = FALSE))

# add motif information to object
obj_all_merged <- AddMotifs(
  object = obj_all_merged,
  genome = BSgenome.Mmusculus.UCSC.mm10,
  pfm = motifs24
)



# Now we can perform ChromVAR analysis using these PWMs

DefaultAssay(obj_all_merged) <- "peaks"

obj_all_merged <- RunChromVAR(
  object = obj_all_merged,
  genome = BSgenome.Mmusculus.UCSC.mm10
)

#Create a data frame of the TF names (can match with TF IDs)
tf_names <- list()
tf_ids  <- list()
tf_family <- list()

for(factor in 1:length(motifs24)) {
  
  check <- motifs24[[factor]]
  tf_names[factor] <- check@name
  tf_ids[factor] <- check@ID
  tf_family[factor] <- check@matrixClass
  
}

tf_names <- unlist(tf_names) 
tf_ids <- unlist(tf_ids) 
tf_family <- unlist(tf_family)

tfs <- data.frame(tf_names, tf_ids, tf_family)
names(tfs) <- c("tf_name", "tf_id", "tf_family")

DefaultAssay(obj_all_merged) <- "chromvar"

# Obtain ChromVAR values for TCF1
# Figure S5G

TCF1_ID <- "MA0769.3"
TCF1_chromvar <- data.frame((as.matrix(obj_all_merged@assays$chromvar$data)[TCF1_ID,]))

cell_genotypes <- as.data.frame(obj_all_merged@meta.data) %>% 
  rownames_to_column("barcode") %>% 
  dplyr::select(c("barcode", "geno_ATAC_cluster_id"))

TCF1_chromvar <- TCF1_chromvar %>% 
  rownames_to_column("barcode") %>% 
  left_join(cell_genotypes) %>% 
  set_colnames(., c("barcode", TCF1_ID,"ATAC_cluster"))


write_csv(TCF1_chromvar, file = paste0(ResultsDir, "TCF1_chromvar.csv"))

# Obtain ChromVAR values for T-bet
# Figure S5G

TBET_ID <- "MA0690.3"
TBET_chromvar <- data.frame((as.matrix(obj_all_merged@assays$chromvar$data)[TBET_ID,]))

cell_genotypes <- as.data.frame(obj_all_merged@meta.data) %>% 
  rownames_to_column("barcode") %>% 
  dplyr::select(c("barcode", "geno_ATAC_cluster_id"))

TBET_chromvar <- TBET_chromvar %>% 
  rownames_to_column("barcode") %>% 
  left_join(cell_genotypes) %>% 
  set_colnames(., c("barcode", TBET_ID,"ATAC_cluster"))


write_csv(TBET_chromvar, file = paste0(ResultsDir, "TBET_chromvar.csv"))
