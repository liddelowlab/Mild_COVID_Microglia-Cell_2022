# This script provides the code used in the analysis of the mouse mild COVID single-cell RNA-sequencing data from Fernández-Castañeda, Lu, Geraghty, & Song et al. 2022
# All code was run using R version 4.1.1; all packages and dependencies along with their versions can be found in the packages.txt file.
# load required libraries
library(ggplot2)
library(dplyr)
library(Seurat)
library(harmony)
library(ggpubr)
library(clustree)
library(speckle)
library(ggbreak)
library(ComplexHeatmap)
library(EnhancedVolcano)
library(clusterProfiler)
library(DOSE)
library(simplifyEnrichment)
library(UpSetR)
library(ggcorrplot)
library(openxlsx)
library(ggsci)
library(khroma)
library(reshape2)
library(scales)
library(MetBrewer)
library(RColorBrewer)
library(circlize)
library(readxl)
library(stringr)

# enable multiprocess with future package
library(future)
plan("multiprocess")
options(future.globals.maxSize= 10000*1024^2)

# save package versions list:
pkgs <- loadedNamespaces()
for (i in 1:length(pkgs)) {
  cat(print(paste(pkgs[i], getNamespaceVersion(pkgs[i]), sep = " ")), file = "packages.txt", append = TRUE)
  cat("\n", file = "packages.txt", append = TRUE)
}

# load Cell Ranger output files and create seurat object
sample1 <- Read10X("relative/path")
sample2 <- Read10X("relative/path")
sample3 <- Read10X("relative/path")
sample4 <- Read10X("relative/path")
sample5 <- Read10X("relative/path")
sample6 <- Read10X("relative/path")
sample7 <- Read10X("relative/path")
sample8 <- Read10X("relative/path")

sample1 <- CreateSeuratObject(sample1)
sample2 <- CreateSeuratObject(sample2)
sample3 <- CreateSeuratObject(sample3)
sample4 <- CreateSeuratObject(sample4)
sample5 <- CreateSeuratObject(sample5)
sample6 <- CreateSeuratObject(sample6)
sample7 <- CreateSeuratObject(sample7)
sample8 <- CreateSeuratObject(sample8)

### add sample metadata
# sample IDs
sample1$orig.ident <- "monje1"
sample2$orig.ident <- "monje2"
sample3$orig.ident <- "monje3"
sample4$orig.ident <- "monje4"
sample5$orig.ident <- "monje5"
sample6$orig.ident <- "monje6"
sample7$orig.ident <- "monje7"
sample8$orig.ident <- "monje8"

# Mild COVID or control
sample1$orig.ident <- "control"
sample2$orig.ident <- "control"
sample3$orig.ident <- "control"
sample4$orig.ident <- "control"
sample5$orig.ident <- "COVID"
sample6$orig.ident <- "COVID"
sample7$orig.ident <- "COVID"
sample8$orig.ident <- "COVID"

monje <- merge(sample1, c(sample2, sample3, sample4, sample5, sample6, sample7, sample8))

# calculate percent UMIs mapping to mitochondrial genome
monje[["percent.mt"]] <- PercentageFeatureSet(monje, pattern = "^mt-")

# filter cells based on QC metrics
monje <- subset(monje, subset = nFeature_RNA > 200 & nFeature_RNA < 5000 & percent.mt < 10)

# Normalize data
monje <- NormalizeData(monje)

# Find top 2000 variable features
monje <- FindVariableFeatures(monje, selection.method = "vst", nfeatures = 2000)

# Scale data
all.genes <- rownames(monje)
monje <- ScaleData(monje, features = all.genes)

# run PCA
monje <- RunPCA(monje, features = VariableFeatures(object = monje))

# cluster cells
monje <- FindNeighbors(monje, dims = 1:20)
monje <- FindClusters(monje, resolution = 0.8) # 28 clusters; we'll manually group these together into cell types based on marker gene expression 

# run UMAP
monje <- RunUMAP(monje, dims = 1:20)

# plot UMAP colored by clusters
DimPlot(monje, group.by = "RNA_snn_res.0.8")

# set cell identities as cluster labels
Idents(monje) <- monje$RNA_snn_res.0.8

# run differential expression testing
cluster.markers <- FindAllMarkers(monje)
cluster.markers <- cluster.markers %>% filter(p_val_adj < 0.05) %>% arrange(cluster, desc(avg_log2FC))

# save DE results
saveRDS(cluster.markers, "original_cluster_markers.rds")

# assign new cell type identities based on marker gene expression
monje <- RenameIdents(monje, 
                      "0"="Macrophages/Microglia",
                      "1"="Macrophages/Microglia",
                      "2"="Neurons",
                      "3"="Macrophages/Microglia",
                      "4"="Pericytes",
                      "5"="Endothelial",
                      "6"="Neurons",
                      "7"="Epithelial",
                      "8"="Epithelial",
                      "9"="Oligodendrocytes",
                      "10"="Macrophages/Microglia",
                      "11"="Astrocytes",
                      "12"="Endothelial",
                      "13"="vSMCs",
                      "14"="Macrophages/Microglia",
                      "15"="Macrophages/Microglia",
                      "16"="pvFibroblasts",  
                      "17"="Endothelial",
                      "18"="Neural progenitors",
                      "19"="T lymphocytes",
                      "20"="Neural progenitors",
                      "21"="Astrocytes",
                      "22"="Neurons",
                      "23"="pvFibroblasts",
                      "24"="Ependymal",
                      "25"="B lymphocytes",
                      "26"="Oligodendrocytes",
                      "27"="Neutrophils",
                      "28"="OPCs")

# order labels for plotting
Idents(monje) <- factor(Idents(monje), levels = c("Macrophages/Microglia",
                                                  "Neurons",
                                                  "Neural progenitors",
                                                  "Astrocytes",
                                                  "Oligodendrocytes",
                                                  "OPCs",
                                                  "Endothelial",
                                                  "Epithelial",
                                                  "Pericytes",
                                                  "vSMCs",
                                                  "pvFibroblasts",
                                                  "Ependymal",
                                                  "T lymphocytes",
                                                  "B lymphocytes",
                                                  "Neutrophils"))

# save labels in another metadata slot
monje$final_cell_type_labels <- Idents(monje)

# get color palette for plotting
discrete_rainbow <- khroma::colour("discrete rainbow")
khroma::plot_scheme(discrete_rainbow(15), colours = TRUE, size = 0.7)

# plot initial cell type umap
tiff("intial_clustering_umap_final.tiff", height = 6, width = 6, res = 300, units = "in")
DimPlot(monje) + scale_color_manual(values = as.vector(discrete_rainbow(15))) + 
  theme(aspect.ratio = 1, axis.title.x = element_blank(), axis.ticks.x = element_blank(), axis.line.x = element_blank(),
        axis.title.y = element_blank(), axis.ticks.y = element_blank(), axis.line.y = element_blank(),
        axis.text.x = element_blank(), axis.text.y = element_blank())
dev.off()

postscript("intial_clustering_umap_final_legend.ps", height = 6, width = 6)
as_ggplot(get_legend(DimPlot(monje) + scale_color_manual(values = as.vector(discrete_rainbow(15))) + 
                       theme(aspect.ratio = 1, axis.title.x = element_blank(), axis.ticks.x = element_blank(), axis.line.x = element_blank(),
                             axis.title.y = element_blank(), axis.ticks.y = element_blank(), axis.line.y = element_blank(),
                             axis.text.x = element_blank(), axis.text.y = element_blank())))
dev.off()

# reverse levels for maintaining order in dotplot
Idents(monje) <- factor(Idents(monje), levels = rev(c("Macrophages/Microglia",
                                                      "Neurons",
                                                      "Neural progenitors",
                                                      "Astrocytes",
                                                      "Oligodendrocytes",
                                                      "OPCs",
                                                      "Endothelial",
                                                      "Epithelial",
                                                      "Pericytes",
                                                      "vSMCs",
                                                      "pvFibroblasts",
                                                      "Ependymal",
                                                      "T lymphocytes",
                                                      "B lymphocytes",
                                                      "Neutrophils")))

# plot dot plot with cell type marker expression
postscript("dimplot_intial_clustering_markers_final.ps", height = 6, width = 8)
DotPlot(monje, features = c("Cx3cr1", "Meg3", "Sox11", "Slc1a3", "Mbp", "Cspg4", "Flt1", "Ecrg4", 
                            "Vtn", "Myh11", "Col1a1", "Tmem212", "Trbc2", "Igkc", "S100a9"), scale = TRUE) + 
  scale_color_distiller(palette = "Blues", direction = 0) + labs(y = NULL, x = NULL, size = "% Expressing") + 
  theme(axis.text.x = element_text(angle = 45, hjust = 1)) + guides(color = guide_colorbar(title = 'Avg. Scaled\nExpression'),
                                                                    size = guide_legend(title = '% Expressing'))
dev.off()

# reorder identities 
Idents(monje) <- factor(Idents(monje), levels = c("Macrophages/Microglia",
                                                  "Neurons",
                                                  "Neural progenitors",
                                                  "Astrocytes",
                                                  "Oligodendrocytes",
                                                  "OPCs",
                                                  "Endothelial",
                                                  "Epithelial",
                                                  "Pericytes",
                                                  "vSMCs",
                                                  "pvFibroblasts",
                                                  "Ependymal",
                                                  "T lymphocytes",
                                                  "B lymphocytes",
                                                  "Neutrophils"))

# save initial clustering seurat object
saveRDS(monje, "monje_full_seurat_object_with_final_cell_type_clusters.rds")

# subset microglia/macrophages for subclustering
mgs <- subset(monje, idents = "Macrophages/Microglia")
mgs # 12,643 cells

# split object by sample identities
mgs.list <- SplitObject(mgs, split.by = "orig.ident")
# find new variable features
features <- SelectIntegrationFeatures(mgs.list, nfeatures = 3000)
# scale data, regressing technical variables
mgs <- ScaleData(mgs, features = features, vars.to.regress = c("nFeature_RNA", "nCount_RNA", "percent.mt"))
# run PCA
mgs <- RunPCA(mgs, features = features)

# integrate between samples using Harmony
set.seed(12)
mgs <- mgs %>% RunHarmony(group.by.vars = "orig.ident", plot_convergence = TRUE)

# plot pre- and post- harmony PC elbow plots
ElbowPlot(mgs, ndims = 50, reduction = "pca")
ElbowPlot(mgs, ndims = 50, reduction = "harmony")

# run UMAP using top 20 harmony-corrected principle components
mgs <- RunUMAP(mgs, dims = 1:20, reduction = "harmony")

# plot UMAPs colored by condition and sample
tiff("mg_subclustering_umap_by_type.tiff", units = "in", width = 6, height = 6, res = 300)
DimPlot(mgs, group.by = "type", cols = c("lightblue", "#8b0000")) + theme(aspect.ratio = 1) + ggtitle("Condition")
dev.off()

tiff("mg_subclustering_umap_by_sample.tiff", units = "in", width = 6, height = 6, res = 300)
DimPlot(mgs, group.by = "orig.ident") + theme(aspect.ratio = 1) + ggsci::scale_color_lancet() + ggtitle("Sample")
dev.off()

# perform clustering
mgs <- FindNeighbors(mgs, dims = 1:20, reduction = "harmony")
mgs <- FindClusters(mgs, resolution = seq(0, 0.5, by = 0.1))

# choosing res = 0.1; plotting UMAP with cluster labels
DimPlot(mgs, group.by = "RNA_snn_res.0.1") + khroma::scale_color_discreterainbow() # 9 clusters
DimPlot(mgs, group.by = "RNA_snn_res.0.1", split.by = "type") + khroma::scale_color_discreterainbow() 

# assign cluster labels to new metadata slot
mgs$clustering1 <- mgs$RNA_snn_res.0.1

# run differential gene expression testing to classify clusters
Idents(mgs) <- mgs$clustering1
mgs.clustering1.markers <- FindAllMarkers(mgs)
mgs.clustering1.markers <- mgs.clustering1.markers %>% arrange(cluster, desc(avg_log2FC)) %>% filter(p_val_adj < 0.05)
# save results
saveRDS(mgs.clustering1.markers, "microglia_clustering_round1_differential_expression_testing_results.rds")

# rename cluster identities
mgs <- RenameIdents(mgs, '0'="Microglia", '1'="Macrophages", "2"="Microglia", "3"="Microglia", "4"="Microglia", "5"="Monocytes", "6" = "Erythrocytes", 
                    "7" = "Dendritic cells", "8"= "Pericytes")

mgs$final_contaminating_cluster_labels <- Idents(mgs)

# plot UMAP with new cell type labels
tiff("mg_subclustering_celltypes_final_panel.tiff", units = "in", width = 6, height = 6, res = 300)
DimPlot(mgs) + ggsci::scale_color_npg() + theme(aspect.ratio = 1) + 
  theme(aspect.ratio = 1, axis.title.x = element_blank(), axis.ticks.x = element_blank(), axis.line.x = element_blank(),
        axis.title.y = element_blank(), axis.ticks.y = element_blank(), axis.line.y = element_blank(),
        axis.text.x = element_blank(), axis.text.y = element_blank()) + ggtitle(NULL)
dev.off()

# plot dotplot with marker gene expression
postscript("dimplot_intial_clustering1_markers_final.ps", height = 6, width = 6)
DotPlot(mgs, features = c("Tmem119", "P2ry12", "Mrc1", "F13a1", "Plac8", "Cd44", "Hbb-bs", "Alas2", 
                          "Cd209a", "S100a11", "Vtn", "Rgs5"), scale = TRUE) + 
  scale_color_distiller(palette = "Blues", direction = 0) + labs(y = NULL, x = NULL, size = "% Expressing") + 
  theme(axis.text.x = element_text(angle = 45, hjust = 1)) + guides(color = guide_colorbar(title = 'Avg. Scaled\nExpression'),
                                                                    size = guide_legend(title = '% Expressing'))
dev.off()

# subset seurat object to subcluster only microglia
mgs2 <- subset(mgs, ident = "Microglia")

saveRDS(mgs2, "microglia_subclustering_seurat_object.rds")

# read in saved object for subclustering analysis
microglia <- readRDS("microglia_subclustering_seurat_object.rds")

# split object by sample 
microglia.list <- SplitObject(microglia, split.by = "orig.ident")

# find variable features within each sample
for(i in 1:length(microglia.list)){
  microglia.list[[i]] <- FindVariableFeatures(microglia.list[[i]], nfeatures = 2000)
}
# get consensus variable features across samples
features <- SelectIntegrationFeatures(microglia.list, nfeatures = 2000)

# scale data, regressing technical variables
microglia <- ScaleData(microglia, features = rownames(microglia), vars.to.regress = c("nFeature_RNA", "nCount_RNA", "percent.mt"))
# run PCA
microglia <- RunPCA(microglia, features = features)

# integrate across samples within each condition using harmony
set.seed(92)
microglia <- microglia %>% RunHarmony(group.by.vars = "orig.ident", 
                                      reference_values = c("monje5", "monje6", "monje7", "monje8"), 
                                      plot_convergence = TRUE, reduction.save = 'harmony1') %>% 
  RunHarmony(group.by.vars = "orig.ident", 
             reference_values = c("monje1", "monje2", "monje3", "monje4"), 
             plot_convergence = TRUE, reduction.save = 'harmony2')

# view PC elbow plots pre- and post- harmony correction
ElbowPlot(microglia, ndims = 50, reduction = "pca")
ElbowPlot(microglia, ndims = 50, reduction = "harmony2")

# run UMAP
microglia <- RunUMAP(microglia, dims = 1:20, reduction = "harmony2", reduction.name = "UMAPh2") 

# plot UMAP
DimPlot(microglia, reduction = "UMAPh2", group.by = "orig.ident")
DimPlot(microglia, reduction = "UMAPh2", group.by = "type")

# run clustering
microglia <- FindNeighbors(microglia, dims = 1:20, reduction = "harmony2")
microglia <- FindClusters(microglia, resolution = seq(0, 3, by = 0.05))

# plot cluster tree to identify clusters stable at multiple resolutions 
clustree(microglia, prefix = "RNA_snn_res.")

# choosing k = 2 to parse small clusters that appear stable at high resolution. Outside four small lineages, other clusters appear unstable across most resolutions, and will be merged
Idents(microglia) <- microglia$RNA_snn_res.2

# re-assign identifiers based on cluster tree and marker gene expression (cluster names changed after DE testing in following lines)
# sInflam cluster corresponds to cell cycle/inflammatory label in paper 
microglia <- RenameIdents(microglia, "0"="Homeostatic",
                          "1"="Homeostatic",
                          "2"="Homeostatic",
                          "3"="Homeostatic",
                          "4"="Homeostatic",
                          "5"="Homeostatic",
                          "6"="Homeostatic",
                          "7"="Homeostatic",
                          "8"="Homeostatic",
                          "9"="Homeostatic",
                          "10"="Homeostatic",
                          "11"="Homeostatic",
                          "12"="Homeostatic",
                          "13"="Homeostatic",
                          "14"="Homeostatic",
                          "15"="Homeostatic",
                          "16"= "ATM-like",
                          "17"="Interferon-responsive",
                          "18"="Chemokine",
                          "19"="sInflam") 

# re-order identities
Idents(microglia) <- factor(Idents(microglia), levels = c("Homeostatic", "ATM-like", "Chemokine", "sInflam", "Interferon-responsive"))
microglia$subcluster_labels <- Idents(microglia)

# plot UMAP with final clusters
DimPlot(microglia, reduction = "UMAPh2")

# run cluster differential expression testing
subcluster.markers <- FindAllMarkers(microglia)
subcluster.markers <- subcluster.markers %>% filter(p_val_adj < 0.05) %>% arrange(cluster, desc(avg_log2FC))

# save results
saveRDS(subcluster.markers, "final_significant_microglia_subcluster_de_results.rds")

# get DE results for all genes across each cluster, for volcano plots
subcluster.de.results.all <- FindAllMarkers(microglia, logfc.threshold = 0, min.pct = 0, min.diff.pct = 0, min.cells.feature = 0, min.cells.group = 0)
subcluster.de.results.all <- subcluster.de.results.all %>% arrange(cluster, desc(avg_log2FC))

# save results
saveRDS(subcluster.de.results.all, "final_microglia_subcluster_de_results_allgenes.rds")

# differential abundance testing
propeller.results <- propeller(clusters = Idents(microglia), sample = microglia$orig.ident, group = microglia$type)

# save differential abundance test results
saveRDS(propeller.results, "propeller_results.rds")

# get cell counts per cluster per sample
table(microglia$orig.ident, Idents(microglia))

# make data frame of proportions of microglia per sample in each cluster
subcluster.counts <- as.data.frame.matrix(table(Idents(microglia), microglia$orig.ident))
subcluster.proportions <- apply(subcluster.counts,2,function(x){x/sum(x)})
prop.df.all <- reshape2::melt(subcluster.proportions)
colnames(prop.df.all) <- c("cluster", "sample", "proportion")
prop.df.all <- prop.df.all %>% mutate(condition = ifelse(sample %in% c("monje1", "monje2", "monje3", "monje4"), "Control", "COVID"))
prop.df.all$cluster <- as.factor(prop.df.all$cluster)
prop.df.all$cluster <- factor(prop.df.all$cluster, levels = c("Interferon-responsive", "sInflam", "Chemokine", "ATM-like", "Homeostatic"))

# calculate median proportion fold change between covid and control for each cluster
prop.medians <- prop.df.all %>% group_by(cluster, condition) %>% summarise(median = median(proportion))
tmp.df <- data.frame(cluster = unique(prop.medians$cluster), log2fc = log2((prop.medians %>% filter(condition == "COVID") %>% .$median) / (prop.medians %>% filter(condition == "Control") %>% .$median)))

# plot subcluster proportions plot
cairo_ps(file = "microglia_subtype_proportions.eps", onefile = FALSE, fallback_resolution = 600, height = 4, width = 6)
ggplot(prop.df.all, aes(y=cluster, x=proportion, color=condition)) +
  geom_jitter(aes(shape = condition, color = condition), 
              position = position_jitterdodge(jitter.width = 0.1, dodge.width = 0.8),
              size = 2, alpha = 0.5) + 
  scale_color_manual(values = c("lightblue", "#8b0000")) + 
  scale_x_break(c(0.06, 0.85), space=.5) + 
  scale_x_continuous(labels = scales::percent) + 
  stat_summary(fun.y = median, fun.ymin = median, fun.ymax = median,
               geom = "crossbar", width = 0.8, position = position_dodge(0.8), show.legend = FALSE) + 
  ggpubr::theme_pubr() + labs(color = NULL, x = "Proportion of microglia per sample", y = NULL) + 
  theme(axis.text.y = element_blank()) + NoLegend()
dev.off()

cairo_ps(file = "microglia_subtype_proportions_legend.eps", onefile = FALSE, fallback_resolution = 600, height = 4, width = 6)
as_ggplot(get_legend(ggplot(prop.df.all, aes(y=cluster, x=proportion, color=condition)) +
                       geom_jitter(aes(shape = condition, color = condition), 
                                   position = position_jitterdodge(jitter.width = 0.1, dodge.width = 0.8),
                                   size = 2, alpha = 0.5) + 
                       scale_color_manual(values = c("lightblue", "#8b0000")) + 
                       #scale_x_break(c(0.053, 0.865)) + 
                       scale_x_break(c(0.06, 0.85)) + 
                       scale_x_continuous(labels = scales::percent) + 
                       stat_summary(fun.y = median, fun.ymin = median, fun.ymax = median,
                                    geom = "crossbar", width = 0.8, position = position_dodge(0.8), show.legend = FALSE) + 
                       ggpubr::theme_pubr() + labs(color = NULL, x = "Proportion of microglia per sample", y = NULL) + 
                       theme(axis.text.y = element_blank())))
dev.off()

# plot median fold change plot
cairo_ps("test_microglia_subcluster_proportion_foldchange_plot.eps", height = 4, width = 2, onefile = FALSE, fallback_resolution = 600)
ggplot(tmp.df, aes(cluster, log2fc, fill = cluster))+geom_bar(stat="identity") + 
  scale_fill_manual(values = rev(pal[1:5])) +   
  theme_pubr() + coord_flip() + scale_y_reverse() + NoLegend() + theme(axis.title.y = element_blank(),
                                                                       axis.text.y = element_blank(),
                                                                       axis.ticks.y = element_blank(),
                                                                       axis.line.y = element_blank()) + 
  labs(y = "Log2 median fold change")
dev.off()

# save proportions and medians dataframes
saveRDS(prop.df.all, "cluster_proportions_dataframe.rds")
saveRDS(prop.medians, "cluster_median_proportions_dataframe.rds")

# plot subcluster marker gene violin plot
VlnPlot(microglia, features = c("P2ry12", "P2ry13", "Tmem119", "Apoe", "Lpl", 
                                "Cd63", "Ccl4", "Ccl3", "Atf3", "Gpr84", "Icam1", 
                                "Vcam1", "Ifit3", "Ifitm3", "Ifit2", "Isg15"), 
        stack = TRUE, flip = TRUE) + NoLegend() + theme(axis.text.x = element_blank()) + 
  labs(x = NULL, y = "Log-normalized expression")

# plot UMAPs overalaid with cluster marker genes
tiff("ccl4_featureplot.tiff", units = "in", height = 4, width = 4, res = 300)
FeaturePlot(microglia, reduction = "UMAPh2", features = "Ccl4", pt.size = 0.25, order = TRUE) + scale_color_distiller(palette = "Blues", direction = 0) + labs(x = "UMAP1", y = "UMAP2") + 
  theme(aspect.ratio = 1, 
        axis.title.x = element_blank(), axis.ticks.x = element_blank(), axis.line.x = element_blank(),
        axis.title.y = element_blank(), axis.ticks.y = element_blank(), axis.line.y = element_blank(),
        axis.text.x = element_blank(), axis.text.y = element_blank())
dev.off()

tiff("ccl3_featureplot.tiff", units = "in", height = 4, width = 4, res = 300)
FeaturePlot(microglia, reduction = "UMAPh2", features = "Ccl3", pt.size = 0.25, order = TRUE) + scale_color_distiller(palette = "Blues", direction = 0) + labs(x = "UMAP1", y = "UMAP2") + 
  theme(aspect.ratio = 1, 
        axis.title.x = element_blank(), axis.ticks.x = element_blank(), axis.line.x = element_blank(),
        axis.title.y = element_blank(), axis.ticks.y = element_blank(), axis.line.y = element_blank(),
        axis.text.x = element_blank(), axis.text.y = element_blank())
dev.off()

tiff("ccl12_featureplot.tiff", units = "in", height = 4, width = 4, res = 300)
FeaturePlot(microglia, reduction = "UMAPh2", features = "Ccl12", pt.size = 0.25, order = TRUE) + scale_color_distiller(palette = "Blues", direction = 0) + labs(x = "UMAP1", y = "UMAP2") + 
  theme(aspect.ratio = 1, 
        axis.title.x = element_blank(), axis.ticks.x = element_blank(), axis.line.x = element_blank(),
        axis.title.y = element_blank(), axis.ticks.y = element_blank(), axis.line.y = element_blank(),
        axis.text.x = element_blank(), axis.text.y = element_blank())
dev.off()

tiff("cxcl10_featureplot.tiff", units = "in", height = 4, width = 4, res = 300)
FeaturePlot(microglia, reduction = "UMAPh2", features = "Cxcl10", pt.size = 0.25, order = TRUE) + scale_color_distiller(palette = "Blues", direction = 0) + labs(x = "UMAP1", y = "UMAP2") + 
  theme(aspect.ratio = 1, 
        axis.title.x = element_blank(), axis.ticks.x = element_blank(), axis.line.x = element_blank(),
        axis.title.y = element_blank(), axis.ticks.y = element_blank(), axis.line.y = element_blank(),
        axis.text.x = element_blank(), axis.text.y = element_blank())
dev.off()

tiff("ccl9_featureplot.tiff", units = "in", height = 4, width = 4, res = 300)
FeaturePlot(microglia, reduction = "UMAPh2", features = "Ccl9", pt.size = 0.25, order = TRUE) + scale_color_distiller(palette = "Blues", direction = 0) + labs(x = "UMAP1", y = "UMAP2") + 
  theme(aspect.ratio = 1, 
        axis.title.x = element_blank(), axis.ticks.x = element_blank(), axis.line.x = element_blank(),
        axis.title.y = element_blank(), axis.ticks.y = element_blank(), axis.line.y = element_blank(),
        axis.text.x = element_blank(), axis.text.y = element_blank())
dev.off()


tiff("tnf_featureplot.tiff", units = "in", height = 4, width = 4, res = 300)
FeaturePlot(microglia, reduction = "UMAPh2", features = "Tnf", pt.size = 0.25, order = TRUE) + scale_color_distiller(palette = "Blues", direction = 0) + labs(x = "UMAP1", y = "UMAP2") + 
  theme(aspect.ratio = 1, 
        axis.title.x = element_blank(), axis.ticks.x = element_blank(), axis.line.x = element_blank(),
        axis.title.y = element_blank(), axis.ticks.y = element_blank(), axis.line.y = element_blank(),
        axis.text.x = element_blank(), axis.text.y = element_blank())
dev.off()

tiff("il1a_featureplot.tiff", units = "in", height = 4, width = 4, res = 300)
FeaturePlot(microglia, reduction = "UMAPh2", features = "Il1a", pt.size = 0.25, order = TRUE) + scale_color_distiller(palette = "Blues", direction = 0) + labs(x = "UMAP1", y = "UMAP2") + 
  theme(aspect.ratio = 1, 
        axis.title.x = element_blank(), axis.ticks.x = element_blank(), axis.line.x = element_blank(),
        axis.title.y = element_blank(), axis.ticks.y = element_blank(), axis.line.y = element_blank(),
        axis.text.x = element_blank(), axis.text.y = element_blank())
dev.off()

tiff("il1b_featureplot.tiff", units = "in", height = 4, width = 4, res = 300)
FeaturePlot(microglia, reduction = "UMAPh2", features = "Il1b", pt.size = 0.25, order = TRUE) + scale_color_distiller(palette = "Blues", direction = 0) + labs(x = "UMAP1", y = "UMAP2") + 
  theme(aspect.ratio = 1, 
        axis.title.x = element_blank(), axis.ticks.x = element_blank(), axis.line.x = element_blank(),
        axis.title.y = element_blank(), axis.ticks.y = element_blank(), axis.line.y = element_blank(),
        axis.text.x = element_blank(), axis.text.y = element_blank())
dev.off()

# get color palette
pal <- MetBrewer::met.brewer("Greek",n=5)

# plot subcluster umap 
tiff("subcluster_umap.tiff", units = "in", height = 4, width = 4, res = 300)
DimPlot(microglia, reduction = "UMAPh2", cols = pal) + 
  theme(aspect.ratio = 1, 
        axis.title.x = element_blank(), axis.ticks.x = element_blank(), axis.line.x = element_blank(),
        axis.title.y = element_blank(), axis.ticks.y = element_blank(), axis.line.y = element_blank(),
        axis.text.x = element_blank(), axis.text.y = element_blank()) + NoLegend() + 
  labs(x = "UMAP1", y = "UMAP2") + NoLegend()
dev.off()

tiff("subcluster_umap_legend.tiff", units = "in", height = 4, width = 4, res = 300)
as_ggplot(get_legend(DimPlot(microglia, reduction = "UMAPh2", cols = pal) + theme(aspect.ratio = 1) + labs(x = "UMAP1", y = "UMAP2")))
dev.off()

# plot subcluster umap split by condition
tiff("subcluster_umap_split_by_type.tiff", units = "in", height = 4, width = 8, res = 300)
DimPlot(microglia, reduction = "UMAPh2", cols = pal, split.by = "type") + labs(x = "UMAP1", y = "UMAP2") + 
  theme(aspect.ratio = 1, 
        axis.title.x = element_blank(), axis.ticks.x = element_blank(), axis.line.x = element_blank(),
        axis.title.y = element_blank(), axis.ticks.y = element_blank(), axis.line.y = element_blank(),
        axis.text.x = element_blank(), axis.text.y = element_blank()) + NoLegend()
dev.off()

# create subcluster DE heatmap using top 20 genes per cluster
top.20 <- subcluster.markers %>%
  group_by(cluster) %>%
  top_n(n = 20, wt = avg_log2FC) %>% .$gene

avg.exp.mat <- as.data.frame(AverageExpression(microglia, features = top.20))
colnames(avg.exp.mat) <- c("Homeostatic", "ATM-like", "Chemokine", "sInflam", "Interferon_responsive")
z.avg.exp.mat <- t(scale(t(log1p(avg.exp.mat))))

ha <- HeatmapAnnotation(Cluster = c("Homeostatic", "ATM-like", "Chemokine", "sInflam", "Interferon_responsive"), 
                        col = list(Cluster = c("Homeostatic"=pal[1], "ATM-like"=pal[2], "Chemokine"=pal[3], "sInflam"=pal[4], "Interferon_responsive"=pal[5])),
                        show_legend = FALSE)

p <- Heatmap(z.avg.exp.mat, name = "z-scored\navg. exp.", 
             col = circlize::colorRamp2(breaks = c(-2, -1, 0, 1, 2), colors = rev(RColorBrewer::brewer.pal(name = "RdBu", n = 5))),
             cluster_rows = FALSE,
             cluster_columns = FALSE,
             top_annotation = ha,
             show_column_names = FALSE,
             row_names_gp = gpar(fontface = "italic")
)

# save plot
tiff("subcluster_markers_heatmap.tiff", units = "in", height = 12, width = 8, res = 300)
p
dev.off()

# differential expression testing between chemokine and all other microglia
chemo.de.results <- FindMarkers(microglia, ident.1 = "Chemokine", logfc.threshold = 0, min.pct = 0, min.diff.pct = 0, min.cells.feature = 0, min.cells.group = 0)
chemo.de.results <- chemo.de.results %>% arrange(desc(avg_log2FC))
sig.chemo.de.results <- chemo.de.results %>% filter(abs(avg_log2FC) > 0.25) %>% filter(p_val_adj < 0.05)

# save results
saveRDS(chemo.de.results, "chemo_de_results_all.rds")
saveRDS(sig.chemo.de.results, "sig_chemo_de_results.rds")

# create volcano plot of chemokine DE results
keyvals <- ifelse(
  (chemo.de.results$avg_log2FC < -0.25) & (chemo.de.results$p_val_adj < 0.05), 'royalblue',
  ifelse((chemo.de.results$avg_log2FC > 0.25) & (chemo.de.results$p_val_adj < 0.05), "#8b0000",
         'gray'))
keyvals[is.na(keyvals)] <- 'gray'
names(keyvals)[keyvals == "#8b0000"] <- 'Up'
names(keyvals)[keyvals == 'gray'] <- 'NS'
names(keyvals)[keyvals == 'royalblue'] <- 'Down'

lab_italics <- paste0("italic('", rownames(chemo.de.results), "')")
selectLab_italics = paste0(
  "italic('",
  c('Ccl4', "Ccl3", "Egr1", "Egr2", "Egr3", "Nes", "Ccl2", "Gadd45b", "Fos", "Atf3", "Gas2l3", "Plaur", "Stk38l", "Gem", "Plekho2", "Tnf", "Fam20c", "Pdgfa", "Csf1", "Cxcl10", "Cd83",
    "Id2", "Slc15a3", "Ccl12", "Il1b", "Il1a", "Sall1", "Tmem119", "P2ry12", "Siglech"),
  "')")

# create plot
tiff("chemo_de_volcano_plot.tiff", res=300, height = 8, width = 10, units = "in")
EnhancedVolcano(chemo.de.results,
                lab = lab_italics,
                x = 'avg_log2FC',
                y = 'p_val_adj',
                caption = NULL,
                selectLab = selectLab_italics,
                xlab = bquote(~Log[2]~ 'fold change'),
                pCutoff = 0.05,
                FCcutoff = 0.25,
                pointSize = 3.0,
                labSize = 6.0,
                labCol = 'black',
                labFace = 'bold',
                parseLabels = TRUE,
                colAlpha = 4/5,
                colCustom = keyvals,
                legendPosition = 'bottom',
                legendLabSize = 14,
                legendIconSize = 4.0,
                drawConnectors = TRUE,
                widthConnectors = 1.0,
                colConnectors = 'black',
                xlim = c(-1, 4.1), 
                gridlines.major = FALSE,
                gridlines.minor = FALSE,
                title = NULL, subtitle = NULL, ylab = "-Log10 p-value") + theme(legend.position = "top")
dev.off()

# save final subclustered microglia seurat object
saveRDS(microglia, "microglia_seurat_object__FINAL.rds")
# save integration features
saveRDS(features, "integration_features.rds")

# DE testing between COVID and control within each microglia subcluster
covid.de.testing.df <- data.frame()
for(i in unique(Idents(microglia))){
  print(i)
  tmp.df <- FindMarkers(microglia, ident.1 = "COVID", ident.2 = "control", group.by = "type", subset.ident = i)
  tmp.df <- tmp.df %>% filter(p_val_adj < 0.05) %>% arrange(desc(avg_log2FC))
  tmp.df$Cluster <- i
  tmp.df$Gene <- rownames(tmp.df)
  covid.de.testing.df <- rbind(covid.de.testing.df, tmp.df)
}
covid.de.testing.df

# save results
saveRDS(covid.de.testing.df, "COVID_DE_testing_results_microglia_subclusters.rds")

# subset homeostatic COVID DE results
homeo.covid.de.results <- covid.de.testing.df %>% filter(Cluster == "Homeostatic")

# homeostatic COVID v control DE testing, all genes (for volcano plot)
homeo.covid.response.all.genes <- FindMarkers(microglia, ident.1 = "COVID", ident.2 = "control", 
                                              group.by = "type", subset.ident = "Homeostatic", 
                                              logfc.threshold = 0, min.pct = 0, min.diff.pct = 0, 
                                              min.cells.feature = 0, min.cells.group = 0)

homeo.covid.response.all.genes <- homeo.covid.response.all.genes %>% arrange(desc(avg_log2FC))


keyvals <- ifelse(
  (homeo.covid.response.all.genes$avg_log2FC < -0.25) & (homeo.covid.response.all.genes$p_val_adj < 0.05), 'royalblue',
  ifelse((homeo.covid.response.all.genes$avg_log2FC > 0.25) & (homeo.covid.response.all.genes$p_val_adj < 0.05), "#8b0000",
         'gray'))
keyvals[is.na(keyvals)] <- 'gray'
names(keyvals)[keyvals == "#8b0000"] <- 'Up'
names(keyvals)[keyvals == 'gray'] <- 'NS'
names(keyvals)[keyvals == 'royalblue'] <- 'Down'

lab_italics <- paste0("italic('", rownames(homeo.covid.response.all.genes), "')")

# homeo covid response volcano plot
tiff("homeo_covid_reponse_de_volcano_plot.tiff", res=300, height = 6, width = 10, units = "in")
EnhancedVolcano(homeo.covid.response.all.genes,
                lab = lab_italics,
                x = 'avg_log2FC',
                y = 'p_val_adj',
                caption = NULL,
                xlab = bquote(~Log[2]~ 'fold change'),
                pCutoff = 0.05,
                FCcutoff = 0.25,
                pointSize = 3.0,
                labSize = 6.0,
                labCol = 'black',
                labFace = 'bold',
                parseLabels = TRUE,
                colAlpha = 4/5,
                colCustom = keyvals,
                legendPosition = 'bottom',
                legendLabSize = 14,
                legendIconSize = 4.0,
                drawConnectors = TRUE,
                widthConnectors = 1.0,
                colConnectors = 'black',
                xlim = c(-0.6, 1.1), 
                gridlines.major = FALSE,
                gridlines.minor = FALSE,
                title = NULL, subtitle = NULL, ylab = "-Log10 p-value") + theme(legend.position = "top")
dev.off()

# plot Cd300lf split UMAP overlaid with expression
tiff("cd300lf_featureplot_split_by_condition.tiff", units = "in", height = 4, width = 8, res = 300)
FeaturePlot(microglia, reduction = "UMAPh2", features = "Cd300lf", split.by = "type", pt.size = 0.25, order = TRUE) & scale_color_distiller(palette = "Blues", direction = 0, lim = c(0, 3.5)) & 
  labs(x = "UMAP1", y = "UMAP2") & 
  theme(aspect.ratio = 1, 
        legend.position = c(1, 0.25),
        axis.title.x = element_blank(), axis.ticks.x = element_blank(), axis.line.x = element_blank(),
        axis.title.y = element_blank(), axis.ticks.y = element_blank(), axis.line.y = element_blank(),
        axis.text.x = element_blank(), axis.text.y = element_blank()) & NoLegend()
dev.off()

tiff("cd300lf_featureplot_split_by_condition_legend.tiff", units = "in", height = 4, width = 4, res = 300)
as_ggplot(get_legend(FeaturePlot(microglia, reduction = "UMAPh2", features = "Cd300lf", pt.size = 0.25, order = TRUE) + scale_color_distiller(palette = "Blues", direction = 0, lim = c(0, 3.5))))
dev.off()

# subset out homeostatic cluster
homeostatic <- subset(microglia, idents = "Homeostatic")

# create stacked violin plot of a few COVID DE genes 
postscript("violinplot_homeostatic_covid_markers.ps", height = 4, width = 16)
VlnPlot(homeostatic, features = c("Cd300lf", "Sfi1", "Sla", "H2-D1", "H2-K1", "B2m", "Apoe", "Ubb"), group.by = "type", stack = TRUE) + NoLegend() + labs(x = "Log-Normalized Expression", y = NULL)
dev.off()

# GO analysis using clusterProfiler
# chemokine cluster GO enrichment analysis

# upregulated genes
sig.chemo.de.results
sig.chemo.enterez.ids = bitr(sig.chemo.de.results %>% filter(avg_log2FC > 0) %>% row.names(.), fromType="SYMBOL", toType="ENTREZID", OrgDb="org.Mm.eg.db")

chemo.go.results <- enrichGO(gene= sig.chemo.enterez.ids$ENTREZID,
                             OrgDb= org.Mm.eg.db,
                             ont= "BP",
                             pAdjustMethod = "BH",
                             pvalueCutoff  = 0.01,
                             qvalueCutoff  = 0.05,
                             readable= TRUE)
head(chemo.go.results)
chemo.go.df <- chemo.go.results@result

# downregulated genes
sig.down.chemo.enterez.ids = bitr(sig.chemo.de.results %>% filter(avg_log2FC < 0) %>% row.names(.), fromType="SYMBOL", toType="ENTREZID", OrgDb="org.Mm.eg.db")

down.chemo.go.results <- enrichGO(gene= sig.down.chemo.enterez.ids$ENTREZID,
                                  OrgDb= org.Mm.eg.db,
                                  ont= "BP",
                                  pAdjustMethod = "BH",
                                  pvalueCutoff  = 0.01,
                                  qvalueCutoff  = 0.05,
                                  readable= TRUE)
head(down.chemo.go.results)
down.chemo.go.df <- down.chemo.go.results@result

# save results
saveRDS(sig.chemo.enterez.ids, "chemokine_cluster_upregulated_DEGs_with_GO_analysis_enterezids.rds")
saveRDS(chemo.go.results, "chemokine_cluster_upregulated_DEGs_GO_analysis_results.rds")
saveRDS(chemo.go.df, "chemokine_cluster_upregulated_DEGs_GO_analysis_results_df.rds")
saveRDS(down.chemo.go.results, "chemokine_cluster_downregulated_DEGs_GO_analysis_results.rds")
saveRDS(down.chemo.go.df, "chemokine_cluster_downregulated_DEGs_GO_analysis_results_df.rds")

# GO analysis of homeostatic DEGs, covid v control
# upregulated genes
sig.covid.enterez.ids = bitr(homeo.covid.de.results %>% filter(avg_log2FC > 0) %>% .$Gene, fromType="SYMBOL", toType="ENTREZID", OrgDb="org.Mm.eg.db")


covid.go.results <- enrichGO(gene= sig.covid.enterez.ids$ENTREZID,
                             OrgDb= org.Mm.eg.db,
                             ont= "BP",
                             pAdjustMethod = "BH",
                             pvalueCutoff  = 0.01,
                             qvalueCutoff  = 0.05,
                             readable= TRUE)
head(covid.go.results)
covid.go.df <- covid.go.results@result

# downregulated genes
sig.down.covid.enterez.ids = bitr(homeo.covid.de.results %>% filter(avg_log2FC < 0) %>% .$Gene, fromType="SYMBOL", toType="ENTREZID", OrgDb="org.Mm.eg.db")

down.covid.go.results <- enrichGO(gene= sig.down.covid.enterez.ids$ENTREZID,
                                  OrgDb= org.Mm.eg.db,
                                  ont= "BP",
                                  pAdjustMethod = "BH",
                                  pvalueCutoff  = 0.01,
                                  qvalueCutoff  = 0.05,
                                  readable= TRUE)
head(down.covid.go.results)
down.covid.go.df <- down.covid.go.results@result

# save results
saveRDS(sig.covid.enterez.ids, "homeostatic_cluster_upregulated_DEGs_in_COVID_v_control_with_GO_analysis_enterezids.rds")
saveRDS(covid.go.results, "homeostatic_cluster_upregulated_DEGs_in_COVID_v_control_GO_analysis_results.rds")
saveRDS(covid.go.df, "homeostatic_cluster_upregulated_DEGs_in_COVID_v_control_GO_analysis_results_df.rds")
saveRDS(down.covid.go.results, "homeostatic_cluster_downregulated_DEGs_in_COVID_v_control_GO_analysis_results.rds")
saveRDS(down.covid.go.df, "homeostatic_cluster_downregulated_DEGs_in_COVID_v_control_GO_analysis_results_df.rds")

# from go results, subset out significant results
chemo.up.final <- chemo.go.df %>% filter(p.adjust < 0.01) 
chemo.down.final <- down.chemo.go.df %>% filter(p.adjust < 0.01) 

covid.up.final <- covid.go.df %>% filter(p.adjust < 0.01) 
covid.down.final <- down.covid.go.df %>% filter(p.adjust < 0.01)

# save significant results
saveRDS(chemo.up.final, "chemokine_cluster_significantly_upregulated_GO_terms.rds")
saveRDS(chemo.down.final, "chemokine_cluster_significantly_downregulated_GO_terms.rds")

saveRDS(covid.up.final, "homeostatic_cluster_significantly_upregulated_GO_terms_in_COVID_v_control.rds")
saveRDS(covid.down.final, "homeostatic_cluster_significantly_downregulated_GO_terms_in_COVID_v_control.rds")


# semantic clustering of significant GO terms enriched in Chemokine and COVID signatures (since many GO terms repeat or describe overlapping categories)
# we'll use this package to cluster the GO terms, and then we'll manually annotate each cluster by scanning the GO terms falling in each cluster ourselves

### Chemokine cluster
# get semantic similarity matrix for GO Biological Process terms significantly upregulated in our Chemokine group
chemo.up.mat = GO_similarity(chemo.up.final$ID, ont = "BP")

# print similarity matrix (with cluster word clouds) and create data frame of GO terms labeled with cluster numbers;
# we'll remove the word clouds in Illustrator and write our own descriptions
tiff("chemo_up_simplified_go_terms_similarity_matrix.tiff", units = "in", height = 6, width = 6, res = 300)
simplified.chemo.up.df = simplifyGO(chemo.up.mat) 
dev.off()
# note: last few clusters have only a few terms in them: only annotating the top 9 clusters manually 

## Homeostatic COVID signature 

# get semantic similarity matrix for GO Biological Process terms significantly upregulated in our Homeostatic group in COVID v control
covid.up.mat = GO_similarity(covid.up.final$ID, ont = "BP")

# print similarity matrix (with cluster word clouds) and create data frame of GO terms labeled with cluster numbers;
# we'll remove the word clouds in Illustrator and write our own descriptions
tiff("COVID_up_simplified_go_terms_similarity_matrix.tiff", units = "in", height = 6, width = 6, res = 300)
simplified.covid.up.df = simplifyGO(covid.up.mat) 
dev.off()

# save GO semantic similarity results
saveRDS(chemo.up.mat, "chemokine_upregulated_go_terms_similarity_matrix.rds")
saveRDS(simplified.chemo.up.df, "chemokine_upregulated_go_term_cluster_df.rds")

saveRDS(covid.up.mat, "homeostatic_COVID_v_control_upregulated_go_terms_similarity_matrix.rds")
saveRDS(simplified.covid.up.df, "homeostatic_COVID_v_control_upregulated_go_term_cluster_df.rds")

# create QC plots to evaluate each sample and cluster 
Idents(microglia) <- microglia$orig.ident
microglia <- RenameIdents(microglia, "monje1"="Control1",
                          "monje2"="Control2",
                          "monje3"="Control3",
                          "monje4"="Control4",
                          "monje5"="COVID1",
                          "monje6"="COVID2",
                          "monje7"="COVID3",
                          "monje8"="COVID4")

microglia$sample <- Idents(microglia)

Idents(microglia) <- microglia$subcluster_labels

# plot sample QC plots - UMAPs overlaid with QC metrics
tiff("microglia_umap_by_sample.tiff", units = "in", width = 6, height = 6, res = 300)
DimPlot(microglia, group.by = "sample", reduction = "UMAPh2", pt.size = 0.25) + theme(aspect.ratio = 1) + ggsci::scale_color_lancet() + ggtitle("Sample") + 
  theme(aspect.ratio = 1, 
        axis.title.x = element_blank(), axis.ticks.x = element_blank(), axis.line.x = element_blank(),
        axis.title.y = element_blank(), axis.ticks.y = element_blank(), axis.line.y = element_blank(),
        axis.text.x = element_blank(), axis.text.y = element_blank())
dev.off()

tiff("genes_detected_featureplot.tiff", units = "in", height = 6, width = 6, res = 300)
FeaturePlot(microglia, reduction = "UMAPh2", features = "nFeature_RNA", pt.size = 0.25, order = TRUE) + scale_color_distiller(palette = "Blues", direction = 0) + labs(x = "UMAP1", y = "UMAP2") + 
  ggtitle("Number of genes detected") + 
  theme(aspect.ratio = 1, 
        axis.title.x = element_blank(), axis.ticks.x = element_blank(), axis.line.x = element_blank(),
        axis.title.y = element_blank(), axis.ticks.y = element_blank(), axis.line.y = element_blank(),
        axis.text.x = element_blank(), axis.text.y = element_blank())
dev.off()

tiff("umis_featureplot.tiff", units = "in", height = 6, width = 6, res = 300)
FeaturePlot(microglia, reduction = "UMAPh2", features = "nCount_RNA", pt.size = 0.25, order = TRUE) + scale_color_distiller(palette = "Blues", direction = 0) + labs(x = "UMAP1", y = "UMAP2") + 
  ggtitle("Number of UMIs") + 
  theme(aspect.ratio = 1, 
        axis.title.x = element_blank(), axis.ticks.x = element_blank(), axis.line.x = element_blank(),
        axis.title.y = element_blank(), axis.ticks.y = element_blank(), axis.line.y = element_blank(),
        axis.text.x = element_blank(), axis.text.y = element_blank())
dev.off()

tiff("percent_mt_featureplot.tiff", units = "in", height = 6, width = 6, res = 300)
FeaturePlot(microglia, reduction = "UMAPh2", features = "percent.mt", pt.size = 0.25, order = TRUE) + scale_color_distiller(palette = "Blues", direction = 0) + labs(x = "UMAP1", y = "UMAP2") + 
  ggtitle("% UMIs from mitochondrial genome") + 
  theme(aspect.ratio = 1, 
        axis.title.x = element_blank(), axis.ticks.x = element_blank(), axis.line.x = element_blank(),
        axis.title.y = element_blank(), axis.ticks.y = element_blank(), axis.line.y = element_blank(),
        axis.text.x = element_blank(), axis.text.y = element_blank())
dev.off()

# plot QC metrics in stacked violin plot
postscript("vln_qc_plots_microglia_sample.ps", height = 6, width = 6)
VlnPlot(microglia, group.by = "sample", features = c("nFeature_RNA", "nCount_RNA", "percent.mt"), stack = TRUE, flip = TRUE) + ggsci::scale_color_lancet() + NoLegend() + labs(x = NULL, y = "Log-Normalized Expression")
dev.off()

# plot microglia subcluster QC plots
postscript("vln_qc_plots_microglia_cluster.ps", height = 6, width = 6)
VlnPlot(microglia, group.by = "subcluster_labels", features = c("nFeature_RNA", "nCount_RNA", "percent.mt"), stack = TRUE, flip = TRUE) + ggsci::scale_color_lancet() + NoLegend() + labs(x = NULL, y = "Log-Normalized Expression")
dev.off()


postscript("vln_microglia_subcluster_markers.ps", height = 20, width = 4)
VlnPlot(microglia, group.by = "subcluster_labels", features = rev(c("P2ry12", "P2ry13", "Tmem119", "Apoe", "Cd63", "Ctsb", "Lpl", "Ccl4", "Ccl3", "Atf3", "Egr1", 
                                                                    "Tlr2", "Tnfaip2", "Mcm3", "Mcm6", "Ifit3", "Ifit2", "Stat1", "Isg15")), stack = TRUE, flip = TRUE) + NoLegend() + labs(x = "Cluster", y = "Log-Normalized Expression") 
dev.off()

# Comparing gene signatures between Mild COVID microglia and previous publications

# extract Chemokin / Homeostatic COVID signatures
chemo.up <- sig.chemo.de.results %>% filter(avg_log2FC > 0.25) %>% filter(p_val_adj < 0.05) %>% row.names(.)
chemo.down <- sig.chemo.de.results %>% filter(avg_log2FC < -0.25) %>% filter(p_val_adj < 0.05) %>% row.names(.)

covid.up <- homeo.covid.de.results %>% filter(avg_log2FC > 0.25) %>% filter(p_val_adj < 0.05) %>% .$Gene
covid.down <- homeo.covid.de.results %>% filter(avg_log2FC < -0.25) %>% filter(p_val_adj < 0.05) %>% .$Gene

# load DAM signature, from Keren-Shaul et al., 2017: Supplementary Table 3 
dam.expression.df <- readxl::read_excel("signature_comparisons/dam/dam_differential_expression.xlsx", na = "NaN")
colnames(dam.expression.df) <- c("Gene", "avg_UMIs_per_homeostatic_mg", "fold_change", "neg_log10_pvalue", "FDR")

# extracting significant differentially expressed genes
dam.expression.df.no.na <- na.omit(dam.expression.df)
dam.up.df <- dam.expression.df %>% filter(FDR > -log10(0.05)) %>% filter(fold_change > 0.25) 
dam.up <- dam.expression.df %>% filter(FDR > -log10(0.05)) %>% filter(fold_change > 0.25) %>% .$Gene
dam.down.df <- dam.expression.df %>% filter(FDR > -log10(0.05)) %>% filter(fold_change < -0.25) 
dam.down <- dam.expression.df %>% filter(FDR > -log10(0.05)) %>% filter(fold_change < -0.25) %>% .$Gene

# load WAM seurat object from Safaiyan et al., 2021 (obtained from authors by request)
wam.10x <- readRDS("signature_comparisons/wam/wam_10x.rds")

# run differential expression testing to identify WAM cluster markers versus other microglia
Idents(wam.10x) <- wam.10x$Population 
wam.10x.markers <- FindAllMarkers(wam.10x)
sig.wam.10x.markers <- wam.10x.markers %>% filter(p_val_adj < 0.05) %>% arrange(cluster, desc(avg_log2FC)) %>% filter(cluster == "WAM")

# extract significant DE genes 
wam.up.df <- sig.wam.10x.markers %>% filter(avg_log2FC > 0.25) 
wam.up <- sig.wam.10x.markers %>% filter(avg_log2FC > 0.25) %>% .$gene
wam.down.df <- sig.wam.10x.markers %>% filter(avg_log2FC < -0.25) 
wam.down <- sig.wam.10x.markers %>% filter(avg_log2FC < -0.25) %>% .$gene

# load LPC count matrices from Hammond et al. 2019
data.files <- list.files("signature_comparisons/hammond/counts/")
setwd("signature_comparisons/hammond/counts")
lpc.data <- lapply(data.files, read.table, header = TRUE, row.names = 1, sep = "\t")
setwd("../..")
names(lpc.data) <- stringr::str_replace(data.files, pattern = ".dge.txt.gz", replacement = "")
lpc.names.1 <- names(lpc.data)
names(lpc.data) <- stringr::str_replace(lpc.names.1, pattern = "GSM......._", replacement = "")
sample_ids <- names(lpc.data)
condition <- gsub("P100_M_", "", gsub("_..$","", gsub("_...$","",sample_ids)))

# create list of seurat objects
lpc.data.s <- lapply(lpc.data, CreateSeuratObject, project = "Hammond")

# add metadata
for (i in seq_along(lpc.data.s)){
  lpc.data.s[[i]]@meta.data[["condition"]] <- condition[i]
  lpc.data.s[[i]]@meta.data[["sample_id"]] <- sample_ids[i]
}

# merge seurat objects into one object
lpc <- merge(lpc.data.s[[1]], c(lpc.data.s[[2]], lpc.data.s[[3]], lpc.data.s[[4]], lpc.data.s[[5]], lpc.data.s[[6]]))

# remove low genes detected cells (following authors thresholds used in Hammond et al. 2019)
lpc <- subset(lpc, subset = nFeature_RNA > 650)

# normalize data
lpc <- NormalizeData(lpc) 

# run differential expression testing between all LPC and saline microglia
Idents(lpc) <- lpc$condition
lpc.markers <- FindMarkers(lpc, ident.1 = "LPC", ident.2 = "SALINE")
lpc.markers <- lpc.markers %>% filter(p_val_adj < 0.05) %>% arrange(desc(avg_log2FC))

# extract significant DEGs
lpc.up.df <- lpc.markers %>% filter(avg_log2FC > 0.25)
lpc.up <- rownames(lpc.up.df)
lpc.down.df <- lpc.markers %>% filter(avg_log2FC < -0.25)
lpc.down <- rownames(lpc.down.df)

# load LDAM signature from Marschallinger et al., 2020 Supplemental Table 3
ldam.expression.df <- readxl::read_excel("signature_comparisons/ldam/41593_2019_566_MOESM2_ESM.xlsx", sheet = 4)
colnames(ldam.expression.df) <- c("Gene", "baseMean", "logFC", "p_value", "p_val_adj", "score")

# extract significant DEGs
ldam.expression.df.no.na <- na.omit(ldam.expression.df)
ldam.up.df <- ldam.expression.df %>% filter(p_val_adj < 0.05) %>% filter(logFC > 0.25) 
ldam.up <- ldam.expression.df %>% filter(p_val_adj < 0.05) %>% filter(logFC > 0.25)  %>% .$Gene

ldam.down.df <- ldam.expression.df %>% filter(p_val_adj < 0.05) %>% filter(logFC < -0.25) 
ldam.down <- ldam.expression.df %>% filter(p_val_adj < 0.05) %>% filter(logFC < -0.25)  %>% .$Gene

## Explore overlap between Chemokine/Homeostatic-COVID and previous gene signatures
up.genes <- list(Chemokine = chemo.up,
                 COVID = covid.up,
                 WAM = wam.up,
                 DAM = dam.up,
                 LDAM = ldam.up,
                 LPC = lpc.up
)

# plot Upset plots demonstrating overlap between gene sets -- Upregulated genes
postscript("upset_upregulated_genes.ps", height = 6, width = 8)
upset(fromList(up.genes), nsets = 6, 
      queries = list(list(query = intersects, params = list("Chemokine", "WAM"), color = pal[3], active = T, query.name = "Chemokine Enriched "), 
                     list(query = intersects, params = list("Chemokine", "DAM"), color = pal[3], active = T, query.name = "Chemokine Enriched "),
                     list(query = intersects, params = list("Chemokine", "LPC"), color = pal[3], active = T, query.name = "Chemokine Enriched "),
                     list(query = intersects, params = list("Chemokine", "DAM", "WAM"), color = pal[3], active = T, query.name = "Chemokine Enriched "),
                     list(query = intersects, params = list("Chemokine", "LPC", "WAM"), color = pal[3], active = T, query.name = "Chemokine Enriched "),
                     list(query = intersects, params = list("Chemokine", "LPC", "DAM"), color = pal[3], active = T, query.name = "Chemokine Enriched "),
                     list(query = intersects, params = list("Chemokine", "LDAM", "WAM"), color = pal[3], active = T, query.name = "Chemokine Enriched "),
                     list(query = intersects, params = list("Chemokine", "LPC", "DAM", "WAM"), color = pal[3], active = T, query.name = "Chemokine Enriched "),
                     list(query = intersects, params = list("Chemokine", "LDAM", "DAM", "WAM"), color = pal[3], active = T, query.name = "Chemokine Enriched "),
                     list(query = intersects, params = list("Chemokine", "LDAM", "LPC", "DAM", "WAM"), color = pal[3], active = T, query.name = "Chemokine Enriched "),
                     list(query = intersects, params = list("COVID", "WAM"), color = "#8b0000", active = T, query.name = "Homeostatic Enriched"),
                     list(query = intersects, params = list("COVID", "DAM", "WAM"), color = "#8b0000", active = T, query.name = "Homeostatic Enriched"),
                     list(query = intersects, params = list("COVID", "DAM"), color = "#8b0000", active = T, query.name = "Homeostatic Enriched"),
                     list(query = intersects, params = list("COVID", "LPC", "DAM", "WAM"), color = "#8b0000", query.name = "Homeostatic Enriched"),
                     
                     list(query = intersects, params = list("COVID", "Chemokine"), color = "gold", active = T, query.name = "Chemokine & Homeostatic Enriched"),
                     list(query = intersects, params = list("COVID", "Chemokine", "WAM"), color = "gold", active = T, query.name = "Chemokine & Homeostatic Enriched")
                     
      ),
      query.legend = "none",
      shade.alpha = 0,
      matrix.dot.alpha = 0
      
)
dev.off()

# plot legend
postscript("upset_example_for_legend.ps", height = 6, width = 8)
upset(fromList(up.genes), nsets = 6, 
      queries = list(list(query = intersects, params = list("Chemokine", "WAM"), color = pal[3], active = T, query.name = "Chemokine Enriched "), 
                     list(query = intersects, params = list("COVID", "WAM"), color = "#8b0000", active = T, query.name = "Homeostatic Enriched"),
                     list(query = intersects, params = list("COVID", "Chemokine", "WAM"), color = "gold", active = T, query.name = "Chemokine & Homeostatic Enriched")
                     
      ),
      query.legend = "top",
      shade.alpha = 0,
      matrix.dot.alpha = 0
      
)
dev.off()

# plot Upset plots demonstrating overlap between gene sets -- Downregulated genes
down.genes <- list(Chemokine = chemo.down,
                   COVID = covid.down,
                   WAM = wam.down,
                   DAM = dam.down,
                   LDAM = ldam.down,
                   LPC = lpc.down
)

postscript("upset_downregulated_genes.ps", height = 6, width = 8)
upset(fromList(down.genes), nsets = 6, 
      queries = list(list(query = intersects, params = list("Chemokine"), color = pal[3], active = T, query.name = "Chemokine Enriched "),
                     list(query = intersects, params = list("Chemokine", "WAM"), color = pal[3], active = T, query.name = "Chemokine Enriched "),
                     list(query = intersects, params = list("Chemokine", "DAM", "WAM"), color = pal[3], active = T, query.name = "Chemokine Enriched "),
                     list(query = intersects, params = list("Chemokine", "LPC", "WAM"), color = pal[3], active = T, query.name = "Chemokine Enriched "),
                     list(query = intersects, params = list("Chemokine", "LPC", "DAM", "WAM"), color = pal[3], active = T, query.name = "Chemokine Enriched "),
                     list(query = intersects, params = list("Chemokine", "LDAM", "LPC", "DAM", "WAM"), color = pal[3], active = T, query.name = "Chemokine Enriched "),
                     
                     list(query = intersects, params = list("COVID"), color = "#8b0000", active = T, query.name = "Homeostatic Enriched"),             
                     list(query = intersects, params = list("COVID", "WAM"), color = "#8b0000", active = T, query.name = "Homeostatic Enriched"),
                     list(query = intersects, params = list("COVID", "DAM"), color = "#8b0000", active = T, query.name = "Homeostatic Enriched"),
                     list(query = intersects, params = list("COVID", "LDAM"), color = "#8b0000", active = T, query.name = "Homeostatic Enriched"),
                     list(query = intersects, params = list("COVID", "LPC"), color = "#8b0000", active = T, query.name = "Homeostatic Enriched"),
                     
                     list(query = intersects, params = list("COVID", "DAM", "WAM"), color = "#8b0000", active = T, query.name = "Homeostatic Enriched"),
                     list(query = intersects, params = list("COVID", "LDAM", "DAM"), color = "#8b0000", active = T, query.name = "Homeostatic Enriched"),
                     list(query = intersects, params = list("COVID", "LPC", "WAM"), color = "#8b0000", query.name = "Homeostatic Enriched"),
                     list(query = intersects, params = list("COVID", "LPC", "DAM", "WAM"), color = "#8b0000", query.name = "Homeostatic Enriched"),
                     list(query = intersects, params = list("COVID", "Chemokine", "LPC", "DAM", "WAM"), color = "gold", active = T, query.name = "Chemokine & Homeostatic Enriched")
                     
      ),
      query.legend = "none",
      shade.alpha = 0,
      matrix.dot.alpha = 0
      
)
dev.off()

# plot dot plots quantifying proportion of Chemo or COVID DEGs upregulated or downregulated in WAM, DAM, and LPC datasets
up.overlap.df <- data.frame(set = c(rep("Chemokine", 4), rep("COVID", 4)), overlap_set = c("WAM", "DAM", "LPC", "LDAM", "WAM", "DAM", "LPC", "LDAM"), 
                            proportion = c(length(intersect(chemo.up, wam.up))/length(chemo.up),
                                           length(intersect(chemo.up, dam.up))/length(chemo.up),
                                           length(intersect(chemo.up, lpc.up))/length(chemo.up),
                                           length(intersect(chemo.up, ldam.up))/length(chemo.up),
                                           length(intersect(covid.up, wam.up))/length(covid.up),
                                           length(intersect(covid.up, dam.up))/length(covid.up),
                                           length(intersect(covid.up, lpc.up))/length(covid.up),
                                           length(intersect(covid.up, ldam.up))/length(covid.up)
                            )
)

down.overlap.df <- data.frame(set = c(rep("Chemokine", 4), rep("COVID", 4)), overlap_set = c("WAM", "DAM", "LPC", "LDAM", "WAM", "DAM", "LPC", "LDAM"), 
                              proportion = c(length(intersect(chemo.down, wam.down))/length(chemo.down),
                                             length(intersect(chemo.down, dam.down))/length(chemo.down),
                                             length(intersect(chemo.down, lpc.down))/length(chemo.down),
                                             length(intersect(chemo.down, ldam.down))/length(chemo.down),
                                             length(intersect(covid.down, wam.down))/length(covid.down),
                                             length(intersect(covid.down, dam.down))/length(covid.down),
                                             length(intersect(covid.down, lpc.down))/length(covid.down),
                                             length(intersect(covid.down, ldam.down))/length(covid.down)
                              )
)

up.overlap.df$overlap_set <- factor(up.overlap.df$overlap_set, levels = c("WAM", "DAM", "LPC", "LDAM"))
up.overlap.df$set <- factor(up.overlap.df$set, levels = c("COVID", "Chemokine"))

down.overlap.df$overlap_set <- factor(down.overlap.df$overlap_set, levels = c("WAM", "DAM", "LPC", "LDAM"))
down.overlap.df$set <- factor(down.overlap.df$set, levels = c("COVID", "Chemokine"))

postscript("overlapping_deg_dotplot_upregulated_genes.ps", height = 4, width = 6)
ggcorrplot(reshape2::acast(up.overlap.df, overlap_set~set, value.var = "proportion"), method = "circle", outline.color = "white") + 
  scale_fill_distiller(palette = "Blues", direction = 0, lim = c(0, max(up.overlap.df$proportion + 0.01)), labels = scales::percent) +
  theme_pubr() + theme(legend.position = "right") + labs(y = NULL, x = "", fill = "% overlapping\nDEGs")
dev.off()

postscript("overlapping_deg_dotplot_downregulated_genes.ps", height = 4, width = 6)
ggcorrplot(reshape2::acast(down.overlap.df, overlap_set~set, value.var = "proportion"), method = "circle", outline.color = "white") + 
  scale_fill_distiller(palette = "Blues", direction = 0, lim = c(0, max(down.overlap.df$proportion + 0.01)), labels = scales::percent) +
  theme_pubr() + theme(legend.position = "right") + labs(y = NULL, x = "", fill = "% overlapping\nDEGs")
dev.off()

# create data frame with overlap statistics - upregulated genes
up.overlap.complete.df <- data.frame(set = c(rep("Chemokine", 4), rep("COVID", 4)), overlap_set = c("WAM", "DAM", "LPC", "LDAM", "WAM", "DAM", "LPC", "LDAM"), 
                                     set_count = c(length(chemo.up),
                                                   length(chemo.up),
                                                   length(chemo.up),
                                                   length(chemo.up),
                                                   length(covid.up),
                                                   length(covid.up),
                                                   length(covid.up),
                                                   length(covid.up)
                                     ),
                                     overlap_set_count = c(length(wam.up),
                                                           length(dam.up),
                                                           length(lpc.up),
                                                           length(ldam.up),
                                                           length(wam.up),
                                                           length(dam.up),
                                                           length(lpc.up),
                                                           length(ldam.up)
                                     ),
                                     n_overlaps = c(length(intersect(chemo.up, wam.up)),
                                                    length(intersect(chemo.up, dam.up)),
                                                    length(intersect(chemo.up, lpc.up)),
                                                    length(intersect(chemo.up, ldam.up)),
                                                    length(intersect(covid.up, wam.up)),
                                                    length(intersect(covid.up, dam.up)),
                                                    length(intersect(covid.up, lpc.up)),
                                                    length(intersect(covid.up, ldam.up))
                                     ),
                                     proportion = c(length(intersect(chemo.up, wam.up))/length(chemo.up),
                                                    length(intersect(chemo.up, dam.up))/length(chemo.up),
                                                    length(intersect(chemo.up, lpc.up))/length(chemo.up),
                                                    length(intersect(chemo.up, ldam.up))/length(chemo.up),
                                                    length(intersect(covid.up, wam.up))/length(covid.up),
                                                    length(intersect(covid.up, dam.up))/length(covid.up),
                                                    length(intersect(covid.up, lpc.up))/length(covid.up),
                                                    length(intersect(covid.up, ldam.up))/length(covid.up)
                                     )
)

# create data frame with overlap statistics - downregulated genes
down.overlap.complete.df <- data.frame(set = c(rep("Chemokine", 4), rep("COVID", 4)), overlap_set = c("WAM", "DAM", "LPC", "LDAM", "WAM", "DAM", "LPC", "LDAM"), 
                                       set_count = c(length(chemo.down),
                                                     length(chemo.down),
                                                     length(chemo.down),
                                                     length(chemo.down),
                                                     length(covid.down),
                                                     length(covid.down),
                                                     length(covid.down),
                                                     length(covid.down)
                                       ),
                                       overlap_set_count = c(length(wam.down),
                                                             length(dam.down),
                                                             length(lpc.down),
                                                             length(ldam.down),
                                                             length(wam.down),
                                                             length(dam.down),
                                                             length(lpc.down),
                                                             length(ldam.down)
                                       ),
                                       n_overlaps = c(length(intersect(chemo.down, wam.down)),
                                                      length(intersect(chemo.down, dam.down)),
                                                      length(intersect(chemo.down, lpc.down)),
                                                      length(intersect(chemo.down, ldam.down)),
                                                      length(intersect(covid.down, wam.down)),
                                                      length(intersect(covid.down, dam.down)),
                                                      length(intersect(covid.down, lpc.down)),
                                                      length(intersect(covid.down, ldam.down))
                                       ),
                                       proportion = c(length(intersect(chemo.down, wam.down))/length(chemo.down),
                                                      length(intersect(chemo.down, dam.down))/length(chemo.down),
                                                      length(intersect(chemo.down, lpc.down))/length(chemo.down),
                                                      length(intersect(chemo.down, ldam.down))/length(chemo.down),
                                                      length(intersect(covid.down, wam.down))/length(covid.down),
                                                      length(intersect(covid.down, dam.down))/length(covid.down),
                                                      length(intersect(covid.down, lpc.down))/length(covid.down),
                                                      length(intersect(covid.down, ldam.down))/length(covid.down)
                                       )
)

# create data frame listing overlapping genes between each set - upregulated genes
up.overlap.genes.df <- data.frame(set = c(rep("Chemokine", 4), rep("COVID", 4)), overlap_set = c("WAM", "DAM", "LPC", "LDAM", "WAM", "DAM", "LPC", "LDAM"), 
                                  overlap = c(paste0(intersect(chemo.up, wam.up), collapse=", "),
                                              paste0(intersect(chemo.up, dam.up), collapse=", "),
                                              paste0(intersect(chemo.up, lpc.up), collapse=", "),
                                              paste0(intersect(chemo.up, ldam.up), collapse=", "),
                                              paste0(intersect(covid.up, wam.up), collapse=", "),
                                              paste0(intersect(covid.up, dam.up), collapse=", "),
                                              paste0(intersect(covid.up, lpc.up), collapse=", "),
                                              paste0(intersect(covid.up, ldam.up), collapse=", ")
                                  )
)
# create data frame listing overlapping genes between each set - downregulated genes
down.overlap.genes.df <- data.frame(set = c(rep("Chemokine", 4), rep("COVID", 4)), overlap_set = c("WAM", "DAM", "LPC", "LDAM", "WAM", "DAM", "LPC", "LDAM"), 
                                    proportion = c(paste0(intersect(chemo.down, wam.down), collapse=", "),
                                                   paste0(intersect(chemo.down, dam.down), collapse=", "),
                                                   paste0(intersect(chemo.down, lpc.down), collapse=", "),
                                                   paste0(intersect(chemo.down, ldam.down), collapse=", "),
                                                   paste0(intersect(covid.down, wam.down), collapse=", "),
                                                   paste0(intersect(covid.down, dam.down), collapse=", "),
                                                   paste0(intersect(covid.down, lpc.down), collapse=", "),
                                                   paste0(intersect(covid.down, ldam.down), collapse=", ")
                                    )
)

# reformat data frames for creation of supplementary table
sig.chemo.de.results$gene <- rownames(sig.chemo.de.results) 
propeller.results$BaselineProp.clusters <- c("Chemokine", "Cell cycle/Inflammatory", "Homeostatic", "Interferon-responsive", "ATM-like")
prop.medians.mat <- reshape2::acast(prop.medians, formula = cluster~condition, value.var = "median")

# create Supplementary Excel file
sheets <- list("KEY" = c(),
               "Sheet1"=cluster.markers, 
               "Sheet1b"=data.frame(cluster = c(0, 1, 2, 3, 4, 5, 6, 7, 8, 9, 10,
                                                11, 12, 13, 14, 15, 16, 17, 18, 19, 20,
                                                21, 22, 23, 24, 25, 26, 27, 28), 
                                    cell_type = c("Macrophages/Microglia",
                                                  "Macrophages/Microglia",
                                                  "Neurons",
                                                  "Macrophages/Microglia",
                                                  "Pericytes",
                                                  "Endothelial",
                                                  "Neurons",
                                                  "Epithelial",
                                                  "Epithelial",
                                                  "Oligodendrocytes",
                                                  "Macrophages/Microglia",
                                                  "Astrocytes",
                                                  "Endothelial",
                                                  "vSMCs",
                                                  "Macrophages/Microglia",
                                                  "Macrophages/Microglia",
                                                  "pvFibroblasts",  
                                                  "Endothelial",
                                                  "Neural progenitors",
                                                  "T lymphocytes",
                                                  "Neural progenitors",
                                                  "Astrocytes",
                                                  "Neurons",
                                                  "pvFibroblasts",
                                                  "Ependymal",
                                                  "B lymphocytes",
                                                  "Oligodendrocytes",
                                                  "Neutrophils",
                                                  "OPCs")),
               "Sheet2"=mgs.clustering1.markers,
               "Sheet2b"=data.frame(cluster = c(0, 1, 2, 3, 4, 5, 6, 7, 8), 
                                    cell_type = c("Microglia", "Macrophages", "Microglia", "Microglia", "Microglia", "Monocytes", "Erythrocytes", 
                                                  "Dendritic cells", "Pericytes")),
               "Sheet3"=subcluster.markers,
               "Sheet4"=prop.df.all,
               "Sheet5"=propeller.results,
               "Sheet6"=sig.chemo.de.results,
               "Sheet7"=covid.de.testing.df,
               "Sheet8"=chemo.go.df,
               "Sheet9"=down.chemo.go.df,
               "Sheet10"=covid.go.df,
               "Sheet11"=down.covid.go.df,
               "Sheet12"=simplified.chemo.up.df,
               "Sheet13"=simplified.covid.up.df,
               "Sheet14"=up.overlap.complete.df,
               "Sheet15"=down.overlap.complete.df,
               "Sheet16"=up.overlap.genes.df,
               "Sheet17"=down.overlap.genes.df
)

# save Excel file
write.xlsx(sheets, file = 'Supplementary_Table_X.xlsx')

# save data
saveRDS(cluster.markers, "first_round_clustering_differential_expression_test_results.rds")
saveRDS(mgs.clustering1.markers, "second_round_clustering_differential_expression_test_results.rds")
saveRDS(subcluster.markers, "microglia_states_differential_expression_test_results.rds")

saveRDS(prop.df.all, "proportion_of_microglia_in_each_state_per_sample_dataframe.rds")
saveRDS(propeller.results, "differential_abundance_test_COVID_v_control_results.rds")

saveRDS(sig.chemo.de.results, "chemokine_cluster_vs_all_other_microglia_differential_expression_test_results_significant_genes_only.rds")
saveRDS(covid.de.testing.df, "covid_vs_control_in_each_microglia_state_differential_expression_test_results_significant_genes_only.rds")

saveRDS(chemo.go.df, "chemokine_cluster_vs_all_other_microglia_upregulated_gene_ontology_results.rds")
saveRDS(down.chemo.go.df, "chemokine_cluster_vs_all_other_microglia_downregulated_gene_ontology_results.rds")

saveRDS(covid.go.df, "homeostatic_cluster_COVID_vs_Control_upregulated_gene_ontology_results.rds")
saveRDS(down.covid.go.df, "homeostatic_cluster_COVID_vs_Control_downregulated_gene_ontology_results.rds")

saveRDS(simplified.chemo.up.df, "chemokine_cluster_vs_all_other_microglia_upregulated_gene_ontology_results_clustering.rds")
saveRDS(simplified.covid.up.df, "homeostatic_cluster_COVID_vs_Control_upregulated_gene_ontology_results_clustering.rds")

saveRDS(up.overlap.df, "proportion_overlapping_degs_up.rds")
saveRDS(down.overlap.df, "proportion_overlapping_degs_down.rds")

saveRDS(up.overlap.complete.df, "overlapping_upregulated_gene_sets_dataframe.rds")
saveRDS(down.overlap.complete.df, "overlapping_downregulated_gene_sets_dataframe.rds")

saveRDS(up.overlap.genes.df, "overlapping_upregulated_gene_sets_gene_lists_dataframe.rds")
saveRDS(down.overlap.genes.df, "overlapping_downregulated_gene_sets_gene_lists_dataframe.rds")
