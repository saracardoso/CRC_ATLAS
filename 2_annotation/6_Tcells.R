project_dir = '/home/scardoso/Documents/PhD/CRC_ATLAS'
source(paste(project_dir, 'utils/modified_plots.R', sep='/'))

Tcells = SeuratDisk::LoadH5Seurat('/home/scardoso/Documents/PhD/CRC_ATLAS/2_annotation/results_Bcells/datasets/Tcells.h5Seurat')





# -----
# - Find First Clusters
# -----

# 1. Run PCA:
Tcells = Seurat::RunPCA(Tcells, assay='integrated')
invisible(gc())

# 2. Choose number of PCs to use (where the elbow occurs):
pct = Tcells[["pca"]]@stdev /sum(Tcells[["pca"]]@stdev) * 100
cumu = cumsum(pct)
co1 = which(cumu > 90 & pct < 5)[1]
co2 = sort(which((pct[1:length(pct) - 1] - pct[2:length(pct)]) > 0.1), decreasing = T)[1] + 1
# 2.1 Number of PCs to use
elbow = min(co1, co2) # 15
# 2.2. Visual representation of the PCs to use
plot_df = data.frame(pct=pct, cumu=cumu, rank=1:length(pct))
ggplot2::ggplot(plot_df, ggplot2::aes(cumu, pct, label = rank, color = rank > elbow)) + 
  ggplot2::geom_text() + 
  ggplot2::geom_vline(xintercept = 90, color = "grey") + 
  ggplot2::geom_hline(yintercept = min(pct[pct > 5]), color = "grey") +
  ggplot2::theme_bw()
# 2.3. Heatmap of the first 18 PCs
Seurat::DimHeatmap(Tcells, dims=1:18, cells=500, balanced=TRUE)

# 3. Find first clusters
# 3.1. Determine the K-nearest neighbor graph
Tcells = Seurat::FindNeighbors(Tcells, dims=1:elbow)
# 3.2. Find clusters:
Tcells = Seurat::FindClusters(Tcells, resolution=seq(0.1, 1, by=.1))
# 3.3. UMAP:
Tcells = Seurat::RunUMAP(Tcells, dims=1:elbow, reduction='pca')
invisible(gc())

# 4. Choose best resolution:
library(ggraph)
clust_tree = clustree::clustree(Tcells, prefix='integrated_snn_res.')
clust_tree
# 4.1. Visualize different resolutions (first resolution was chosen to be the best for now):
clusters_01 = Seurat::DimPlot(Tcells, reduction="umap", group.by='integrated_snn_res.0.1', label=TRUE, label.size=6, pt.size=.2)
clusters_02 = Seurat::DimPlot(Tcells, reduction="umap", group.by='integrated_snn_res.0.2', label=TRUE, label.size=6, pt.size=.2)
clusters_03 = Seurat::DimPlot(Tcells, reduction="umap", group.by='integrated_snn_res.0.3', label=TRUE, label.size=6, pt.size=.2)
clusters_04 = Seurat::DimPlot(Tcells, reduction="umap", group.by='integrated_snn_res.0.4', label=TRUE, label.size=6, pt.size=.2)
((clusters_01 | clusters_02) / (clusters_03 | clusters_04)) | clust_tree
# 4.2. Plot Metrics
metrics =  c("nUMI", "nGene", "S.Score", "G2M.Score", "percent.mitochondrial_RNA")
feature_plots(Tcells, metrics, ncol=3, with_dimplot=TRUE, dimplot_group='integrated_snn_res.0.2')
# 4.2. Feature Plots:
feature_plots(Tcells, c('rna_CD3E', 'rna_CD3D', 'rna_CD3G', 'rna_CD8A', 'rna_CD4', 'rna_TRDC', 'rna_TRGC1', 'rna_TRGC2'),
              ncol=3, with_dimplot=TRUE, dimplot_group='integrated_snn_res.0.2')
feature_plots(Tcells, c('rna_KLRB1', 'rna_CXCL13', 'rna_IL2RA', 'rna_IL17A'),
              ncol=2, with_dimplot=TRUE, dimplot_group='integrated_snn_res.0.2')
# 4.4. Resolution 0.2 is chosen

# 5. First Annotations:
# Cluster 0     --> CD4 (and some CD8)
# Cluster 1     --> CD8
# Cluster 2     --> Regulatory Tcells
# Cluster 3     --> gd Tcells + NK + ILCs
# Cluster 4     --> CD8 + NKT + CD4, expressing CXCL13
# Cluster 5     --> CD4 (mostly Th17)
# Cluster 6     --> Proliferative Tcells
# Cluster 7     --> CD4

# 6. Save object:
SeuratDisk::SaveH5Seurat(Tcells, paste(project_dir, '2_annotation/results_Bcells/datasets/Tcells_clustered.h5Seurat', sep='/'))



# ---
# - Check SIGMA for further clusterability (res.0.2)
# ---

# 1. Get data processed for SIGMA
# 1.1. Subset cells to half of the original
n_cells = round(dim(Tcells@assays$integrated@data)[2] / 4)
sampled_cells = colnames(Tcells@assays$integrated@data)[sample(1:dim(Tcells@assays$integrated@data)[2], n_cells)]
# 1.2. Check if subsetting still resembles unsubsetted dataset
clusters_02 | Seurat::DimPlot(subset(Tcells, cells=sampled_cells), reduction="umap", group.by='integrated_snn_res.0.2',
                              label=TRUE, label.size=6, pt.size=.2)
# 1.3. Subset dataset
Tcells_subset = subset(Tcells, cells=sampled_cells)
invisible(gc())
# 1.4. Get integrated dataset (because of the integration between datasets)
sigma_data = as.matrix(SeuratObject::GetAssayData(Tcells_subset, assay='integrated', slot='data'))
invisible(gc())
# 1.5. Get vector with resolution 0.1 clusters and dataset
clusters_02_vector = Tcells_subset@meta.data[colnames(sigma_data), 'integrated_snn_res.0.2']
dataset_vector = Tcells_subset@meta.data[colnames(sigma_data), 'dataset']

# 2. Run SIGMA
SIGMA_all_res02 = SIGMA::sigma_funct(sigma_data, clusters_02_vector)
# 2.1. Save SIGMA result
saveRDS(SIGMA_all_res02, paste(project_dir, '2_annotation/results_Bcells/datasets/SIGMA_all_res02.Rdata', sep='/'))

# 3. Check results
# 3.1. Clusterability of all clusters
SIGMA::plot_sigma(SIGMA_all_res02) # All clusters have a sigma value greater than 0.90.
SIGMA_all_res02$maximum_measure # min (Cluster 4) = 0.9424960; max (Cluster 0) = 0.9761377
# 3.2. Check that clusterability is not due to datasets
dtSIGMA_C0 = SIGMA::plot_singular_vectors(SIGMA_all_res02, '0', colour=dataset_vector[Tcells_subset@meta.data$integrated_snn_res.0.2=='0']) +
  ggplot2::ggtitle('Cluster 0')
dtSIGMA_C1 = SIGMA::plot_singular_vectors(SIGMA_all_res02, '1', colour=dataset_vector[Tcells_subset@meta.data$integrated_snn_res.0.2=='1']) +
  ggplot2::ggtitle('Cluster 1')
dtSIGMA_C2 = SIGMA::plot_singular_vectors(SIGMA_all_res02, '2', colour=dataset_vector[Tcells_subset@meta.data$integrated_snn_res.0.2=='2']) +
  ggplot2::ggtitle('Cluster 2')
dtSIGMA_C3 = SIGMA::plot_singular_vectors(SIGMA_all_res02, '3', colour=dataset_vector[Tcells_subset@meta.data$integrated_snn_res.0.2=='3']) +
  ggplot2::ggtitle('Cluster 3')
dtSIGMA_C4 = SIGMA::plot_singular_vectors(SIGMA_all_res02, '4', colour=dataset_vector[Tcells_subset@meta.data$integrated_snn_res.0.2=='4']) +
  ggplot2::ggtitle('Cluster 4')
dtSIGMA_C5 = SIGMA::plot_singular_vectors(SIGMA_all_res02, '5', colour=dataset_vector[Tcells_subset@meta.data$integrated_snn_res.0.2=='5']) +
  ggplot2::ggtitle('Cluster 5')
dtSIGMA_C6 = SIGMA::plot_singular_vectors(SIGMA_all_res02, '6', colour=dataset_vector[Tcells_subset@meta.data$integrated_snn_res.0.2=='6']) +
  ggplot2::ggtitle('Cluster 6')
dtSIGMA_C7 = SIGMA::plot_singular_vectors(SIGMA_all_res02, '7', colour=dataset_vector[Tcells_subset@meta.data$integrated_snn_res.0.2=='7']) +
  ggplot2::ggtitle('Cluster 7')
patchwork::wrap_plots(list(dtSIGMA_C0, dtSIGMA_C1, dtSIGMA_C2, dtSIGMA_C3, dtSIGMA_C4, dtSIGMA_C5, dtSIGMA_C6,dtSIGMA_C7), ncol=2)

# Remove SIGMA objects to continue analysis:
remove(n_cells, sampled_cells, Tcells_subset, sigma_data, clusters_02_vector, dataset_vector, SIGMA_all_res02)
invisible(gc())





# -----
# - Sub-cluster cluster 0 from resolution 0.2 separately
# -----
#Tcells = SeuratDisk::LoadH5Seurat(paste(project_dir, '2_annotation/results_Bcells/datasets/Tcells_clustered.h5Seurat', sep='/'))


# 1. Sub-cluster:
Tcells_0 = subset(Tcells, cells=rownames(Tcells@meta.data)[Tcells$integrated_snn_res.0.2=='0'])
(Seurat::DimPlot(Tcells, group.by = 'integrated_snn_res.0.2', label=T, label.size=4) + ggplot2::theme_minimal() + Seurat::NoLegend()) |
  (Seurat::DimPlot(Tcells_0, group.by = 'integrated_snn_res.0.2', label=T, label.size=4) + ggplot2::theme_minimal() + Seurat::NoLegend())
# 1.1. PCA
Tcells_0 = Seurat::RunPCA(Tcells_0, assay='integrated')
# 1.2. Choose number of PCs to use
pct = Tcells_0[["pca"]]@stdev /sum(Tcells_0[["pca"]]@stdev) * 100
cumu = cumsum(pct)
co1 = which(cumu > 90 & pct < 5)[1]
co2 = sort(which((pct[1:length(pct) - 1] - pct[2:length(pct)]) > 0.1), decreasing = T)[1] + 1
elbow = min(co1, co2) # 5
# 1.3. Visual representation of the PCs to use
plot_df = data.frame(pct=pct, cumu=cumu, rank=1:length(pct))
ggplot2::ggplot(plot_df, ggplot2::aes(cumu, pct, label = rank, color = rank > elbow)) + 
  ggplot2::geom_text() + 
  ggplot2::geom_vline(xintercept = 90, color = "grey") + 
  ggplot2::geom_hline(yintercept = min(pct[pct > 5]), color = "grey") +
  ggplot2::theme_bw()
# 1.4. Heatmap of the first 6 PCs
Seurat::DimHeatmap(Tcells_0, dims=1:6, cells=500, balanced=TRUE)
# 1.5. Determine the K-nearest neighbor graph
Tcells_0 = Seurat::FindNeighbors(Tcells_0, dims=1:elbow)
# 1.6. Find clusters:
Tcells_0 = Seurat::FindClusters(Tcells_0, resolution=seq(0.1, 1, by=.1))
# 1.7. UMAP:
Tcells_0 = Seurat::RunUMAP(Tcells_0, dims=1:elbow, reduction='pca')
invisible(gc())
# 1.8. Save object
SeuratDisk::SaveH5Seurat(Tcells_0, paste(project_dir, '2_annotation/results_Bcells/datasets/Tcells_0.h5Seurat', sep='/'))


# 2. Choose best resolution using clustree:
library(ggraph)
clust_tree = clustree::clustree(Tcells_0, prefix='integrated_snn_res.')
clust_tree
# 2.1. Visualize different resolutions (first resolution was chosen to be the best for now):
clusters_01_0 = Seurat::DimPlot(Tcells_0, reduction="umap", group.by='integrated_snn_res.0.1', label=TRUE, label.size=6, pt.size=.2)
clusters_02_0 = Seurat::DimPlot(Tcells_0, reduction="umap", group.by='integrated_snn_res.0.2', label=TRUE, label.size=6, pt.size=.2)
clusters_03_0 = Seurat::DimPlot(Tcells_0, reduction="umap", group.by='integrated_snn_res.0.3', label=TRUE, label.size=6, pt.size=.2)
clusters_04_0 = Seurat::DimPlot(Tcells_0, reduction="umap", group.by='integrated_snn_res.0.4', label=TRUE, label.size=6, pt.size=.2)
((clusters_01_0 | clusters_02_0) / (clusters_03_0 | clusters_04_0)) | clust_tree
# 2.2. Feature Plots to assess known gene markers:
feature_plots(Tcells_0, c('rna_CCR7', 'rna_SELL', 'rna_PASK'),
              ncol=2, with_dimplot=TRUE, dimplot_group='integrated_snn_res.0.2') # Naive, CM
feature_plots(Tcells_0, c('rna_KLRB1', 'rna_S100A4', 'rna_PTPRCAP', 'rna_PFN1', 'rna_GZMA'),
              ncol=2, with_dimplot=TRUE, dimplot_group='integrated_snn_res.0.2') # EM
feature_plots(Tcells_0, c('rna_STAT4', 'rna_IL12RB2', 'rna_IFNG', 'rna_GATA3', 'rna_STAT6', 'rna_IL4', 'rna_CXCL13', 'rna_IL17A'),
              ncol=3, with_dimplot=TRUE, dimplot_group='integrated_snn_res.0.2')
# 2.3. Violin Plots to assess known markers:
violin_plots(Tcells_0, c('rna_CCR7', 'rna_SELL', 'rna_PASK', 'rna_KLRB1', 'rna_S100A4', 'rna_PTPRCAP', 'rna_PFN1', 'rna_GZMA'),
             ncol=3, with_dimplot=TRUE, group.by='integrated_snn_res.0.2')
# 2.4. Resolution 0.2 will be chosen


# 3. Find markers
# 3.1. Calculate all markers
Seurat::Idents(Tcells_0) = 'integrated_snn_res.0.2'
Tcells_0_markers = Seurat::FindAllMarkers(Tcells_0, assay='RNA', slot='data', logfc.threshold=.5, only.pos=TRUE)
View(Tcells_0_markers)
write.csv(Tcells_0_markers, paste(project_dir, '2_annotation/results_Bcells/markers/markers_C0_res02.csv', sep='/'))
# 3.2. Get top 20 markers for each sub-cluster
Tcells_0_markers_top20 = c()
for(clust in as.character(unique(Tcells_0_markers$cluster))){
  clust_markers = Tcells_0_markers[Tcells_0_markers$cluster==clust,]
  n_markers = dim(clust_markers)[1]
  order_by_avg_log2FC = order(clust_markers$avg_log2FC, decreasing=T)
  if(n_markers > 20) top20_clust_markers = clust_markers[order_by_avg_log2FC[1:20],]
  else top20_clust_markers = clust_markers[order_by_avg_log2FC,]
  Tcells_0_markers_top20 = rbind(Tcells_0_markers_top20, top20_clust_markers)
}
View(Tcells_0_markers_top20)
write.csv(Tcells_0_markers_top20, paste(project_dir, '2_annotation/results_Bcells/markers/markers_C0_res02_top20.csv', sep='/'))


# 4. Annotation:
# Sub-cluster 0        --> Effector Memory (EM) CD4
# Sub-cluster 1        --> Naive CD4
# Sub-cluster 2        --> EM CD4
# Sub-cluster 3        --> Central Memory (CM) CD4


# 5. Check previous SIGMA result and colour cluster 0 by new sub-clusters
SIGMA_all_res02 = readRDS(paste(project_dir, '2_annotation/results_Bcells/datasets/SIGMA_all_res02.Rdata', sep='/'))
cluster_0_cells = rownames(Tcells_0@meta.data)
subclusters_0 = Tcells_0$integrated_snn_res.0.2[intersect(colnames(SIGMA_all_res02$input_parameters$expr), cluster_0_cells)]
subclusters_firstAnnotations = as.character(subclusters_0)
subclusters_firstAnnotations[subclusters_firstAnnotations%in%c('0','2')] = 'EM CD4'
subclusters_firstAnnotations[subclusters_firstAnnotations=='1'] = 'Naive CD4'
subclusters_firstAnnotations[subclusters_firstAnnotations=='3'] = 'CM CD4'
SIGMA::plot_singular_vectors(SIGMA_all_res02, '0', colour=subclusters_firstAnnotations)


# 6. Annotate sub-clusters:
# Sub-cluster 0        --> Effector Memory (EM) CD4
# Sub-cluster 1        --> Naive CD4
# Sub-cluster 2        --> EM CD4
# Sub-cluster 3        --> Central Memory (CM) CD4
# 6.1. Insert the annotations in the metadata:
Tcells_0[['Final_Annotation']] = rep('', length(Tcells_0$integrated_snn_res.0.2))
CM_cells = rownames(Tcells_0@meta.data)[Tcells_0$integrated_snn_res.0.2=='3']
Tcells_0@meta.data[CM_cells, 'Final_Annotation'] = 'CM CD4 Tcells'
naive_cells = rownames(Tcells_0@meta.data)[Tcells_0$integrated_snn_res.0.2=='1']
Tcells_0@meta.data[naive_cells, 'Final_Annotation'] = 'Naive CD4 Tcells'
EM_cells = rownames(Tcells_0@meta.data)[Tcells_0$integrated_snn_res.0.2%in%c('0', '2')]
Tcells_0@meta.data[EM_cells, 'Final_Annotation'] = 'EM CD4 Tcells'
# 7.2. Visualize annotations:
Seurat::DimPlot(Tcells_0, group.by='Final_Annotation', pt.size=0.5, label=T, label.size=4) +
  ggplot2::ggtitle('') + ggplot2::theme_minimal() + Seurat::NoLegend()
# 7.3. Save changes:
SeuratDisk::SaveH5Seurat(Tcells_0, paste(project_dir, '2_annotation/results_Bcells/datasets/Tcells_0.h5Seurat', sep='/'), overwrite=TRUE)





# -----
# - Sub-cluster cluster 1 from resolution 0.2 separately
# -----
#Tcells = SeuratDisk::LoadH5Seurat(paste(project_dir, '2_annotation/results_Bcells/datasets/Tcells_clustered.h5Seurat', sep='/'))


# 1. Sub-cluster:
Tcells_1 = subset(Tcells, cells=rownames(Tcells@meta.data)[Tcells$integrated_snn_res.0.2=='1'])
(Seurat::DimPlot(Tcells, group.by = 'integrated_snn_res.0.2', label=T, label.size=4) + ggplot2::theme_minimal() + Seurat::NoLegend()) |
  (Seurat::DimPlot(Tcells_1, group.by = 'integrated_snn_res.0.2', label=T, label.size=4) + ggplot2::theme_minimal() + Seurat::NoLegend())
# 1.1. PCA
Tcells_1 = Seurat::RunPCA(Tcells_1, assay='integrated')
# 1.2. Choose number of PCs to use
pct = Tcells_1[["pca"]]@stdev /sum(Tcells_1[["pca"]]@stdev) * 100
cumu = cumsum(pct)
co1 = which(cumu > 90 & pct < 5)[1]
co2 = sort(which((pct[1:length(pct) - 1] - pct[2:length(pct)]) > 0.1), decreasing = T)[1] + 1
elbow = min(co1, co2) # 11
# 1.3. Visual representation of the PCs to use
plot_df = data.frame(pct=pct, cumu=cumu, rank=1:length(pct))
ggplot2::ggplot(plot_df, ggplot2::aes(cumu, pct, label = rank, color = rank > elbow)) + 
  ggplot2::geom_text() + 
  ggplot2::geom_vline(xintercept = 90, color = "grey") + 
  ggplot2::geom_hline(yintercept = min(pct[pct > 5]), color = "grey") +
  ggplot2::theme_bw()
# 1.4. Heatmap of the first 12 PCs
Seurat::DimHeatmap(Tcells_1, dims=1:12, cells=500, balanced=TRUE)
# 1.5. Determine the K-nearest neighbor graph
Tcells_1 = Seurat::FindNeighbors(Tcells_1, dims=1:elbow)
# 1.6. Find clusters:
Tcells_1 = Seurat::FindClusters(Tcells_1, resolution=seq(0.1, 1, by=.1))
# 1.7. UMAP:
Tcells_1 = Seurat::RunUMAP(Tcells_1, dims=1:elbow, reduction='pca')
invisible(gc())
# 1.8. Save object
SeuratDisk::SaveH5Seurat(Tcells_1, paste(project_dir, '2_annotation/results_Bcells/datasets/Tcells_1.h5Seurat', sep='/'))


# 2. Choose best resolution using clustree:
library(ggraph)
clust_tree = clustree::clustree(Tcells_1, prefix='integrated_snn_res.')
invisible(gc())
clust_tree
# 2.1. Visualize different resolutions (first resolution was chosen to be the best for now):
clusters_01_1 = Seurat::DimPlot(Tcells_1, reduction="umap", group.by='integrated_snn_res.0.1', label=TRUE, label.size=6, pt.size=.2)
clusters_06_1 = Seurat::DimPlot(Tcells_1, reduction="umap", group.by='integrated_snn_res.0.6', label=TRUE, label.size=6, pt.size=.2)
clusters_08_1 = Seurat::DimPlot(Tcells_1, reduction="umap", group.by='integrated_snn_res.0.8', label=TRUE, label.size=6, pt.size=.2)
clusters_09_1 = Seurat::DimPlot(Tcells_1, reduction="umap", group.by='integrated_snn_res.0.9', label=TRUE, label.size=6, pt.size=.2)
((clusters_01_1 | clusters_06_1) / (clusters_08_1 | clusters_09_1)) | clust_tree
# 2.2. Feature Plots to assess known gene markers:
feature_plots(Tcells_1, c('rna_CD8A', 'rna_CD8B', 'rna_TRAC', 'rna_TRBC1', 'rna_TRBC2',
                          'rna_TRDC', 'rna_TRGC1', 'rna_TRGC2'),
              ncol=3, with_dimplot=TRUE, dimplot_group='integrated_snn_res.0.9', point_size = .5)
feature_plots(Tcells_1, c('rna_CD8A', 'rna_CD8B', 'rna_TRAV1-2', 'rna_IFNG',
                          'rna_LAG3', 'rna_PDCD1', 'rna_TBX21'),
              ncol=3, with_dimplot=TRUE, dimplot_group='integrated_snn_res.0.9', point_size = .5) # CD8aa IELs or classical MAIT
feature_plots(Tcells_1, c('rna_TBX21', 'rna_IFNG', 'rna_STAT4', 'rna_IL12RB2'),
              ncol=3, with_dimplot=TRUE, dimplot_group='integrated_snn_res.0.9', point_size = .5) # Th1 like
feature_plots(Tcells_1, c('rna_GATA3', 'rna_IL5', 'rna_IL13', 'rna_STAT6', 'rna_IL4', 'rna_CSF2'),
              ncol=3, with_dimplot=TRUE, dimplot_group='integrated_snn_res.0.9', point_size = .5) # Th2 like
feature_plots(Tcells_1, c('rna_RORC', 'rna_IL17A', 'rna_IL17B', 'rna_IL17C', 'rna_IL17D', 'rna_IL17F',
                          'rna_IL22', 'rna_CCR6', 'rna_KLRB1', 'rna_RORA'),
              ncol=3, with_dimplot=TRUE, dimplot_group='integrated_snn_res.0.9', point_size = .5) # Th17 like
feature_plots(Tcells_1, c('rna_NCR1', 'rna_NCAM1', 'rna_TYROBP', 'rna_FGFBP2', 'rna_KLRD1', 'rna_KLRF1', 'rna_KLRB1', 'rna_CX3CR1',
                          'rna_FCGR3A', 'rna_XCL1', 'rna_XCL2', 'rna_GNLY', 'rna_NKG7', 'rna_PRF1'),
              ncol=3, with_dimplot=TRUE, dimplot_group='integrated_snn_res.0.9', point_size = .5) # NK like
feature_plots(Tcells_1, c('rna_KLRB1', 'rna_GNLY', 'rna_NKG7', 'rna_GZMA', 'rna_GZMB', 'rna_GZMH', 'rna_GZMK', 'rna_GZMM', 'rna_PRF1',
                          'rna_IFNG'),
              ncol=3, with_dimplot=TRUE, dimplot_group='integrated_snn_res.0.9', point_size = .5) # Citotoxicity
feature_plots(Tcells_1, c('rna_HAVCR2', 'rna_LAG3', 'rna_PDCD1', 'rna_CTLA4', 'rna_TIGIT', 'rna_BTLA', 'rna_KLRC1'),
              ncol=3, with_dimplot=TRUE, dimplot_group='integrated_snn_res.0.9', point_size = .5) # Inhibitory markers
feature_plots(Tcells_1, c('rna_CCR7', 'rna_SELL', 'rna_LEF1', 'rna_TCF7', 'rna_PASK', 'rna_ZNF683',
                          'rna_ANXA1', 'rna_ANKRD28', 'rna_IL7R', 'rna_CD69', 'rna_CD40LG'),
              ncol=3, with_dimplot=TRUE, dimplot_group='integrated_snn_res.0.9', point_size = .5) # naive, CM, TRM, and EM
# 2.4. Resolution 0.9 will be chosen


# 3. Find markers
# 3.1. Calculate all markers
Seurat::Idents(Tcells_1) = 'integrated_snn_res.0.9'
Tcells_1_markers = Seurat::FindAllMarkers(Tcells_1, assay='RNA', slot='data', logfc.threshold=.5, only.pos=TRUE)
View(Tcells_1_markers)
write.csv(Tcells_1_markers, paste(project_dir, '2_annotation/results_Bcells/markers/markers_C1_res09.csv', sep='/'))
# 3.2. Get top 20 markers for each sub-cluster
Tcells_1_markers_top20 = c()
for(clust in as.character(unique(Tcells_1_markers$cluster))){
  clust_markers = Tcells_1_markers[Tcells_1_markers$cluster==clust,]
  n_markers = dim(clust_markers)[1]
  order_by_avg_log2FC = order(clust_markers$avg_log2FC, decreasing=T)
  if(n_markers > 20) top20_clust_markers = clust_markers[order_by_avg_log2FC[1:20],]
  else top20_clust_markers = clust_markers[order_by_avg_log2FC,]
  Tcells_1_markers_top20 = rbind(Tcells_1_markers_top20, top20_clust_markers)
}
View(Tcells_1_markers_top20)
write.csv(Tcells_1_markers_top20, paste(project_dir, '2_annotation/results_Bcells/markers/markers_C1_res09_top20.csv', sep='/'))


# 4. Annotate sub-clusters:
# Sub-cluster 0                         --> CD160+ CD8 Tcells
# Sub-cluster 1                         --> cytotoxic CD8 Tcells
# Sub-cluster 2                         --> DN Tcells
# Sub-cluster 3, 4, 5, 11, 13, 14, 16   --> EM CD8 Tcells
# Sub-cluster 6, 8, 12                  --> CD8aa IELs
# Sub-cluster 7                         --> TRM CD8 Tcells
# Sub-cluster 9                         --> cytotoxic CD8aa IELs
# Sub-cluster 10                        --> DN Tcells? NKT cells?
# Sub-cluster 15                        --> Naive CD8 Tcells
# 4.1. Insert the annotations in the metadata:
Tcells_1[['Final_Annotation']] = rep('', length(Tcells_1$integrated_snn_res.0.9))
cd160_cells = rownames(Tcells_1@meta.data)[Tcells_1$integrated_snn_res.0.9=='0']
Tcells_1@meta.data[cd160_cells, 'Final_Annotation'] = 'CD160+ CD8 Tcells'
CTL_cells = rownames(Tcells_1@meta.data)[Tcells_1$integrated_snn_res.0.9=='1']
Tcells_1@meta.data[CTL_cells, 'Final_Annotation'] = 'cytotoxic CD8 Tcells'
dn_cells = rownames(Tcells_1@meta.data)[Tcells_1$integrated_snn_res.0.9=='2']
Tcells_1@meta.data[dn_cells, 'Final_Annotation'] = 'DN Tcells'
EM_cells = rownames(Tcells_1@meta.data)[Tcells_1$integrated_snn_res.0.9%in%c('3', '4', '5', '11', '13', '14', '16')]
Tcells_1@meta.data[EM_cells, 'Final_Annotation'] = 'EM CD8 Tcells'
cd8aa_cells = rownames(Tcells_1@meta.data)[Tcells_1$integrated_snn_res.0.9%in%c('6', '8', '12')]
Tcells_1@meta.data[cd8aa_cells, 'Final_Annotation'] = 'CD8aa IELs'
TRM_cells = rownames(Tcells_1@meta.data)[Tcells_1$integrated_snn_res.0.9=='7']
Tcells_1@meta.data[TRM_cells, 'Final_Annotation'] = 'TRM CD8 Tcells'
citcd8aa_cells = rownames(Tcells_1@meta.data)[Tcells_1$integrated_snn_res.0.9=='9']
Tcells_1@meta.data[citcd8aa_cells, 'Final_Annotation'] = 'cytotoxic CD8aa IELs'
DNNKT_cells = rownames(Tcells_1@meta.data)[Tcells_1$integrated_snn_res.0.9=='10']
Tcells_1@meta.data[DNNKT_cells, 'Final_Annotation'] = 'DN Tcells? NKT cells?'
naive_cells = rownames(Tcells_1@meta.data)[Tcells_1$integrated_snn_res.0.9=='15']
Tcells_1@meta.data[naive_cells, 'Final_Annotation'] = 'Naive CD8 Tcells'
# 4.2. Check previous SIGMA result and colour cluster 0 by new sub-clusters
SIGMA_all_res02 = readRDS(paste(project_dir, '2_annotation/results_Bcells/datasets/SIGMA_all_res02.Rdata', sep='/'))
cluster_1_cells = rownames(Tcells_1@meta.data)
subclusters_1 = Tcells_1$Final_Annotation[intersect(colnames(SIGMA_all_res02$input_parameters$expr), cluster_1_cells)]
SIGMA::plot_singular_vectors(SIGMA_all_res02, '1', colour=as.character(subclusters_1))
# 4.3. Visualize annotations:
Seurat::DimPlot(Tcells_1, group.by='Final_Annotation', pt.size=0.5, label=F) +
  ggplot2::ggtitle('') + ggplot2::theme_minimal()# + Seurat::NoLegend()
# 4.4. Save changes:
SeuratDisk::SaveH5Seurat(Tcells_1, paste(project_dir, '2_annotation/results_Bcells/datasets/Tcells_1.h5Seurat', sep='/'), overwrite=TRUE)





# -----
# - Sub-cluster cluster 2 from resolution 0.1 separately
# -----
#Tcells = SeuratDisk::LoadH5Seurat(paste(project_dir, '2_annotation/results_Tcells/datasets/Tcells_clustered.h5Seurat', sep='/'))


# 1. Sub-cluster:
Tcells_2 = subset(Tcells, cells=rownames(Tcells@meta.data)[Tcells$integrated_snn_res.0.2=='2'])
(Seurat::DimPlot(Tcells, group.by = 'integrated_snn_res.0.2', label=T, label.size=4) + ggplot2::theme_minimal() + Seurat::NoLegend()) |
  (Seurat::DimPlot(Tcells_2, group.by = 'integrated_snn_res.0.2', label=T, label.size=4) + ggplot2::theme_minimal() + Seurat::NoLegend())
# 1.1. PCA
Tcells_2 = Seurat::RunPCA(Tcells_2, assay='integrated')
# 1.2. Choose number of PCs to use
pct = Tcells_2[["pca"]]@stdev /sum(Tcells_2[["pca"]]@stdev) * 100
cumu = cumsum(pct)
co1 = which(cumu > 90 & pct < 5)[1]
co2 = sort(which((pct[1:length(pct) - 1] - pct[2:length(pct)]) > 0.1), decreasing = T)[1] + 1
elbow = min(co1, co2) # 9
# 1.3. Visual representation of the PCs to use
plot_df = data.frame(pct=pct, cumu=cumu, rank=1:length(pct))
ggplot2::ggplot(plot_df, ggplot2::aes(cumu, pct, label = rank, color = rank > elbow)) + 
  ggplot2::geom_text() + 
  ggplot2::geom_vline(xintercept = 90, color = "grey") + 
  ggplot2::geom_hline(yintercept = min(pct[pct > 5]), color = "grey") +
  ggplot2::theme_bw()
# 1.4. Heatmap of the first 12 PCs
Seurat::DimHeatmap(Tcells_2, dims=1:12, cells=500, balanced=TRUE)
# 1.5. Determine the K-nearest neighbor graph
Tcells_2 = Seurat::FindNeighbors(Tcells_2, dims=1:elbow)
# 1.6. Find clusters:
Tcells_2 = Seurat::FindClusters(Tcells_2, resolution=seq(0.1, 1, by=.1))
# 1.7. UMAP:
Tcells_2 = Seurat::RunUMAP(Tcells_2, dims=1:elbow, reduction='pca')
invisible(gc())
# 1.8. Save object
SeuratDisk::SaveH5Seurat(Tcells_2, paste(project_dir, '2_annotation/results_Tcells/datasets/Tcells_2.h5Seurat', sep='/'), overwrite=TRUE)


# 2. Choose best resolution using clustree:
library(ggraph)
clust_tree = clustree::clustree(Tcells_2, prefix='integrated_snn_res.')
invisible(gc())
clust_tree
# 2.1. Visualize different resolutions (first resolution was chosen to be the best for now):
clusters_01_0 = Seurat::DimPlot(Tcells_2, reduction="umap", group.by='integrated_snn_res.0.1', label=TRUE, label.size=6, pt.size=.5)
clusters_02_0 = Seurat::DimPlot(Tcells_2, reduction="umap", group.by='integrated_snn_res.0.2', label=TRUE, label.size=6, pt.size=.5)
clusters_03_0 = Seurat::DimPlot(Tcells_2, reduction="umap", group.by='integrated_snn_res.0.3', label=TRUE, label.size=6, pt.size=.5)
clusters_05_0 = Seurat::DimPlot(Tcells_2, reduction="umap", group.by='integrated_snn_res.0.5', label=TRUE, label.size=6, pt.size=.5)
((clusters_01_0 | clusters_02_0) / (clusters_03_0 | clusters_05_0)) | clust_tree
# 2.2. Feature Plots to assess known gene markers:
feature_plots(Tcells_2, c('rna_CD3D', 'rna_CD3E', 'rna_CD3G', 'rna_CD8A', 'rna_CD8B', 'rna_CD4', 'rna_TRAC', 'rna_TRBC1', 'rna_TRBC2',
                          'rna_TRDC', 'rna_TRGC1'),
              ncol=3, with_dimplot=TRUE, dimplot_group='integrated_snn_res.0.2', point_size = .5)
feature_plots(Tcells_2, c('rna_CCR7', 'rna_SELL', 'rna_PASK'),
              ncol=2, with_dimplot=TRUE, dimplot_group='integrated_snn_res.0.2') # Naive, CM
feature_plots(Tcells_2, c('rna_KLRB1', 'rna_S100A4', 'rna_PTPRCAP', 'rna_PFN1', 'rna_GZMA'),
              ncol=2, with_dimplot=TRUE, dimplot_group='integrated_snn_res.0.2') # EM
feature_plots(Tcells_2, c('rna_IL2RA', 'rna_IL2RB', 'rna_FOXP3'),
              ncol=2, with_dimplot=TRUE, dimplot_group='integrated_snn_res.0.2') # Regulatory CD4
feature_plots(Tcells_2, c('rna_STAT4', 'rna_IL12RB2', 'rna_IFNG', 'rna_GATA3', 'rna_STAT6', 'rna_IL4', 'rna_CXCL13', 'rna_IL17A'),
              ncol=3, with_dimplot=TRUE, dimplot_group='integrated_snn_res.0.2') # Other Thelper
# 2.3. Violin Plots to assess known markers:
violin_plots(Tcells_2, c('rna_CCR7', 'rna_SELL', 'rna_PASK', 'rna_KLRB1', 'rna_S100A4', 'rna_PTPRCAP', 'rna_PFN1', 'rna_GZMA',
                         'rna_IL2RA', 'rna_IL2RB', 'rna_FOXP3'),
             ncol=3, with_dimplot=TRUE, group.by='integrated_snn_res.0.2')
# 2.4. Resolution 0.2will be chosen


# 3. Find markers
# 3.1. Calculate all markers
Seurat::Idents(Tcells_2) = 'integrated_snn_res.0.2'
Tcells_2_markers = Seurat::FindAllMarkers(Tcells_2, assay='RNA', slot='data', logfc.threshold=.5, only.pos=TRUE)
View(Tcells_2_markers)
write.csv(Tcells_2_markers, paste(project_dir, '2_annotation/results_Tcells/markers/markers_C2_res02.csv', sep='/'))
# 3.2. Get top 20 markers for each sub-cluster
Tcells_2_markers_top20 = c()
for(clust in as.character(unique(Tcells_2_markers$cluster))){
  clust_markers = Tcells_2_markers[Tcells_2_markers$cluster==clust,]
  n_markers = dim(clust_markers)[1]
  order_by_avg_log2FC = order(clust_markers$avg_log2FC, decreasing=T)
  if(n_markers > 20) top20_clust_markers = clust_markers[order_by_avg_log2FC[1:20],]
  else top20_clust_markers = clust_markers[order_by_avg_log2FC,]
  Tcells_2_markers_top20 = rbind(Tcells_2_markers_top20, top20_clust_markers)
}
View(Tcells_2_markers_top20)
write.csv(Tcells_2_markers_top20, paste(project_dir, '2_annotation/results_Tcells/markers/markers_C2_res02_top20.csv', sep='/'))


# 4. Annotate sub-clusters:
Tcells_2[['Final_Annotation']] = rep('Regulatory Tcells', length(Tcells_2$integrated_snn_res.0.2))
# 4.1. Visualize annotations:
Seurat::DimPlot(Tcells_2, group.by='Final_Annotation', pt.size=0.5, label=T, label.size=4) +
  ggplot2::ggtitle('') + ggplot2::theme_minimal() + Seurat::NoLegend()
# 4.2. Save changes:
SeuratDisk::SaveH5Seurat(Tcells_2, paste(project_dir, '2_annotation/results_Tcells/datasets/Tcells_2.h5Seurat', sep='/'), overwrite=TRUE)





# -----
# - Sub-cluster cluster 3 from resolution 0.1 separately
# -----
#Tcells = SeuratDisk::LoadH5Seurat(paste(project_dir, '2_annotation/results_Tcells/datasets/Tcells_clustered.h5Seurat', sep='/'))


# 1. Sub-cluster:
Tcells_3 = subset(Tcells, cells=rownames(Tcells@meta.data)[Tcells$integrated_snn_res.0.2=='3'])
(Seurat::DimPlot(Tcells, group.by = 'integrated_snn_res.0.2', label=T, label.size=4) + ggplot2::theme_minimal() + Seurat::NoLegend()) |
  (Seurat::DimPlot(Tcells_3, group.by = 'integrated_snn_res.0.2', label=T, label.size=4) + ggplot2::theme_minimal() + Seurat::NoLegend())
# 1.1. PCA
Tcells_3 = Seurat::RunPCA(Tcells_3, assay='integrated')
# 1.2. Choose number of PCs to use
pct = Tcells_3[["pca"]]@stdev /sum(Tcells_3[["pca"]]@stdev) * 100
cumu = cumsum(pct)
co1 = which(cumu > 90 & pct < 5)[1]
co2 = sort(which((pct[1:length(pct) - 1] - pct[2:length(pct)]) > 0.1), decreasing = T)[1] + 1
elbow = min(co1, co2) # 11
# 1.3. Visual representation of the PCs to use
plot_df = data.frame(pct=pct, cumu=cumu, rank=1:length(pct))
ggplot2::ggplot(plot_df, ggplot2::aes(cumu, pct, label = rank, color = rank > elbow)) + 
  ggplot2::geom_text() + 
  ggplot2::geom_vline(xintercept = 90, color = "grey") + 
  ggplot2::geom_hline(yintercept = min(pct[pct > 5]), color = "grey") +
  ggplot2::theme_bw()
# 1.4. Heatmap of the first 12 PCs
Seurat::DimHeatmap(Tcells_3, dims=1:12, cells=500, balanced=TRUE)
# 1.5. Determine the K-nearest neighbor graph
Tcells_3 = Seurat::FindNeighbors(Tcells_3, dims=1:elbow)
# 1.6. Find clusters:
Tcells_3 = Seurat::FindClusters(Tcells_3, resolution=seq(0.1, 1, by=.1))
# 1.7. UMAP:
Tcells_3 = Seurat::RunUMAP(Tcells_3, dims=1:elbow, reduction='pca')
invisible(gc())
# 1.8. Save object
SeuratDisk::SaveH5Seurat(Tcells_3, paste(project_dir, '2_annotation/results_Tcells/datasets/Tcells_3.h5Seurat', sep='/'), overwrite=TRUE)


# 2. Choose best resolution using clustree:
library(ggraph)
clust_tree = clustree::clustree(Tcells_3, prefix='integrated_snn_res.')
invisible(gc())
clust_tree
# 2.1. Visualize different resolutions (first resolution was chosen to be the best for now):
clusters_02_0 = Seurat::DimPlot(Tcells_3, reduction="umap", group.by='integrated_snn_res.0.2', label=TRUE, label.size=6, pt.size=.5)
clusters_04_0 = Seurat::DimPlot(Tcells_3, reduction="umap", group.by='integrated_snn_res.0.4', label=TRUE, label.size=6, pt.size=.5)
clusters_05_0 = Seurat::DimPlot(Tcells_3, reduction="umap", group.by='integrated_snn_res.0.5', label=TRUE, label.size=6, pt.size=.5)
clusters_07_0 = Seurat::DimPlot(Tcells_3, reduction="umap", group.by='integrated_snn_res.0.7', label=TRUE, label.size=6, pt.size=.5)
((clusters_02_0 | clusters_04_0) / (clusters_05_0 | clusters_07_0)) | clust_tree
# 2.2. Feature Plots to assess known gene markers:
feature_plots(Tcells_3, c('rna_CD3D', 'rna_CD3E', 'rna_CD3G', 'rna_CD8A', 'rna_CD8B',
                          'rna_TRDC', 'rna_TRGC1', 'rna_TRGC2'),
              ncol=3, with_dimplot=TRUE, dimplot_group='integrated_snn_res.0.7', point_size = .5)
feature_plots(Tcells_3, c('rna_CD8A', 'rna_CD8B', 'rna_TRAC', 'rna_TRBC1', 'rna_TRBC2'),
              ncol=3, with_dimplot=TRUE, dimplot_group='integrated_snn_res.0.7', point_size = .5)
feature_plots(Tcells_3, c('rna_CD8A', 'rna_CD8B', 'rna_TRAV1-2', 'rna_IFNG',
                          'rna_LAG3', 'rna_PDCD1', 'rna_TBX21'),
              ncol=3, with_dimplot=TRUE, dimplot_group='integrated_snn_res.0.7', point_size = .5) # CD8aa IELs or classical MAIT
feature_plots(Tcells_3, c('rna_RORC', 'rna_LTA', 'rna_LTB'),
              ncol=2, with_dimplot=TRUE, dimplot_group='integrated_snn_res.0.7', point_size = .5) # ILC (LTi)
feature_plots(Tcells_3, c('rna_NCR1', 'rna_NCAM1', 'rna_TYROBP', 'rna_FGFBP2', 'rna_KLRD1', 'rna_KLRF1', 'rna_KLRB1', 'rna_CX3CR1',
                          'rna_FCGR3A', 'rna_XCL1', 'rna_XCL2', 'rna_GNLY', 'rna_NKG7', 'rna_PRF1'),
              ncol=3, with_dimplot=TRUE, dimplot_group='integrated_snn_res.0.7', point_size = .5) # NK like
feature_plots(Tcells_3, c('rna_KLRB1', 'rna_GNLY', 'rna_NKG7', 'rna_GZMA', 'rna_GZMB', 'rna_GZMH', 'rna_GZMK', 'rna_GZMM', 'rna_PRF1',
                          'rna_IFNG'),
              ncol=3, with_dimplot=TRUE, dimplot_group='integrated_snn_res.0.7', point_size = .5) # Citotoxicity
feature_plots(Tcells_3, c('rna_HAVCR2', 'rna_LAG3', 'rna_PDCD1', 'rna_CTLA4', 'rna_TIGIT', 'rna_BTLA', 'rna_KLRC1'),
              ncol=3, with_dimplot=TRUE, dimplot_group='integrated_snn_res.0.7', point_size = .5) # Inhibitory markers
feature_plots(Tcells_3, c('rna_HSPA6', 'rna_HSPA1B', 'rna_HSPA1A', 'rna_HSPB1', 'rna_HSP90AA1',
                          'rna_HSPD1', 'rna_HSPH1'),
              ncol=3, with_dimplot=TRUE, dimplot_group='integrated_snn_res.0.7', point_size = .5) # Heat-shock proteins
# 2.3. Violin Plots to assess known markers:
violin_plots(Tcells_3, c('rna_CD3D', 'rna_CD3E', 'rna_CD3G', 'rna_CD8A', 'rna_CD8B',
                         'rna_TRDC', 'rna_TRGC1', 'rna_TRGC2'),
             ncol=3, with_dimplot=TRUE, group.by='integrated_snn_res.0.7')
violin_plots(Tcells_3, c('rna_CD8A', 'rna_CD8B', 'rna_TRAC', 'rna_TRBC1', 'rna_TRBC2'),
             ncol=3, with_dimplot=TRUE, group.by='integrated_snn_res.0.7')
violin_plots(Tcells_3, c('rna_CD8A', 'rna_CD8B', 'rna_TRAV1-2', 'rna_IFNG',
                         'rna_LAG3', 'rna_PDCD1', 'rna_TBX21'),
             ncol=3, with_dimplot=TRUE, group.by='integrated_snn_res.0.7') # CD8aa IELs or classical MAIT
violin_plots(Tcells_3, c('rna_KLRB1', 'rna_GNLY', 'rna_NKG7', 'rna_GZMA', 'rna_GZMB', 'rna_GZMH', 'rna_GZMK', 'rna_GZMM', 'rna_PRF1',
                         'rna_IFNG'),
             ncol=3, with_dimplot=TRUE, group.by='integrated_snn_res.0.7') # Citotoxicity
violin_plots(Tcells_3, c('rna_NCR1', 'rna_NCAM1', 'rna_TYROBP', 'rna_FGFBP2', 'rna_KLRD1', 'rna_KLRF1', 'rna_KLRB1', 'rna_CX3CR1',
                         'rna_FCGR3A', 'rna_XCL1', 'rna_XCL2', 'rna_GNLY', 'rna_NKG7', 'rna_PRF1'),
             ncol=4, with_dimplot=TRUE, group.by='integrated_snn_res.0.7') # NK like
violin_plots(Tcells_3, c('rna_HAVCR2', 'rna_LAG3', 'rna_PDCD1', 'rna_CTLA4', 'rna_TIGIT', 'rna_BTLA', 'rna_KLRC1'),
             ncol=3, with_dimplot=TRUE, group.by='integrated_snn_res.0.7') # Inhibitory markers
violin_plots(Tcells_3, c('rna_HSPA6', 'rna_HSPA1B', 'rna_HSPA1A', 'rna_HSPB1', 'rna_HSP90AA1',
                          'rna_HSPD1', 'rna_HSPH1'),
              ncol=3, with_dimplot=TRUE, group.by='integrated_snn_res.0.7') # Heat-shock proteins
# Resolution 0.7 will be used.


# 3. Find markers
# 3.1. Calculate all markers
Seurat::Idents(Tcells_3) = 'integrated_snn_res.0.7'
Tcells_3_markers = Seurat::FindAllMarkers(Tcells_3, assay='RNA', slot='data', logfc.threshold=.5, only.pos=TRUE)
View(Tcells_3_markers)
write.csv(Tcells_3_markers, paste(project_dir, '2_annotation/results_Tcells/markers/markers_C3_res07.csv', sep='/'))
# 3.2. Get top 20 markers for each sub-cluster
Tcells_3_markers_top20 = c()
for(clust in as.character(unique(Tcells_3_markers$cluster))){
  clust_markers = Tcells_3_markers[Tcells_3_markers$cluster==clust,]
  n_markers = dim(clust_markers)[1]
  order_by_avg_log2FC = order(clust_markers$avg_log2FC, decreasing=T)
  if(n_markers > 20) top20_clust_markers = clust_markers[order_by_avg_log2FC[1:20],]
  else top20_clust_markers = clust_markers[order_by_avg_log2FC,]
  Tcells_3_markers_top20 = rbind(Tcells_3_markers_top20, top20_clust_markers)
}
View(Tcells_3_markers_top20)
write.csv(Tcells_3_markers_top20, paste(project_dir, '2_annotation/results_Tcells/markers/markers_C3_res08_top20.csv', sep='/'))

# 4. Annotate sub-clusters:
# Sub-clusters 0, 1, 3                      --> NK cells
# Sub-cluster 2, 4, 5, 6, 10, 11, 12, 15    --> gd Tcells
# Sub-cluster 7                             --> ILCs
# Sub-cluster 8, 9                          --> CD8 CTL? NKT?
# Sub-cluster 13                            --> NKT
# Sub-cluster 14                            --> Cytotoxic CD8aa IELs? NKT?
# Sub-cluster 16, 17                        --> CD8? NKT?
# Sub-cluster 18                            --> Heat-shock proteins
# 4.1. Insert the annotations in the metadata:
Tcells_3[['Final_Annotation']] = rep('', length(Tcells_3$integrated_snn_res.0.7))
NK_cells = rownames(Tcells_3@meta.data)[Tcells_3$integrated_snn_res.0.7%in%c('0', '1', '3')]
Tcells_3@meta.data[NK_cells, 'Final_Annotation'] = 'NK cells'
gdT_cells = rownames(Tcells_3@meta.data)[Tcells_3$integrated_snn_res.0.7%in%c('2', '4', '5', '6', '10', '11', '12', '15')]
Tcells_3@meta.data[gdT_cells, 'Final_Annotation'] = 'gd Tcells'
ILC_cells = rownames(Tcells_3@meta.data)[Tcells_3$integrated_snn_res.0.7=='7']
Tcells_3@meta.data[ILC_cells, 'Final_Annotation'] = 'ILCs'
CTLNKT_cells = rownames(Tcells_3@meta.data)[Tcells_3$integrated_snn_res.0.7%in%c('8', '9')]
Tcells_3@meta.data[CTLNKT_cells, 'Final_Annotation'] = 'CD8 CTL? NKT?'
NKT_cells = rownames(Tcells_3@meta.data)[Tcells_3$integrated_snn_res.0.7=='13']
Tcells_3@meta.data[NKT_cells, 'Final_Annotation'] = 'NKT'
CD8aaNKT_cells = rownames(Tcells_3@meta.data)[Tcells_3$integrated_snn_res.0.7=='14']
Tcells_3@meta.data[CD8aaNKT_cells, 'Final_Annotation'] = 'Cytotoxic CD8aa IELs? NKT?'
CD8NKT_cells = rownames(Tcells_3@meta.data)[Tcells_3$integrated_snn_res.0.7%in%c('16', '17')]
Tcells_3@meta.data[CD8NKT_cells, 'Final_Annotation'] = 'CD8? NKT?'
HSP_cells = rownames(Tcells_3@meta.data)[Tcells_3$integrated_snn_res.0.7=='18']
Tcells_3@meta.data[HSP_cells, 'Final_Annotation'] = 'Heat-shock proteins overexpression'
# 4.2. Check previous SIGMA result and colour cluster 0 by new sub-clusters
SIGMA_all_res02 = readRDS(paste(project_dir, '2_annotation/results_Tcells/datasets/SIGMA_all_res02.Rdata', sep='/'))
cluster_3_cells = rownames(Tcells_3@meta.data)
subclusters_3 = Tcells_3$Final_Annotation[intersect(colnames(SIGMA_all_res02$input_parameters$expr), cluster_3_cells)]
SIGMA::plot_singular_vectors(SIGMA_all_res02, '3', colour=as.character(subclusters_3))
# 4.2. Visualize annotations:
Seurat::DimPlot(Tcells_3, group.by='Final_Annotation', pt.size=0.8, label=F) +#T, label.size=4) +
  ggplot2::ggtitle('') + ggplot2::theme_minimal()# + Seurat::NoLegend()
# 4.3. Save changes:
SeuratDisk::SaveH5Seurat(Tcells_3, paste(project_dir, '2_annotation/results_Tcells/datasets/Tcells_3.h5Seurat', sep='/'), overwrite=TRUE)


# -----
# - Sub-cluster cluster 4 from resolution 0.2 separately
# -----
#Tcells = SeuratDisk::LoadH5Seurat(paste(project_dir, '2_annotation/results_Tcells/datasets/Tcells_clustered.h5Seurat', sep='/'))


# 1. Sub-cluster:
Tcells_4 = subset(Tcells, cells=rownames(Tcells@meta.data)[Tcells$integrated_snn_res.0.2=='4'])
(Seurat::DimPlot(Tcells, group.by = 'integrated_snn_res.0.2', label=T, label.size=4) + ggplot2::theme_minimal() + Seurat::NoLegend()) |
  (Seurat::DimPlot(Tcells_4, group.by = 'integrated_snn_res.0.2', label=T, label.size=4) + ggplot2::theme_minimal() + Seurat::NoLegend())
# 1.1. PCA
Tcells_4 = Seurat::RunPCA(Tcells_4, assay='integrated')
# 1.2. Choose number of PCs to use
pct = Tcells_4[["pca"]]@stdev /sum(Tcells_4[["pca"]]@stdev) * 100
cumu = cumsum(pct)
co1 = which(cumu > 90 & pct < 5)[1]
co2 = sort(which((pct[1:length(pct) - 1] - pct[2:length(pct)]) > 0.1), decreasing = T)[1] + 1
elbow = min(co1, co2) # 8
# 1.3. Visual representation of the PCs to use
plot_df = data.frame(pct=pct, cumu=cumu, rank=1:length(pct))
ggplot2::ggplot(plot_df, ggplot2::aes(cumu, pct, label = rank, color = rank > elbow)) + 
  ggplot2::geom_text() + 
  ggplot2::geom_vline(xintercept = 90, color = "grey") + 
  ggplot2::geom_hline(yintercept = min(pct[pct > 5]), color = "grey") +
  ggplot2::theme_bw()
# 1.4. Heatmap of the first 9 PCs
Seurat::DimHeatmap(Tcells_4, dims=1:9, cells=500, balanced=TRUE)
# 1.5. Determine the K-nearest neighbor graph
Tcells_4 = Seurat::FindNeighbors(Tcells_4, dims=1:elbow)
# 1.6. Find clusters:
Tcells_4 = Seurat::FindClusters(Tcells_4, resolution=seq(0.1, 1, by=.1))
# 1.7. UMAP:
Tcells_4 = Seurat::RunUMAP(Tcells_4, dims=1:elbow, reduction='pca')
invisible(gc())
# 1.8. Save object
SeuratDisk::SaveH5Seurat(Tcells_4, paste(project_dir, '2_annotation/results_Tcells/datasets/Tcells_4.h5Seurat', sep='/'), overwrite=TRUE)


# 2. Choose best resolution using clustree:
library(ggraph)
clust_tree = clustree::clustree(Tcells_4, prefix='integrated_snn_res.')
clust_tree
# 2.1. Visualize different resolutions (first resolution was chosen to be the best for now):
clusters_01_0 = Seurat::DimPlot(Tcells_4, reduction="umap", group.by='integrated_snn_res.0.1', label=TRUE, label.size=6, pt.size=.5)
clusters_03_0 = Seurat::DimPlot(Tcells_4, reduction="umap", group.by='integrated_snn_res.0.3', label=TRUE, label.size=6, pt.size=.5)
clusters_05_0 = Seurat::DimPlot(Tcells_4, reduction="umap", group.by='integrated_snn_res.0.5', label=TRUE, label.size=6, pt.size=.5)
clusters_09_0 = Seurat::DimPlot(Tcells_4, reduction="umap", group.by='integrated_snn_res.0.9', label=TRUE, label.size=6, pt.size=.5)
((clusters_01_0 | clusters_03_0) / (clusters_05_0 | clusters_09_0)) | clust_tree
# 2.2. Feature Plots to assess known gene markers:
feature_plots(Tcells_4, c('rna_CD3D', 'rna_CD3E', 'rna_CD3G', 'rna_CD8A', 'rna_CD8B', 'rna_CD4', 'rna_TRAC', 'rna_TRBC1', 'rna_TRBC2'),
              ncol=4, with_dimplot=TRUE, dimplot_group='integrated_snn_res.0.9', point_size = .5)
feature_plots(Tcells_4, c('rna_CCR7', 'rna_SELL', 'rna_PASK', 'rna_CXCL13', 'rna_ZNF683',
                           'rna_CXCL13', 'rna_CD40LG'), point_size = 1,
              ncol=4, with_dimplot=TRUE, dimplot_group='integrated_snn_res.0.9') # Naive, CM, RM
feature_plots(Tcells_4, c('rna_KLRB1', 'rna_GNLY', 'rna_NKG7', 'rna_GZMA', 'rna_GZMB', 'rna_GZMH', 'rna_GZMK', 'rna_GZMM', 'rna_PRF1',
                          'rna_IFNG'),
              ncol=3, with_dimplot=TRUE, dimplot_group='integrated_snn_res.0.9', point_size = .5) # Citotoxicity
# 2.3. Violin Plots to assess known markers:
violin_plots(Tcells_4, c('rna_CD3D', 'rna_CD3E', 'rna_CD3G', 'rna_CD8A', 'rna_CD8B', 'rna_CD4', 'rna_TRAC', 'rna_TRBC1', 'rna_TRBC2',
                         'rna_CCR7', 'rna_SELL', 'rna_PASK', 'rna_CXCL13', 'rna_ZNF683', 'rna_CXCL13', 'rna_CD40LG'),
             ncol=4, with_dimplot=TRUE, group.by='integrated_snn_res.0.9')
violin_plots(Tcells_4, c('rna_KLRB1', 'rna_GNLY', 'rna_NKG7', 'rna_GZMA', 'rna_GZMB', 'rna_GZMH', 'rna_GZMK', 'rna_GZMM', 'rna_PRF1',
                         'rna_IFNG', 'rna_NCR1', 'rna_NCAM1', 'rna_TYROBP', 'rna_FGFBP2', 'rna_KLRD1', 'rna_KLRF1', 'rna_CX3CR1',
                         'rna_FCGR3A', 'rna_XCL1', 'rna_XCL2'),
             ncol=5, with_dimplot=TRUE, group.by='integrated_snn_res.0.9') # Citotoxicity + NK like
# 2.4. Resolution 0.9 will be considered

# Find Markers
# 3.1. Calculate all markers
Seurat::Idents(Tcells_4) = 'integrated_snn_res.0.9'
Tcells_4_markers = Seurat::FindAllMarkers(Tcells_4, assay='RNA', slot='data', logfc.threshold=.5, only.pos=TRUE)
View(Tcells_4_markers)
write.csv(Tcells_4_markers, paste(project_dir, '2_annotation/results_Tcells/markers/markers_C4_res09.csv', sep='/'))
# 3.2. Get top 20 markers for each sub-cluster
Tcells_4_markers_top20 = c()
for(clust in as.character(unique(Tcells_4_markers$cluster))){
  clust_markers = Tcells_4_markers[Tcells_4_markers$cluster==clust,]
  n_markers = dim(clust_markers)[1]
  order_by_avg_log2FC = order(clust_markers$avg_log2FC, decreasing=T)
  if(n_markers > 20) top20_clust_markers = clust_markers[order_by_avg_log2FC[1:20],]
  else top20_clust_markers = clust_markers[order_by_avg_log2FC,]
  Tcells_4_markers_top20 = rbind(Tcells_4_markers_top20, top20_clust_markers)
}
View(Tcells_4_markers_top20)
write.csv(Tcells_4_markers_top20, paste(project_dir, '2_annotation/results_Tcells/markers/markers_C4_res09_top20.csv', sep='/'))

# 4. Annotate cells:
# Sub-clusters 0, 1             --> CXCL13+ CD8 Tcells?
# Sub-cluster 2, 7, 11, 13      --> Follicular CD4 Tcells
# Sub-cluster 3, 10             --> TRM CD8 Tcells
# Sub-cluster 4, 5, 8, 9, 12    --> CXCL13+ CD8 Tcells
# Sub-cluster 6                 --> Follicular CD4 Tcells?
# 4.1. Insert the annotations in the metadata:
Tcells_4[['Final_Annotation']] = rep('', length(Tcells_4$integrated_snn_res.0.9))
fh_cells = rownames(Tcells_4@meta.data)[Tcells_4$integrated_snn_res.0.9%in%c('2', '7', '11', '13')]
Tcells_4@meta.data[fh_cells, 'Final_Annotation'] = 'Follicular CD4 Tcells'
fh2_cells = rownames(Tcells_4@meta.data)[Tcells_4$integrated_snn_res.0.9=='6']
Tcells_4@meta.data[fh2_cells, 'Final_Annotation'] = 'Follicular CD4 Tcells?'
RM_cells = rownames(Tcells_4@meta.data)[Tcells_4$integrated_snn_res.0.9%in%c('3', '10')]
Tcells_4@meta.data[RM_cells, 'Final_Annotation'] = 'TRM CD8 Tcells'
CD8CXCL13_cells = rownames(Tcells_4@meta.data)[Tcells_4$integrated_snn_res.0.9%in%c('4', '5', '8', '9', '12')]
Tcells_4@meta.data[CD8CXCL13_cells, 'Final_Annotation'] = 'CXCL13+ CD8 Tcells'
CD8CXCL132_cells = rownames(Tcells_4@meta.data)[Tcells_4$integrated_snn_res.0.9%in%c('0', '1')]
Tcells_4@meta.data[CD8CXCL132_cells, 'Final_Annotation'] = 'CXCL13+ CD8 Tcells?'
# 4.2. Check previous SIGMA result and colour cluster 0 by new sub-clusters
SIGMA_all_res02 = readRDS(paste(project_dir, '2_annotation/results_Tcells/datasets/SIGMA_all_res02.Rdata', sep='/'))
cluster_4_cells = rownames(Tcells_4@meta.data)
subclusters_4 = Tcells_4$Final_Annotation[intersect(colnames(SIGMA_all_res02$input_parameters$expr), cluster_4_cells)]
SIGMA::plot_singular_vectors(SIGMA_all_res02, '4', colour=as.character(subclusters_4))
# 4.3. Visualize annotations:
Seurat::DimPlot(Tcells_4, group.by='Final_Annotation', pt.size=0.5, label=T, label.size=4) +
  ggplot2::ggtitle('') + ggplot2::theme_minimal() + Seurat::NoLegend()
# 4.4. Save changes:
SeuratDisk::SaveH5Seurat(Tcells_4, paste(project_dir, '2_annotation/results_Tcells/datasets/Tcells_4.h5Seurat', sep='/'), overwrite=TRUE)


# -----
# - Sub-cluster cluster 5 from resolution 0.2 separately
# -----
#Tcells = SeuratDisk::LoadH5Seurat(paste(project_dir, '2_annotation/results_Tcells/datasets/Tcells_clustered.h5Seurat', sep='/'))


# 1. Sub-cluster:
Tcells_5 = subset(Tcells, cells=rownames(Tcells@meta.data)[Tcells$integrated_snn_res.0.2=='5'])
(Seurat::DimPlot(Tcells, group.by = 'integrated_snn_res.0.2', label=T, label.size=4) + ggplot2::theme_minimal() + Seurat::NoLegend()) |
  (Seurat::DimPlot(Tcells_5, group.by = 'integrated_snn_res.0.2', label=T, label.size=4) + ggplot2::theme_minimal() + Seurat::NoLegend())
# 1.1. PCA
Tcells_5 = Seurat::RunPCA(Tcells_5, assay='integrated')
# 1.2. Choose number of PCs to use
pct = Tcells_5[["pca"]]@stdev /sum(Tcells_5[["pca"]]@stdev) * 100
cumu = cumsum(pct)
co1 = which(cumu > 90 & pct < 5)[1]
co2 = sort(which((pct[1:length(pct) - 1] - pct[2:length(pct)]) > 0.1), decreasing = T)[1] + 1
elbow = min(co1, co2) # 8
# 1.3. Visual representation of the PCs to use
plot_df = data.frame(pct=pct, cumu=cumu, rank=1:length(pct))
ggplot2::ggplot(plot_df, ggplot2::aes(cumu, pct, label = rank, color = rank > elbow)) + 
  ggplot2::geom_text() + 
  ggplot2::geom_vline(xintercept = 90, color = "grey") + 
  ggplot2::geom_hline(yintercept = min(pct[pct > 5]), color = "grey") +
  ggplot2::theme_bw()
# 1.4. Heatmap of the first 8 PCs
Seurat::DimHeatmap(Tcells_5, dims=1:9, cells=500, balanced=TRUE)
# 1.5. Determine the K-nearest neighbor graph
Tcells_5 = Seurat::FindNeighbors(Tcells_5, dims=1:elbow)
# 1.6. Find clusters:
Tcells_5 = Seurat::FindClusters(Tcells_5, resolution=seq(0.1, 1, by=.1))
# 1.7. UMAP:
Tcells_5 = Seurat::RunUMAP(Tcells_5, dims=1:elbow, reduction='pca')
invisible(gc())
# 1.8. Save object
SeuratDisk::SaveH5Seurat(Tcells_5, paste(project_dir, '2_annotation/results_Tcells/datasets/Tcells_5.h5Seurat', sep='/'), overwrite=TRUE)


# 2. Choose best resolution using clustree:
library(ggraph)
clust_tree = clustree::clustree(Tcells_5, prefix='integrated_snn_res.')
clust_tree
# 2.1. Visualize different resolutions (first resolution was chosen to be the best for now):
clusters_01_0 = Seurat::DimPlot(Tcells_5, reduction="umap", group.by='integrated_snn_res.0.1', label=TRUE, label.size=6, pt.size=.5)
clusters_03_0 = Seurat::DimPlot(Tcells_5, reduction="umap", group.by='integrated_snn_res.0.3', label=TRUE, label.size=6, pt.size=.5)
clusters_07_0 = Seurat::DimPlot(Tcells_5, reduction="umap", group.by='integrated_snn_res.0.7', label=TRUE, label.size=6, pt.size=.5)
clusters_08_0 = Seurat::DimPlot(Tcells_5, reduction="umap", group.by='integrated_snn_res.0.8', label=TRUE, label.size=6, pt.size=.5)
((clusters_01_0 | clusters_03_0) / (clusters_07_0 | clusters_08_0)) | clust_tree
# 2.2. Feature Plots to assess known gene markers:
feature_plots(Tcells_5, c('rna_CD3D', 'rna_CD3E', 'rna_CD3G', 'rna_CD8A', 'rna_CD8B', 'rna_CD4', 'rna_TRAC', 'rna_TRBC1', 'rna_TRBC2'),
              ncol=4, with_dimplot=TRUE, dimplot_group='integrated_snn_res.0.7', point_size = .5)
feature_plots(Tcells_5, c('rna_CCR7', 'rna_SELL', 'rna_PASK', 'rna_ZNF683'), point_size = 1,
              ncol=3, with_dimplot=TRUE, dimplot_group='integrated_snn_res.0.7') # Naive, CM, RM
feature_plots(Tcells_5, c('rna_KLRB1', 'rna_S100A4', 'rna_PTPRCAP', 'rna_PFN1', 'rna_GZMA'),
              ncol=2, with_dimplot=TRUE, dimplot_group='integrated_snn_res.0.7') # EM
feature_plots(Tcells_5, c('rna_RORC', 'rna_IL17A', 'rna_IL17F',
                          'rna_IL22', 'rna_CCR6', 'rna_KLRB1', 'rna_RORA'),
              ncol=3, with_dimplot=TRUE, dimplot_group='integrated_snn_res.0.7', point_size = .5) # Th17 like
feature_plots(Tcells_5, c('rna_TBX21', 'rna_IFNG', 'rna_STAT4', 'rna_IL12RB2'),
              ncol=3, with_dimplot=TRUE, dimplot_group='integrated_snn_res.0.7', point_size = .5) # Th1 like
feature_plots(Tcells_5, c('rna_GATA3', 'rna_IL5', 'rna_IL13', 'rna_STAT6', 'rna_IL4', 'rna_CSF2'),
              ncol=3, with_dimplot=TRUE, dimplot_group='integrated_snn_res.0.7', point_size = .5) # Th2 like
feature_plots(Tcells_5, c('rna_KLRB1', 'rna_GNLY', 'rna_NKG7', 'rna_GZMA', 'rna_GZMB', 'rna_GZMH', 'rna_GZMK', 'rna_GZMM', 'rna_PRF1',
                          'rna_IFNG'),
              ncol=3, with_dimplot=TRUE, dimplot_group='integrated_snn_res.0.7', point_size = .5) # Citotoxicity
feature_plots(Tcells_5, c('rna_TMIGD2'),
              ncol=2, with_dimplot=TRUE, dimplot_group='integrated_snn_res.0.7', point_size = .5) # ~Proliferation
# 2.3. Violin Plots to assess known markers:
violin_plots(Tcells_5, c('rna_CD3D', 'rna_CD3E', 'rna_CD3G', 'rna_CD8A', 'rna_CD8B', 'rna_CD4', 'rna_TRAC', 'rna_TRBC1', 'rna_TRBC2',
                         'rna_CCR7', 'rna_SELL', 'rna_PASK', 'rna_RORC', 'rna_IL17A', 'rna_IL17F', 'rna_IL22', 'rna_CCR6', 'rna_RORA',
                         'rna_KLRB1', 'rna_GNLY', 'rna_NKG7', 'rna_GZMA', 'rna_GZMB', 'rna_GZMH', 'rna_GZMK', 'rna_GZMM', 'rna_PRF1',
                         'rna_IFNG'),
             ncol=5, with_dimplot=TRUE, group.by='integrated_snn_res.0.7')
# 2.4. Resolution 0.7 will be considered

# Find Markers
# 3.1. Calculate all markers
Seurat::Idents(Tcells_5) = 'integrated_snn_res.0.7'
Tcells_5_markers = Seurat::FindAllMarkers(Tcells_5, assay='RNA', slot='data', logfc.threshold=.5, only.pos=TRUE)
View(Tcells_5_markers)
write.csv(Tcells_5_markers, paste(project_dir, '2_annotation/results_Tcells/markers/markers_C5_res07.csv', sep='/'))
# 3.2. Get top 20 markers for each sub-cluster
Tcells_5_markers_top20 = c()
for(clust in as.character(unique(Tcells_5_markers$cluster))){
  clust_markers = Tcells_5_markers[Tcells_5_markers$cluster==clust,]
  n_markers = dim(clust_markers)[1]
  order_by_avg_log2FC = order(clust_markers$avg_log2FC, decreasing=T)
  if(n_markers > 20) top20_clust_markers = clust_markers[order_by_avg_log2FC[1:20],]
  else top20_clust_markers = clust_markers[order_by_avg_log2FC,]
  Tcells_5_markers_top20 = rbind(Tcells_5_markers_top20, top20_clust_markers)
}
View(Tcells_5_markers_top20)
write.csv(Tcells_5_markers_top20, paste(project_dir, '2_annotation/results_Tcells/markers/markers_C5_res07_top20.csv', sep='/'))

# 4. Annotate cells:
# Sub-clusters 0, 1, 2, 3, 5, 8, 9    --> CD4 ??
# Sub-cluster 4, 7                    --> IL17A+ CD4
# Sub-cluster 6                       --> IL22+ CD4
# Sub-cluster 10                      --> IL17A+ IL17F+ CD4
# Sub-cluster 11                      --> IL17A+ IL22+ CD4
# Sub-cluster 12                      --> IL17F+ CD4
# 4.1. Insert the annotations in the metadata:
Tcells_5[['Final_Annotation']] = rep('', length(Tcells_5$integrated_snn_res.0.7))
CD4_cells = rownames(Tcells_5@meta.data)[Tcells_5$integrated_snn_res.0.7%in%c('0', '1', '2', '3', '5', '8', '9')]
Tcells_5@meta.data[CD4_cells, 'Final_Annotation'] = 'CD4 ??'
IL17A_cells = rownames(Tcells_5@meta.data)[Tcells_5$integrated_snn_res.0.7%in%c('4', '7')]
Tcells_5@meta.data[IL17A_cells, 'Final_Annotation'] = 'IL17A+ CD4'
IL22_cells = rownames(Tcells_5@meta.data)[Tcells_5$integrated_snn_res.0.7=='6']
Tcells_5@meta.data[IL22_cells, 'Final_Annotation'] = 'IL22+ CD4'
IL17AF_cells = rownames(Tcells_5@meta.data)[Tcells_5$integrated_snn_res.0.7=='10']
Tcells_5@meta.data[IL17AF_cells, 'Final_Annotation'] = 'IL17A+ IL17F+ CD4'
IL17A22_cells = rownames(Tcells_5@meta.data)[Tcells_5$integrated_snn_res.0.7=='11']
Tcells_5@meta.data[IL17A22_cells, 'Final_Annotation'] = 'IL17A+ IL22+ CD4'
IL17F_cells = rownames(Tcells_5@meta.data)[Tcells_5$integrated_snn_res.0.7=='12']
Tcells_5@meta.data[IL17F_cells, 'Final_Annotation'] = 'IL17F+ CD4'
# 4.2. Check previous SIGMA result and colour cluster 0 by new sub-clusters
SIGMA_all_res02 = readRDS(paste(project_dir, '2_annotation/results_Tcells/datasets/SIGMA_all_res02.Rdata', sep='/'))
cluster_5_cells = rownames(Tcells_5@meta.data)
subclusters_5 = Tcells_5$Final_Annotation[intersect(colnames(SIGMA_all_res02$input_parameters$expr), cluster_5_cells)]
SIGMA::plot_singular_vectors(SIGMA_all_res02, '5', colour=as.character(subclusters_5))
# 4.3. Visualize annotations:
Seurat::DimPlot(Tcells_5, group.by='Final_Annotation', pt.size=0.5, label=T, label.size=4) +
  ggplot2::ggtitle('') + ggplot2::theme_minimal() + Seurat::NoLegend()
# 4.4. Save changes:
SeuratDisk::SaveH5Seurat(Tcells_5, paste(project_dir, '2_annotation/results_Tcells/datasets/Tcells_5.h5Seurat', sep='/'), overwrite=TRUE)





# -----
# - Sub-cluster cluster 6 from resolution 0.1 separately
# -----
#Tcells = SeuratDisk::LoadH5Seurat(paste(project_dir, '2_annotation/results_Tcells/datasets/Tcells_clustered.h5Seurat', sep='/'))


# 1. Sub-cluster:
Tcells_6 = subset(Tcells, cells=rownames(Tcells@meta.data)[Tcells$integrated_snn_res.0.2=='6'])
(Seurat::DimPlot(Tcells, group.by = 'integrated_snn_res.0.2', label=T, label.size=4) + ggplot2::theme_minimal() + Seurat::NoLegend()) |
  (Seurat::DimPlot(Tcells_6, group.by = 'integrated_snn_res.0.2', label=T, label.size=4) + ggplot2::theme_minimal() + Seurat::NoLegend())
# 1.1. PCA
Tcells_6 = Seurat::RunPCA(Tcells_6, assay='integrated')
# 1.2. Choose number of PCs to use
pct = Tcells_6[["pca"]]@stdev /sum(Tcells_6[["pca"]]@stdev) * 100
cumu = cumsum(pct)
co1 = which(cumu > 90 & pct < 5)[1]
co2 = sort(which((pct[1:length(pct) - 1] - pct[2:length(pct)]) > 0.1), decreasing = T)[1] + 1
elbow = min(co1, co2) # 8
# 1.3. Visual representation of the PCs to use
plot_df = data.frame(pct=pct, cumu=cumu, rank=1:length(pct))
ggplot2::ggplot(plot_df, ggplot2::aes(cumu, pct, label = rank, color = rank > elbow)) + 
  ggplot2::geom_text() + 
  ggplot2::geom_vline(xintercept = 90, color = "grey") + 
  ggplot2::geom_hline(yintercept = min(pct[pct > 5]), color = "grey") +
  ggplot2::theme_bw()
# 1.4. Heatmap of the first 7 PCs
Seurat::DimHeatmap(Tcells_6, dims=1:9, cells=500, balanced=TRUE)
# 1.5. Determine the K-nearest neighbor graph
Tcells_6 = Seurat::FindNeighbors(Tcells_6, dims=1:elbow)
# 1.6. Find clusters:
Tcells_6 = Seurat::FindClusters(Tcells_6, resolution=seq(0.1, 1, by=.1))
# 1.7. UMAP:
Tcells_6 = Seurat::RunUMAP(Tcells_6, dims=1:elbow, reduction='pca')
invisible(gc())
# 1.8. Save object
SeuratDisk::SaveH5Seurat(Tcells_6, paste(project_dir, '2_annotation/results_Tcells/datasets/Tcells_6.h5Seurat', sep='/'), overwrite=TRUE)


# 2. Choose best resolution using clustree:
library(ggraph)
clust_tree = clustree::clustree(Tcells_6, prefix='integrated_snn_res.')
clust_tree
# 2.1. Visualize different resolutions (first resolution was chosen to be the best for now):
clusters_01_0 = Seurat::DimPlot(Tcells_6, reduction="umap", group.by='integrated_snn_res.0.1', label=TRUE, label.size=6, pt.size=.5)
clusters_03_0 = Seurat::DimPlot(Tcells_6, reduction="umap", group.by='integrated_snn_res.0.3', label=TRUE, label.size=6, pt.size=.5)
clusters_06_0 = Seurat::DimPlot(Tcells_6, reduction="umap", group.by='integrated_snn_res.0.6', label=TRUE, label.size=6, pt.size=.5)
clusters_07_0 = Seurat::DimPlot(Tcells_6, reduction="umap", group.by='integrated_snn_res.0.7', label=TRUE, label.size=6, pt.size=.5)
((clusters_01_0 | clusters_03_0) / (clusters_06_0 | clusters_07_0)) | clust_tree
# 2.2. Feature Plots to assess known gene markers:
feature_plots(Tcells_6, c('rna_CD3D', 'rna_CD3E', 'rna_CD3G', 'rna_CD8A', 'rna_CD8B', 'rna_CD4', 'rna_TRAC', 'rna_TRBC1', 'rna_TRBC2'),
              ncol=4, with_dimplot=TRUE, dimplot_group='integrated_snn_res.0.3', point_size = .5)
feature_plots(Tcells_6, c('rna_MKI67', 'rna_PCNA', 'S.Score', 'G2M.Score'),
              ncol=3, with_dimplot=TRUE, dimplot_group='integrated_snn_res.0.3', point_size = .5) # ~Proliferation
# 2.3. Violin Plots to assess known markers:
violin_plots(Tcells_6, c('rna_CD3D', 'rna_CD3E', 'rna_CD3G', 'rna_CD8A', 'rna_CD8B', 'rna_CD4', 'rna_TRAC', 'rna_TRBC1', 'rna_TRBC2',
                         'rna_MKI67', 'rna_PCNA', 'S.Score', 'G2M.Score'),
             ncol=5, with_dimplot=TRUE, group.by='integrated_snn_res.0.3')
violin_plots(Tcells, c('rna_MKI67', 'rna_PCNA', 'S.Score', 'G2M.Score'),
             ncol=3, with_dimplot=TRUE, group.by='integrated_snn_res.0.2')
# 2.4. Resolution 0.3 will be considered

# Find Markers
# 3.1. Calculate all markers
Seurat::Idents(Tcells_6) = 'integrated_snn_res.0.3'
Tcells_6_markers = Seurat::FindAllMarkers(Tcells_6, assay='RNA', slot='data', logfc.threshold=.5, only.pos=TRUE)
View(Tcells_6_markers)
write.csv(Tcells_6_markers, paste(project_dir, '2_annotation/results_Tcells/markers/markers_C6_res03.csv', sep='/'))
# 3.2. Get top 20 markers for each sub-cluster
Tcells_6_markers_top20 = c()
for(clust in as.character(unique(Tcells_6_markers$cluster))){
  clust_markers = Tcells_6_markers[Tcells_6_markers$cluster==clust,]
  n_markers = dim(clust_markers)[1]
  order_by_avg_log2FC = order(clust_markers$avg_log2FC, decreasing=T)
  if(n_markers > 20) top20_clust_markers = clust_markers[order_by_avg_log2FC[1:20],]
  else top20_clust_markers = clust_markers[order_by_avg_log2FC,]
  Tcells_6_markers_top20 = rbind(Tcells_6_markers_top20, top20_clust_markers)
}
View(Tcells_6_markers_top20)
write.csv(Tcells_6_markers_top20, paste(project_dir, '2_annotation/results_Tcells/markers/markers_C6_res03_top20.csv', sep='/'))

# 4. Annotate cells:
# Sub-cluster 0    --> proliferative CD4
# Sub-cluster 1    --> proliferative CD8
# Sub-cluster 2    --> proliferative CD8
# Sub-cluster 3    --> Mast cell Markers
# 4.1. Insert the annotations in the metadata:
Tcells_6[['Final_Annotation']] = rep('', length(Tcells_6$integrated_snn_res.0.3))
prolCD4_cells = rownames(Tcells_6@meta.data)[Tcells_6$integrated_snn_res.0.3=='0']
Tcells_6@meta.data[prolCD4_cells, 'Final_Annotation'] = 'Proliferative CD4'
prolCD8_cells = rownames(Tcells_6@meta.data)[Tcells_6$integrated_snn_res.0.3%in%c('1', '2')]
Tcells_6@meta.data[prolCD8_cells, 'Final_Annotation'] = 'Proliferative CD8'
mast_cells = rownames(Tcells_6@meta.data)[Tcells_6$integrated_snn_res.0.3=='3']
Tcells_6@meta.data[mast_cells, 'Final_Annotation'] = 'Mast cell Markers'
# 4.2. Check previous SIGMA result and colour cluster 0 by new sub-clusters
SIGMA_all_res02 = readRDS(paste(project_dir, '2_annotation/results_Tcells/datasets/SIGMA_all_res02.Rdata', sep='/'))
cluster_6_cells = rownames(Tcells_6@meta.data)
subclusters_6 = Tcells_6$Final_Annotation[intersect(colnames(SIGMA_all_res02$input_parameters$expr), cluster_6_cells)]
SIGMA::plot_singular_vectors(SIGMA_all_res02, '6', colour=as.character(subclusters_6))
# 4.3. Visualize annotations:
Seurat::DimPlot(Tcells_6, group.by='Final_Annotation', pt.size=0.5, label=T, label.size=3) +
  ggplot2::ggtitle('') + ggplot2::theme_minimal() + Seurat::NoLegend()
# 4.4. Save changes:
SeuratDisk::SaveH5Seurat(Tcells_6, paste(project_dir, '2_annotation/results_Tcells/datasets/Tcells_6.h5Seurat', sep='/'), overwrite=TRUE)





# -----
# - Annotate cluster 7 from resolution 0.2
# -----
#Tcells = SeuratDisk::LoadH5Seurat(paste(project_dir, '2_annotation/results_Tcells/datasets/Tcells_clustered.h5Seurat', sep='/'))
# This cluster will not be further sub-clustered, as clusterability was 0.


# 1. Sub-cluster:
Tcells_7 = subset(Tcells, cells=rownames(Tcells@meta.data)[Tcells$integrated_snn_res.0.2=='7'])
(Seurat::DimPlot(Tcells, group.by = 'integrated_snn_res.0.2', label=T, label.size=4) + ggplot2::theme_minimal() + Seurat::NoLegend()) |
  (Seurat::DimPlot(Tcells_7, group.by = 'integrated_snn_res.0.2', label=T, label.size=4) + ggplot2::theme_minimal() + Seurat::NoLegend())

# 2. Visualize Cluster 7:
# 2.1. Feature Plots to assess known gene markers:
feature_plots(Tcells_7, c('rna_CD3D', 'rna_CD3E', 'rna_CD3G', 'rna_CD8A', 'rna_CD8B', 'rna_CD4', 'rna_TRAC', 'rna_TRBC1', 'rna_TRBC2'),
              ncol=4, with_dimplot=TRUE, dimplot_group='integrated_snn_res.0.2', point_size = .5)
feature_plots(Tcells_7, c('rna_CCR7', 'rna_SELL', 'rna_PASK', 'rna_ZNF683'),
              ncol=3, with_dimplot=TRUE, dimplot_group='integrated_snn_res.0.2', point_size = .5) # Naive, CM, RM
# 2.3. Resolution 0.3 will be considered

# 3. Try merge cluster 7 and 0 and see what happens when clustering:
Tcells_07 = subset(Tcells, cells=rownames(Tcells@meta.data)[Tcells$integrated_snn_res.0.2%in%c('0', '7')])
Tcells_07[['previous_0_7']] = Tcells_07$integrated_snn_res.0.2
(Seurat::DimPlot(Tcells, group.by = 'integrated_snn_res.0.2', label=T, label.size=4) + ggplot2::theme_minimal() + Seurat::NoLegend()) |
  (Seurat::DimPlot(Tcells_07, group.by = 'previous_0_7', label=T, label.size=4) + ggplot2::theme_minimal() + Seurat::NoLegend())
# 3.1. PCA
Tcells_07 = Seurat::RunPCA(Tcells_07, assay='integrated')
# 3.2. Choose number of PCs to use
pct = Tcells_07[["pca"]]@stdev /sum(Tcells_07[["pca"]]@stdev) * 100
cumu = cumsum(pct)
co1 = which(cumu > 90 & pct < 5)[1]
co2 = sort(which((pct[1:length(pct) - 1] - pct[2:length(pct)]) > 0.1), decreasing = T)[1] + 1
elbow = min(co1, co2) # 5
# 3.3. Visual representation of the PCs to use
plot_df = data.frame(pct=pct, cumu=cumu, rank=1:length(pct))
ggplot2::ggplot(plot_df, ggplot2::aes(cumu, pct, label = rank, color = rank > elbow)) + 
  ggplot2::geom_text() + 
  ggplot2::geom_vline(xintercept = 90, color = "grey") + 
  ggplot2::geom_hline(yintercept = min(pct[pct > 5]), color = "grey") +
  ggplot2::theme_bw()
# 3.4. Heatmap of the first 6 PCs
Seurat::DimHeatmap(Tcells_07, dims=1:6, cells=500, balanced=TRUE)
# 3.5. Determine the K-nearest neighbor graph
Tcells_07 = Seurat::FindNeighbors(Tcells_07, dims=1:elbow)
# 3.6. Find clusters:
Tcells_07 = Seurat::FindClusters(Tcells_07, resolution=seq(0.1, 1, by=.1))
# 3.7. UMAP:
Tcells_07 = Seurat::RunUMAP(Tcells_07, dims=1:elbow, reduction='pca')
invisible(gc())
# 3.8. Save object
SeuratDisk::SaveH5Seurat(Tcells_07, paste(project_dir, '2_annotation/results_Tcells/datasets/Tcells_07.h5Seurat', sep='/'), overwrite=TRUE)


# 4. Choose best resolution using clustree:
library(ggraph)
clust_tree = clustree::clustree(Tcells_07, prefix='integrated_snn_res.')
clust_tree
# 4.1. Visualize different resolutions (first resolution was chosen to be the best for now):
clusters_01_0 = Seurat::DimPlot(Tcells_07, reduction="umap", group.by='integrated_snn_res.0.1', label=TRUE, label.size=6, pt.size=.5)
clusters_02_0 = Seurat::DimPlot(Tcells_07, reduction="umap", group.by='integrated_snn_res.0.2', label=TRUE, label.size=6, pt.size=.5)
clusters_09_0 = Seurat::DimPlot(Tcells_07, reduction="umap", group.by='integrated_snn_res.0.9', label=TRUE, label.size=6, pt.size=.5)
clusters_07_orig = Seurat::DimPlot(Tcells_07, reduction="umap", group.by='previous_0_7', label=TRUE, label.size=6, pt.size=.5)
((clusters_01_0 | clusters_02_0) / (clusters_09_0 | clusters_07_orig)) | clust_tree
# 4.2. Feature Plots to assess known gene markers:
feature_plots(Tcells_07, c('rna_CD3D', 'rna_CD3E', 'rna_CD3G', 'rna_CD8A', 'rna_CD8B', 'rna_CD4', 'rna_TRAC', 'rna_TRBC1', 'rna_TRBC2'),
              ncol=4, with_dimplot=TRUE, dimplot_group='integrated_snn_res.0.1', point_size = .5)
feature_plots(Tcells_07, c('rna_CCR7', 'rna_SELL', 'rna_PASK', 'rna_ZNF683'),
              ncol=3, with_dimplot=TRUE, dimplot_group='integrated_snn_res.0.1', point_size = .5) # Naive, CM, RM
feature_plots(Tcells_07, c('rna_KLRB1', 'rna_S100A4', 'rna_PTPRCAP', 'rna_PFN1', 'rna_GZMA'),
              ncol=2, with_dimplot=TRUE, dimplot_group='integrated_snn_res.0.1') # EM
feature_plots(Tcells_07, c('rna_STAT4', 'rna_IL12RB2', 'rna_IFNG', 'rna_GATA3', 'rna_STAT6', 'rna_IL4', 'rna_CXCL13', 'rna_IL17A'),
              ncol=3, with_dimplot=TRUE, dimplot_group='integrated_snn_res.0.1')
# 2.3. Violin Plots to assess known markers:
violin_plots(Tcells_07, c('rna_CCR7', 'rna_SELL', 'rna_PASK', 'rna_KLRB1', 'rna_S100A4', 'rna_PTPRCAP', 'rna_PFN1', 'rna_GZMA'),
             ncol=3, with_dimplot=TRUE, group.by='integrated_snn_res.0.1')

# Find Markers
# 3.1. Calculate all markers
Seurat::Idents(Tcells_07) = 'integrated_snn_res.0.1'
Tcells_7_markers = Seurat::FindAllMarkers(Tcells_07, assay='RNA', slot='data', logfc.threshold=.5, only.pos=TRUE)
View(Tcells_7_markers)
write.csv(Tcells_7_markers, paste(project_dir, '2_annotation/results_Tcells/markers/markers_C07_res01.csv', sep='/'))
# 3.2. Get top 20 markers for each sub-cluster
Tcells_7_markers_top20 = c()
for(clust in as.character(unique(Tcells_7_markers$cluster))){
  clust_markers = Tcells_7_markers[Tcells_7_markers$cluster==clust,]
  n_markers = dim(clust_markers)[1]
  order_by_avg_log2FC = order(clust_markers$avg_log2FC, decreasing=T)
  if(n_markers > 20) top20_clust_markers = clust_markers[order_by_avg_log2FC[1:20],]
  else top20_clust_markers = clust_markers[order_by_avg_log2FC,]
  Tcells_7_markers_top20 = rbind(Tcells_7_markers_top20, top20_clust_markers)
}
View(Tcells_7_markers_top20)
write.csv(Tcells_7_markers_top20, paste(project_dir, '2_annotation/results_Tcells/markers/markers_C07_res01_top20.csv', sep='/'))

# 4. Annotate cells:
# Sub-cluster 0    --> Effector Memory (EM) CD4
# Sub-cluster 1    --> Naive CD4
# Sub-cluster 2    --> Central Memory (CM) CD4
# 4.1. Insert the annotations in the metadata:
Tcells_07[['Final_Annotation']] = rep('', length(Tcells_07$integrated_snn_res.0.1))
EM_cells = rownames(Tcells_07@meta.data)[Tcells_07$integrated_snn_res.0.1=='0']
Tcells_07@meta.data[EM_cells, 'Final_Annotation'] = 'EM CD4'
naive_cells = rownames(Tcells_07@meta.data)[Tcells_07$integrated_snn_res.0.1=='1']
Tcells_07@meta.data[naive_cells, 'Final_Annotation'] = 'Naive CD4'
CM_cells = rownames(Tcells_07@meta.data)[Tcells_07$integrated_snn_res.0.1=='2']
Tcells_07@meta.data[CM_cells, 'Final_Annotation'] = 'CM CD4'
# 4.2. Check previous SIGMA result and colour cluster 0 by new sub-clusters
SIGMA_all_res02 = readRDS(paste(project_dir, '2_annotation/results_Tcells/datasets/SIGMA_all_res02.Rdata', sep='/'))
cluster_0_cells = rownames(Tcells_07@meta.data)[Tcells_07$previous_0_7=='0']
cluster_7_cells = rownames(Tcells_07@meta.data)[Tcells_07$previous_0_7=='7']
subclusters_0 = Tcells_07$Final_Annotation[intersect(colnames(SIGMA_all_res02$input_parameters$expr), cluster_0_cells)]
subclusters_7 = Tcells_07$Final_Annotation[intersect(colnames(SIGMA_all_res02$input_parameters$expr), cluster_7_cells)]
SIGMA::plot_singular_vectors(SIGMA_all_res02, '0', colour=as.character(subclusters_0))
SIGMA::plot_singular_vectors(SIGMA_all_res02, '7', colour=as.character(subclusters_7))
# 4.3. Visualize annotations:
Seurat::DimPlot(Tcells_07, group.by='Final_Annotation', pt.size=0.5, label=T, label.size=3) +
  ggplot2::ggtitle('') + ggplot2::theme_minimal() + Seurat::NoLegend()
# 4.4. Compare these annotations with previously done with only 0:
Tcells_0 = SeuratDisk::LoadH5Seurat(paste(project_dir, '2_annotation/results_Tcells/datasets/Tcells_0.h5Seurat', sep='/'), overwrite=TRUE)
Tcells_0[['Final_annotation2']] = Tcells_07$Final_Annotation[Tcells_07$previous_0_7=='0']
Seurat::DimPlot(Tcells_0, group.by='Final_Annotation', pt.size=0.5, label=T, label.size=3) |
  Seurat::DimPlot(Tcells_0, group.by='Final_annotation2', pt.size=0.5, label=T, label.size=3)
feature_plots(Tcells_0, c('rna_CCR7', 'rna_SELL', 'rna_PASK', 'rna_ZNF683'),
              ncol=3, with_dimplot=TRUE, dimplot_group='Final_annotation2', point_size = .5) # Naive, CM, RM
# 4.5. Save changes:
SeuratDisk::SaveH5Seurat(Tcells_07, paste(project_dir, '2_annotation/results_Tcells/datasets/Tcells_07.h5Seurat', sep='/'), overwrite=TRUE)





# -----
# - Analyse all clusters together
# -----
#Tcells = SeuratDisk::LoadH5Seurat(paste(project_dir, '2_annotation/results_Tcells/datasets/Tcells_clustered.h5Seurat', sep='/'))


# 1. Change metadata according to clusters performed individually:
Tcells[['Temporary_Annotation']] = rep('', dim(Tcells@meta.data)[1])
Tcells[['big_Tcell_clust']] = rep('', dim(Tcells@meta.data)[1])
Tcells[['orig_clusts']] = rep('', dim(Tcells@meta.data)[1])
# 2.1. C1
Tcells_1 = SeuratDisk::LoadH5Seurat(paste(project_dir, '2_annotation/results_Tcells/datasets/Tcells_1.h5Seurat', sep='/'), overwrite=TRUE)
Tcells@meta.data[rownames(Tcells_1@meta.data), 'Temporary_Annotation'] = paste('C1', Tcells_1$Final_Annotation, sep='_')
Tcells@meta.data[rownames(Tcells_1@meta.data), 'big_Tcell_clust'] = rep('1', dim(Tcells_1@meta.data)[1])
Tcells@meta.data[rownames(Tcells_1@meta.data), 'orig_clusts'] = paste('C1', Tcells_1$integrated_snn_res.0.9, sep='_')
# 2.2. C2
Tcells_2 = SeuratDisk::LoadH5Seurat(paste(project_dir, '2_annotation/results_Tcells/datasets/Tcells_2.h5Seurat', sep='/'), overwrite=TRUE)
Tcells@meta.data[rownames(Tcells_2@meta.data), 'Temporary_Annotation'] = paste('C2', Tcells_2$Final_Annotation, sep='_')
Tcells@meta.data[rownames(Tcells_2@meta.data), 'big_Tcell_clust'] = rep('2', dim(Tcells_2@meta.data)[1])
Tcells@meta.data[rownames(Tcells_2@meta.data), 'orig_clusts'] =  paste('C2', rep('6', dim(Tcells_2@meta.data)[1]), sep='_')
# 2.3. C3
Tcells_3 = SeuratDisk::LoadH5Seurat(paste(project_dir, '2_annotation/results_Tcells/datasets/Tcells_3.h5Seurat', sep='/'), overwrite=TRUE)
Tcells@meta.data[rownames(Tcells_3@meta.data), 'Temporary_Annotation'] = paste('C3', Tcells_3$Final_Annotation, sep='_')
Tcells@meta.data[rownames(Tcells_3@meta.data), 'big_Tcell_clust'] = rep('3', dim(Tcells_3@meta.data)[1])
Tcells@meta.data[rownames(Tcells_3@meta.data), 'orig_clusts'] = paste('C3', Tcells_3$integrated_snn_res.0.7, sep='_')
# 2.4. C4
Tcells_4 = SeuratDisk::LoadH5Seurat(paste(project_dir, '2_annotation/results_Tcells/datasets/Tcells_4.h5Seurat', sep='/'), overwrite=TRUE)
Tcells@meta.data[rownames(Tcells_4@meta.data), 'Temporary_Annotation'] = paste('C4', Tcells_4$Final_Annotation, sep='_')
Tcells@meta.data[rownames(Tcells_4@meta.data), 'big_Tcell_clust'] = rep('4', dim(Tcells_4@meta.data)[1])
Tcells@meta.data[rownames(Tcells_4@meta.data), 'orig_clusts'] = paste('C4', Tcells_4$integrated_snn_res.0.9, sep='_')
# 2.5. C5
Tcells_5 = SeuratDisk::LoadH5Seurat(paste(project_dir, '2_annotation/results_Tcells/datasets/Tcells_5.h5Seurat', sep='/'), overwrite=TRUE)
Tcells@meta.data[rownames(Tcells_5@meta.data), 'Temporary_Annotation'] = paste('C5', Tcells_5$Final_Annotation, sep='_')
Tcells@meta.data[rownames(Tcells_5@meta.data), 'big_Tcell_clust'] = rep('5', dim(Tcells_5@meta.data)[1])
Tcells@meta.data[rownames(Tcells_5@meta.data), 'orig_clusts'] = paste('C5', Tcells_5$integrated_snn_res.0.7, sep='_')
# 2.6. C6
Tcells_6 = SeuratDisk::LoadH5Seurat(paste(project_dir, '2_annotation/results_Tcells/datasets/Tcells_6.h5Seurat', sep='/'), overwrite=TRUE)
Tcells@meta.data[rownames(Tcells_6@meta.data), 'Temporary_Annotation'] = paste('C6', Tcells_6$Final_Annotation, sep='_')
Tcells@meta.data[rownames(Tcells_6@meta.data), 'big_Tcell_clust'] = rep('6', dim(Tcells_6@meta.data)[1])
Tcells@meta.data[rownames(Tcells_6@meta.data), 'orig_clusts'] = paste('C6', Tcells_6$integrated_snn_res.0.3, sep='_')
# 2.7. C07
Tcells_07 = SeuratDisk::LoadH5Seurat(paste(project_dir, '2_annotation/results_Tcells/datasets/Tcells_07.h5Seurat', sep='/'), overwrite=TRUE)
Tcells@meta.data[rownames(Tcells_07@meta.data), 'Temporary_Annotation'] = paste('C07', Tcells_07$Final_Annotation, sep='_')
Tcells@meta.data[rownames(Tcells_07@meta.data), 'big_Tcell_clust'] = rep('07', dim(Tcells_07@meta.data)[1])
Tcells@meta.data[rownames(Tcells_07@meta.data), 'orig_clusts'] = paste('C07', Tcells_07$integrated_snn_res.0.1, sep='_')
# 2.6. Remove [integrated...] metadata variables
Tcells@meta.data = Tcells@meta.data[, !colnames(Tcells@meta.data) %in%
                                                c(grep('^in', colnames(Tcells@meta.data), value=T), 'seurat_clusters')]
invisible(gc())
# 2.7. Remove individual datasets:
remove(Tcells_1, Tcells_2, Tcells_3, Tcells_4, Tcells_5, Tcells_6, Tcells_07)
invisible(gc())
# 2.8. Save this temporary Tcells dataset:
SeuratDisk::SaveH5Seurat(Tcells, paste(project_dir, '2_annotation/results_Tcells/datasets/Tcells_temporaryAnnots.h5Seurat', sep='/'),
                         overwrite=T)


# 3. Visualize the clusters:
orig_clusts = Seurat::DimPlot(Tcells, group.by='orig_clusts', label=T, label.size=4, pt.size=.5) +
  ggplot2::theme_minimal() + Seurat::NoLegend()
Temporary_Annotation = Seurat::DimPlot(Tcells, group.by='Temporary_Annotation', label=T, label.size=4, pt.size=.5) +
  ggplot2::theme_minimal() + Seurat::NoLegend()
orig_clusts | Temporary_Annotation

indiv_plots = list()
for(clust in unique(Tcells$big_Tcell_clust)){
  if(clust=='2') next
  indiv_plots[[clust]] = Seurat::DimPlot(subset(Tcells, cells=rownames(Tcells@meta.data)[Tcells$big_Tcell_clust==clust]),
                                         group.by='Temporary_Annotation', label=T, label.size=4, pt.size=.5) +
    ggplot2::theme_minimal() + ggplot2::ggtitle('') + Seurat::NoLegend()
}
patchwork::wrap_plots(indiv_plots, ncol=3)
indiv_plots$`1`
indiv_plots$`3`
indiv_plots$`4`

# 4. Calculate similarity between clusters
# 4.1 Calculate average expression profiles of each cluster:
n_genes = dim(Tcells@assays$RNA@data)[1]
n_clusters = length(unique(Tcells@meta.data[, 'orig_clusts']))
average_profiles = matrix(rep(0, n_genes * n_clusters), nrow = n_genes)
colnames(average_profiles) = unique(Tcells@meta.data[, 'orig_clusts'])
rownames(average_profiles) = rownames(Tcells@assays$RNA@data)
for(clust in unique(Tcells@meta.data[, 'orig_clusts'])){
  message('Cluster:', clust)
  clust_cells = rownames(Tcells@meta.data)[Tcells@meta.data[, 'orig_clusts'] == clust]
  clust_matrix = as.matrix(Tcells@assays$RNA@data[,clust_cells])
  average_profiles[,clust] = rowMeans(clust_matrix)
  invisible(gc())
}
# 4.2. Calculate distance matrices between clusters:
dist_average_clusts = dist(t(average_profiles), method='euclidean')
# 4.3. Visualize distance matrices:
dist_average_clusts_matrix = as.matrix(dist_average_clusts)
#dist_average_clusts_matrix[upper.tri(dist_average_clusts_matrix)] <- NA
pheatmap::pheatmap(dist_average_clusts_matrix, cluster_rows=F, cluster_cols=F, na_col="white", display_numbers=TRUE,
                   main='Clusters distance based on average profiles')
# 4.4. Create dendogram:
#hclust_average_clusters = hclust(dist_average_clusts, method='complete')
#plot(hclust_average_clusters, main='Clustering based on average profiles, cells named by original clusters')
op <- par(mfrow = c(2, 1))
hc.result = hclust_average_clusters
plot(hc.result, main='Clustering based on average profiles', xlab='original clusters')
samps_cols = unique(Tcells@meta.data[,c('orig_clusts', 'Temporary_Annotation')])
hc.result$labels = samps_cols$Temporary_Annotation
plot(hc.result, main='', xlab='Temporary Annotation')
par(op)


# 5. Cluster C5 and C07 cells
Tcells_C5C07 = subset(Tcells, cells=rownames(Tcells@meta.data)[Tcells$big_Tcell_clust%in%c('5', '07')])
(Seurat::DimPlot(Tcells, group.by = 'big_Tcell_clust', label=T, label.size=4) + ggplot2::theme_minimal() + Seurat::NoLegend()) |
  (Seurat::DimPlot(Tcells_C5C07, group.by = 'big_Tcell_clust', label=T, label.size=4) + ggplot2::theme_minimal() + Seurat::NoLegend())
# 5.1. PCA
Tcells_C5C07 = Seurat::RunPCA(Tcells_C5C07, assay='integrated')
# 5.2. Choose number of PCs to use
pct = Tcells_C5C07[["pca"]]@stdev /sum(Tcells_C5C07[["pca"]]@stdev) * 100
cumu = cumsum(pct)
co1 = which(cumu > 90 & pct < 5)[1]
co2 = sort(which((pct[1:length(pct) - 1] - pct[2:length(pct)]) > 0.1), decreasing = T)[1] + 1
elbow = min(co1, co2) # 8
# 5.3. UMAP:
Tcells_C5C07 = Seurat::RunUMAP(Tcells_C5C07, dims=1:elbow, reduction='pca')
invisible(gc())
# 5.4. Visualize clusters:
Seurat::DimPlot(Tcells_C5C07, reduction="umap", group.by='orig_clusts', label=TRUE, label.size=3, pt.size=.5, repel = T) |
  Seurat::DimPlot(Tcells_C5C07, reduction="umap", group.by='Temporary_Annotation', label=TRUE, label.size=3, pt.size=.5, repel = T)
# 6.5. Visualize specific genes:
feature_plots(Tcells_C5C07, c('rna_CXCL13', 'rna_ZNF683', 'rna_CCR7', 'rna_SELL', 'rna_PASK', 'rna_KLRB1', 'rna_S100A4', 'rna_PTPRCAP',
                              'rna_PFN1', 'rna_GZMA'),
              ncol=3, with_dimplot=TRUE, dimplot_group='orig_clusts')
violin_plots(Tcells_C5C07, c('rna_CCR7', 'rna_SELL', 'rna_PASK', 'rna_KLRB1', 'rna_S100A4', 'rna_PTPRCAP', 'rna_PFN1', 'rna_GZMA'),
             ncol=3, with_dimplot=TRUE, group.by='orig_clusts')
# Decision: C5_0, C5_1, C5_2, C5_3, C5_5, C5_8, C5_9 will be considered a 'second' type of EM Tcells (EM_2). C07_0 will be EM_1.  


# 6. Cluster C5, C07, 4 cells
Tcells_C5C074 = subset(Tcells, cells=rownames(Tcells@meta.data)[Tcells$big_Tcell_clust%in%c('5', '07', '4')])
(Seurat::DimPlot(Tcells, group.by = 'big_Tcell_clust', label=T, label.size=4) + ggplot2::theme_minimal() + Seurat::NoLegend()) |
  (Seurat::DimPlot(Tcells_C5C074, group.by = 'big_Tcell_clust', label=T, label.size=4) + ggplot2::theme_minimal() + Seurat::NoLegend())
# 6.1. PCA
Tcells_C5C074 = Seurat::RunPCA(Tcells_C5C074, assay='integrated')
# 6.2. Choose number of PCs to use
pct = Tcells_C5C074[["pca"]]@stdev /sum(Tcells_C5C074[["pca"]]@stdev) * 100
cumu = cumsum(pct)
co1 = which(cumu > 90 & pct < 5)[1]
co2 = sort(which((pct[1:length(pct) - 1] - pct[2:length(pct)]) > 0.1), decreasing = T)[1] + 1
elbow = min(co1, co2) # 8
# 6.3. UMAP:
Tcells_C5C074 = Seurat::RunUMAP(Tcells_C5C074, dims=1:elbow, reduction='pca')
invisible(gc())
# 6.4.Visualize clusters:
Seurat::DimPlot(Tcells_C5C074, reduction="umap", group.by='orig_clusts', label=TRUE, label.size=3, pt.size=.5, repel = T) |
  Seurat::DimPlot(Tcells_C5C074, reduction="umap", group.by='Temporary_Annotation', label=TRUE, label.size=3, pt.size=.5, repel = T)
# 6.5. Visualize specific genes:
feature_plots(Tcells_C5C074, c('rna_CXCL13', 'rna_ZNF683', 'rna_CCR7', 'rna_SELL', 'rna_PASK', 'rna_KLRB1', 'rna_S100A4', 'rna_PTPRCAP',
                              'rna_PFN1', 'rna_GZMA', 'rna_CD4', 'rna_CD8A', 'rna_CD8B'),
              ncol=3, with_dimplot=TRUE, dimplot_group='Temporary_Annotation')
feature_plots(Tcells_C5C074, c('rna_CXCL13', 'rna_ZNF683', 'rna_CD4', 'rna_CD8A', 'rna_CD8B'),
              ncol=3, with_dimplot=TRUE, dimplot_group='Temporary_Annotation')
violin_plots(Tcells_C5C074, c('rna_CXCL13', 'rna_ZNF683', 'rna_CCR7', 'rna_SELL', 'rna_PASK', 'rna_KLRB1', 'rna_S100A4', 'rna_PTPRCAP',
                              'rna_PFN1', 'rna_GZMA', 'rna_CD4', 'rna_CD8A', 'rna_CD8B'),
             ncol=4, with_dimplot=TRUE, group.by='Temporary_Annotation')
violin_plots(Tcells_C5C074, c('rna_CXCL13', 'rna_ZNF683', 'rna_CCR7', 'rna_SELL', 'rna_PASK', 'rna_KLRB1', 'rna_S100A4', 'rna_PTPRCAP',
                              'rna_PFN1', 'rna_GZMA', 'rna_CD4', 'rna_CD8A', 'rna_CD8B'),
             ncol=3, with_dimplot=TRUE, group.by='orig_clusts')
# Decision: C4_6 really is a Follicular CD4 Tcell cluster. C4_0 and C4_1 really are a CXCL13+ CD8 Tcells cluster.


# 7. Cluster C5, C07, C4_2, C4_6, C4_7, C4_11, C4_13 cells
Tcells_C5C02s4 = subset(Tcells, cells=rownames(Tcells@meta.data)[Tcells$big_Tcell_clust%in%c('5', '07') |
                                                                  Tcells$orig_clusts%in%c('C4_2', 'C4_6', 'C4_7', 'C4_11', 'C4_13')])
(Seurat::DimPlot(Tcells, group.by = 'big_Tcell_clust', label=T, label.size=4) + ggplot2::theme_minimal() + Seurat::NoLegend()) |
  (Seurat::DimPlot(Tcells_C5C02s4, group.by = 'big_Tcell_clust', label=T, label.size=4) + ggplot2::theme_minimal() + Seurat::NoLegend())
# 7.1. PCA
Tcells_C5C02s4 = Seurat::RunPCA(Tcells_C5C02s4, assay='integrated')
# 7.2. Choose number of PCs to use
pct = Tcells_C5C02s4[["pca"]]@stdev /sum(Tcells_C5C02s4[["pca"]]@stdev) * 100
cumu = cumsum(pct)
co1 = which(cumu > 90 & pct < 5)[1]
co2 = sort(which((pct[1:length(pct) - 1] - pct[2:length(pct)]) > 0.1), decreasing = T)[1] + 1
elbow = min(co1, co2) # 8
# 7.3. UMAP:
Tcells_C5C02s4 = Seurat::RunUMAP(Tcells_C5C02s4, dims=1:elbow, reduction='pca')
invisible(gc())
# 7.4. Visualize clusters:
Seurat::DimPlot(Tcells_C5C02s4, reduction="umap", group.by='orig_clusts', label=TRUE, label.size=3, pt.size=.5, repel = T) |
  Seurat::DimPlot(Tcells_C5C02s4, reduction="umap", group.by='Temporary_Annotation', label=TRUE, label.size=3, pt.size=.5, repel = T)
# 7.5. Visualize specific genes:
feature_plots(Tcells_C5C02s4, c('rna_CXCL13', 'rna_ZNF683', 'rna_CCR7', 'rna_SELL', 'rna_PASK', 'rna_KLRB1', 'rna_S100A4', 'rna_PTPRCAP',
                               'rna_PFN1', 'rna_GZMA', 'rna_CD4', 'rna_CD8A', 'rna_CD8B'),
              ncol=3, with_dimplot=TRUE, dimplot_group='Temporary_Annotation')
feature_plots(Tcells_C5C02s4, c('rna_CXCL13', 'rna_ZNF683', 'rna_CD4', 'rna_CD8A', 'rna_CD8B'),
              ncol=3, with_dimplot=TRUE, dimplot_group='Temporary_Annotation')
violin_plots(Tcells_C5C02s4, c('rna_CXCL13', 'rna_ZNF683', 'rna_CCR7', 'rna_SELL', 'rna_PASK', 'rna_KLRB1', 'rna_S100A4', 'rna_PTPRCAP',
                              'rna_PFN1', 'rna_GZMA', 'rna_CD4', 'rna_CD8A', 'rna_CD8B'),
             ncol=4, with_dimplot=TRUE, group.by='Temporary_Annotation')
violin_plots(Tcells_C5C02s4, c('rna_CXCL13', 'rna_ZNF683', 'rna_CCR7', 'rna_SELL', 'rna_PASK', 'rna_KLRB1', 'rna_S100A4', 'rna_PTPRCAP',
                              'rna_PFN1', 'rna_GZMA', 'rna_CD4', 'rna_CD8A', 'rna_CD8B'),
             ncol=3, with_dimplot=TRUE, group.by='orig_clusts')
# Final check that C4_6 really is a Follicular CD4 Tcell cluster and that C5_0, C5_1, C5_2, C5_3, C5_5, C5_8, C5_9 are EM, as they cluster
# accordingly


# 10. Cluster C1, C3, C4_0, C4_1, C4_3, C4_10 cells
Tcells_C1C3s4 = subset(Tcells, cells=rownames(Tcells@meta.data)[(Tcells$big_Tcell_clust%in%c('1', '3') |
                                                                   Tcells$orig_clusts%in%c('C4_0', 'C4_1', 'C4_3', 'C4_10'))
                                                                & (!Tcells$orig_clusts%in%c('C3_18', 'C3_2', 'C3_4', 'C3_5', 'C3_6',
                                                                                            'C3_10', 'C3_11', 'C3_12', 'C3_15'))])
(Seurat::DimPlot(Tcells, group.by = 'big_Tcell_clust', label=T, label.size=4) + ggplot2::theme_minimal() + Seurat::NoLegend()) |
  (Seurat::DimPlot(Tcells_C1C3s4, group.by = 'big_Tcell_clust', label=T, label.size=4) + ggplot2::theme_minimal() + Seurat::NoLegend())
# 10.1. PCA
Tcells_C1C3s4 = Seurat::RunPCA(Tcells_C1C3s4, assay='integrated')
# 10.2. Choose number of PCs to use
pct = Tcells_C1C3s4[["pca"]]@stdev /sum(Tcells_C1C3s4[["pca"]]@stdev) * 100
cumu = cumsum(pct)
co1 = which(cumu > 90 & pct < 5)[1]
co2 = sort(which((pct[1:length(pct) - 1] - pct[2:length(pct)]) > 0.1), decreasing = T)[1] + 1
elbow = min(co1, co2) # 8
# 10.3. UMAP:
Tcells_C1C3s4 = Seurat::RunUMAP(Tcells_C1C3s4, dims=1:elbow, reduction='pca')
invisible(gc())
# 10.4. Visualize clusters:
Seurat::DimPlot(Tcells_C1C3s4, reduction="umap", group.by='orig_clusts', label=TRUE, label.size=3, pt.size=.5, repel = T) |
  Seurat::DimPlot(Tcells_C1C3s4, reduction="umap", group.by='Temporary_Annotation', label=TRUE, label.size=3, pt.size=.5, repel = T)
# 10.5. Highlight specific groups of cells:
NK_cells = Seurat::DimPlot(Tcells_C1C3s4, reduction="umap",
                           cells.highlight = rownames(Tcells_C1C3s4@meta.data)[Tcells_C1C3s4$Temporary_Annotation=='C3_NK cells'],
                           pt.size=.5, repel=T) + ggplot2::ggtitle('NK')
NKT_cells = Seurat::DimPlot(Tcells_C1C3s4, reduction="umap",
                           cells.highlight = rownames(Tcells_C1C3s4@meta.data)[Tcells_C1C3s4$Temporary_Annotation=='C3_NKT'],
                          pt.size=.5, repel=T) + ggplot2::ggtitle('NKT')
pNKT_cells = Seurat::DimPlot(Tcells_C1C3s4, reduction="umap",
                             cells.highlight = rownames(Tcells_C1C3s4@meta.data)[Tcells_C1C3s4$Temporary_Annotation%in%
                                                                                   c('C3_CD8 CTL? NKT?', 'C3_CD8? NKT?',
                                                                                     'C3_Cytotoxic CD8aa IELs? NKT?',
                                                                                     'C1_DN Tcells? NKT cells?')],
                             pt.size=.5, repel=T) + ggplot2::ggtitle('NKT?')
cyto_cells = Seurat::DimPlot(Tcells_C1C3s4, reduction="umap",
                             cells.highlight = rownames(Tcells_C1C3s4@meta.data)[Tcells_C1C3s4$Temporary_Annotation%in%
                                                                                   c('C1_cytotoxic CD8 Tcells',
                                                                                     'C1_cytotoxic CD8aa IELs')],
                             pt.size=.5, repel=T) + ggplot2::ggtitle('Cytotoxic CD8 (aa and ab)')
EM_cells = Seurat::DimPlot(Tcells_C1C3s4, reduction="umap",
                             cells.highlight = rownames(Tcells_C1C3s4@meta.data)[Tcells_C1C3s4$Temporary_Annotation%in%
                                                                                   c('C1_EM CD8 Tcells')],
                             pt.size=.5, repel=T) + ggplot2::ggtitle('EM CD8')
temp_anots = Seurat::DimPlot(Tcells_C1C3s4, reduction="umap", group.by='Temporary_Annotation',
                             label=TRUE, label.size=3, pt.size=.5, repel=T) + Seurat::NoLegend()
patchwork::wrap_plots(list(temp_anots, NK_cells, NKT_cells, pNKT_cells, cyto_cells, EM_cells), ncol=2)
# 10.6. Visualize specific genes:
feature_plots(Tcells_C1C3s4, c('rna_CXCL13', 'rna_ZNF683', 'rna_CCR7', 'rna_SELL', 'rna_PASK', 'rna_CD8A', 'rna_CD8B'),
              ncol=3, with_dimplot=TRUE, dimplot_group='Temporary_Annotation')
violin_plots(Tcells_C1C3s4, c('rna_CD3E', 'rna_CD3D', 'rna_CD3G', 'rna_CD4', 'rna_CD8A', 'rna_CD8B', 'rna_TRAC', 'rna_TRBC1', 'rna_TRBC2'),
             ncol=4, with_dimplot=TRUE, group.by='Temporary_Annotation')
violin_plots(Tcells_C1C3s4, c('rna_CXCL13', 'rna_ZNF683', 'rna_CCR7', 'rna_SELL', 'rna_PASK', 'rna_KLRB1', 'rna_S100A4', 'rna_PTPRCAP',
                              'rna_PFN1', 'rna_GZMA'),
             ncol=4, with_dimplot=TRUE, group.by='Temporary_Annotation')
violin_plots(Tcells_C1C3s4, c('rna_KLRB1', 'rna_GNLY', 'rna_NKG7', 'rna_GZMA', 'rna_GZMB', 'rna_GZMH', 'rna_GZMK', 'rna_GZMM', 'rna_PRF1',
                              'rna_IFNG'),
             ncol=4, with_dimplot=TRUE, group.by='Temporary_Annotation')
violin_plots(Tcells_C1C3s4, c('rna_NCR1', 'rna_NCAM1', 'rna_TYROBP', 'rna_FGFBP2', 'rna_KLRD1', 'rna_KLRF1', 'rna_KLRB1', 'rna_CX3CR1',
                              'rna_FCGR3A', 'rna_XCL1', 'rna_XCL2', 'rna_GNLY', 'rna_NKG7', 'rna_PRF1'),
             ncol=4, with_dimplot=TRUE, group.by='Temporary_Annotation')
# Decision: C1_7, C4_3, C4_10 are TRM CD8 Tcells. The following are NKT cells: C3_9, C3_13, C3_14 (NKT_1), C3_8, C3_16, C3_17 (NKT_2),
#           C1_10 (NKT_3)





# -----
# - Set Final Annotations
# -----
#Tcells = SeuratDisk::LoadH5Seurat(paste(project_dir, '2_annotation/results_Tcells/datasets/Tcells_temporaryAnnots.h5Seurat', sep='/'))


# 1. Remove unwanted cells (C6_3: Cells with mast cell markers; C3_18: cells with too many heat-shock protein markers):
Tcells = subset(Tcells, cells=rownames(Tcells@meta.data)[!Tcells$orig_clusts%in%c('C6_3', 'C3_18')])


# 2. Annotate cells - Level 6
Tcells[['Annotation_Level_6']] = rep('Unkown', dim(Tcells@meta.data)[1])
# 2.1. Proliferative CD4 Tcells
Tcells@meta.data[Tcells$orig_clusts=='C6_0', 'Annotation_Level_6'] = 'Proliferative CD4 Tcells'
# 2.2. Regulatory CD4 Tcells
Tcells@meta.data[Tcells$orig_clusts=='C2_6', 'Annotation_Level_6'] = 'Regulatory CD4 Tcells'
# 2.3. Naive CD4 Tcells
Tcells@meta.data[Tcells$orig_clusts=='C07_1', 'Annotation_Level_6'] = 'Naive CD4 Tcells'
# 2.4. CM CD4 Tcells
Tcells@meta.data[Tcells$orig_clusts=='C07_2', 'Annotation_Level_6'] = 'CM CD4 Tcells'
# 2.5. EM_1 CD4 Tcells
Tcells@meta.data[Tcells$orig_clusts=='C07_0', 'Annotation_Level_6'] = 'EM_1 CD4 Tcells'
# 2.6. EM_2 CD4 Tcells
Tcells@meta.data[Tcells$orig_clusts%in%c('C5_0', 'C5_1', 'C5_2', 'C5_3', 'C5_5', 'C5_8', 'C5_9'), 'Annotation_Level_6'] = 'EM_2 CD4 Tcells'
# 2.7. Follicular_1 CD4 Tcells
Tcells@meta.data[Tcells$orig_clusts=='C4_6', 'Annotation_Level_6'] = 'Follicular_1 CD4 Tcells'
# 2.8. Follicular_2 CD4 Tcells
Tcells@meta.data[Tcells$orig_clusts%in%c('C4_2', 'C4_7', 'C4_11', 'C4_13'), 'Annotation_Level_6'] = 'Follicular_2 CD4 Tcells'
# 2.9. IL22+ CD4 Tcells
Tcells@meta.data[Tcells$orig_clusts=='C5_6', 'Annotation_Level_6'] = 'IL22+ CD4 Tcells'
# 2.10. IL17A+ IL17F- IL22- CD4 Tcells
Tcells@meta.data[Tcells$orig_clusts%in%c('C5_4', 'C5_7'), 'Annotation_Level_6'] = 'IL17A+ IL17F- IL22- CD4 Tcells'
# 2.11. IL17A+ IL17F+ IL22- CD4 Tcells
Tcells@meta.data[Tcells$orig_clusts=='C5_10', 'Annotation_Level_6'] = 'IL17A+ IL17F+ IL22- CD4 Tcells'
# 2.12. IL17A+ IL17F- IL22+ CD4 Tcells
Tcells@meta.data[Tcells$orig_clusts=='C5_11', 'Annotation_Level_6'] = 'IL17A+ IL17F- IL22+ CD4 Tcells'
# 2.13. IL17A- IL17F+ IL22- CD4 Tcells
Tcells@meta.data[Tcells$orig_clusts=='C5_12', 'Annotation_Level_6'] = 'IL17A- IL17F+ IL22- CD4 Tcells'
# 2.14. Proliferative CD8 Tcells
Tcells@meta.data[Tcells$orig_clusts%in%c('C6_1', 'C6_2'), 'Annotation_Level_6'] = 'Proliferative CD8 Tcells'
# 2.15. Naive CD8 Tcells
Tcells@meta.data[Tcells$orig_clusts=='C1_15', 'Annotation_Level_6'] = 'Naive CD8 Tcells'
# 2.16. Memory CD160+ CD8 Tcells
Tcells@meta.data[Tcells$orig_clusts=='C1_0', 'Annotation_Level_6'] = 'Memory CD160+ CD8 Tcells'
# 2.17. EM CD8 Tcells
Tcells@meta.data[Tcells$orig_clusts%in%c('C1_3', 'C1_4', 'C1_5', 'C1_11', 'C1_13', 'C1_14', 'C1_16'), 'Annotation_Level_6'] = 'EM CD8 Tcells'
# 2.18. TRM CD8 Tcells
Tcells@meta.data[Tcells$orig_clusts%in%c('C1_7', 'C4_3', 'C4_10'), 'Annotation_Level_6'] = 'TRM CD8 Tcells'
# 2.19. CXCL13+ CD8 Tcells
Tcells@meta.data[Tcells$orig_clusts%in%c('C4_0', 'C4_1', 'C4_4', 'C4_5', 'C4_8', 'C4_9', 'C4_12'), 'Annotation_Level_6'] = 'CXCL13+ CD8 Tcells'
# 2.20. Cytotoxic CD8 Tcells
Tcells@meta.data[Tcells$orig_clusts=='C1_1', 'Annotation_Level_6'] = 'Cytotoxic CD8 Tcells'
# 2.20. Double-Negative CD8 Tcells
Tcells@meta.data[Tcells$orig_clusts=='C1_2', 'Annotation_Level_6'] = 'Double-Negative Tcells'
# 2.21. gdTcells
Tcells@meta.data[Tcells$orig_clusts%in%c('C3_2', 'C3_4', 'C3_5', 'C3_6', 'C3_10', 'C3_11', 'C3_12', 'C3_15'), 'Annotation_Level_6'] = 'gdTcells'
# 2.22. NKT_1 cells
Tcells@meta.data[Tcells$orig_clusts%in%c('C3_9', 'C3_13', 'C3_14'), 'Annotation_Level_6'] = 'NKT_1 cells'
# 2.23. NKT_2 cells
Tcells@meta.data[Tcells$orig_clusts%in%c('C3_8', 'C3_16', 'C3_17'), 'Annotation_Level_6'] = 'NKT_2 cells'
# 2.24. NKT_3 cells
Tcells@meta.data[Tcells$orig_clusts=='C1_10', 'Annotation_Level_6'] = 'NKT_3 cells'
# 2.25. Cytotoxic CD8aa IELs
Tcells@meta.data[Tcells$orig_clusts=='C1_9', 'Annotation_Level_6'] = 'Cytotoxic CD8aa IELs'
# 2.26. Non-Cytotoxic CD8aa IELs
Tcells@meta.data[Tcells$orig_clusts%in%c('C1_6', 'C1_8', 'C1_12'), 'Annotation_Level_6'] = 'Non-Cytotoxic CD8aa IELs'
# 2.27. NK cells
Tcells@meta.data[Tcells$orig_clusts%in%c('C3_0', 'C3_1', 'C3_3'), 'Annotation_Level_6'] = 'NK cells'
# 2.27. LTi cells
Tcells@meta.data[Tcells$orig_clusts=='C3_7', 'Annotation_Level_6'] = 'LTi cells'
# 2.28. Order the levels:
Tcells$Annotation_Level_6 = factor(Tcells$Annotation_Level_6,
                                   levels=c('Proliferative CD4 Tcells', 'Regulatory CD4 Tcells', 'Naive CD4 Tcells', 'CM CD4 Tcells',
                                            'EM_1 CD4 Tcells', 'EM_2 CD4 Tcells', 'Follicular_1 CD4 Tcells', 'Follicular_2 CD4 Tcells',
                                            'IL22+ CD4 Tcells', 'IL17A+ IL17F- IL22- CD4 Tcells', 'IL17A+ IL17F+ IL22- CD4 Tcells',
                                            'IL17A+ IL17F- IL22+ CD4 Tcells', 'IL17A- IL17F+ IL22- CD4 Tcells', 'Proliferative CD8 Tcells',
                                            'Naive CD8 Tcells', 'Memory CD160+ CD8 Tcells', 'EM CD8 Tcells', 'TRM CD8 Tcells',
                                            'CXCL13+ CD8 Tcells', 'Cytotoxic CD8 Tcells', 'Double-Negative Tcells', 'gdTcells',
                                            'NKT_1 cells', 'NKT_2 cells', 'NKT_3 cells', 'Cytotoxic CD8aa IELs', 'Non-Cytotoxic CD8aa IELs',
                                            'NK cells', 'LTi cells'))


# 3. Annotate cells - Level 5
Tcells[['Annotation_Level_5']] = as.character(Tcells$Annotation_Level_6)
# 3.1. EM CD4 Tcells
Tcells@meta.data[Tcells$Annotation_Level_6%in%c('EM_1 CD4 Tcells', 'EM_2 CD4 Tcells'), 'Annotation_Level_5'] = 'EM CD4 Tcells'
# 3.2. Follicular CD4 Tcells
Tcells@meta.data[Tcells$Annotation_Level_6%in%c('Follicular_1 CD4 Tcells', 'Follicular_2 CD4 Tcells'),
                 'Annotation_Level_5'] = 'Follicular CD4 Tcells'
# 3.3. IL17+ CD4 Tcells
Tcells@meta.data[Tcells$Annotation_Level_6%in%c('IL17A+ IL17F- IL22- CD4 Tcells', 'IL17A+ IL17F+ IL22- CD4 Tcells',
                                                'IL17A+ IL17F- IL22+ CD4 Tcells', 'IL17A- IL17F+ IL22- CD4 Tcells'),
                 'Annotation_Level_5'] = 'IL17+ CD4 Tcells'
# 3.4. NKT cells
Tcells@meta.data[Tcells$Annotation_Level_6%in%c('NKT_1 cells', 'NKT_2 cells', 'NKT_3 cells'), 'Annotation_Level_5'] = 'NKT cells'
# 3.5. CD8aa IELs
Tcells@meta.data[Tcells$Annotation_Level_6%in%c('Cytotoxic CD8aa IELs', 'Non-Cytotoxic CD8aa IELs'), 'Annotation_Level_5'] = 'CD8aa IELs'
# 3.6. Order the levels:
Tcells$Annotation_Level_5 = factor(Tcells$Annotation_Level_5,
                                   levels=c('Proliferative CD4 Tcells', 'Regulatory CD4 Tcells', 'Naive CD4 Tcells', 'CM CD4 Tcells',
                                            'EM CD4 Tcells', 'Follicular CD4 Tcells', 'IL22+ CD4 Tcells', 'IL17+ CD4 Tcells',
                                            'Proliferative CD8 Tcells', 'Naive CD8 Tcells', 'Memory CD160+ CD8 Tcells', 'EM CD8 Tcells',
                                            'TRM CD8 Tcells', 'CXCL13+ CD8 Tcells', 'Cytotoxic CD8 Tcells', 'Double-Negative Tcells',
                                            'gdTcells', 'NKT cells', 'CD8aa IELs', 'NK cells', 'LTi cells'))


# 4. Annotate cells - Level 4
Tcells[['Annotation_Level_4']] = as.character(Tcells$Annotation_Level_5)
# 4.1. Memory CD4 Tcells
Tcells@meta.data[Tcells$Annotation_Level_5%in%c('EM CD4 Tcells', 'CM CD4 Tcells'), 'Annotation_Level_4'] = 'Memory CD4 Tcells'
# 4.2. Memory CD8 Tcells
Tcells@meta.data[Tcells$Annotation_Level_5%in%c('EM CD8 Tcells', 'Memory CD160+ CD8 Tcells', 'TRM CD8 Tcells'),
                 'Annotation_Level_4'] = 'Memory CD8 Tcells'
# 4.3. Order the levels:
Tcells$Annotation_Level_4 = factor(Tcells$Annotation_Level_4,
                                   levels=c('Proliferative CD4 Tcells', 'Regulatory CD4 Tcells', 'Naive CD4 Tcells', 'Memory CD4 Tcells',
                                            'Follicular CD4 Tcells', 'IL22+ CD4 Tcells', 'IL17+ CD4 Tcells', 'Proliferative CD8 Tcells',
                                            'Naive CD8 Tcells', 'Memory CD8 Tcells', 'CXCL13+ CD8 Tcells', 'Cytotoxic CD8 Tcells',
                                            'Double-Negative Tcells', 'gdTcells', 'NKT cells', 'CD8aa IELs', 'NK cells', 'LTi cells'))


# 5. Annotate cells - Level 3
Tcells[['Annotation_Level_3']] = as.character(Tcells$Annotation_Level_4)
# 5.1. CD4 Tcells
Tcells@meta.data[Tcells$Annotation_Level_4%in%c('Proliferative CD4 Tcells', 'Regulatory CD4 Tcells', 'Naive CD4 Tcells',
                                                'Memory CD4 Tcells', 'Follicular CD4 Tcells', 'IL22+ CD4 Tcells', 'IL17+ CD4 Tcells'),
                 'Annotation_Level_3'] = 'CD4 Tcells'
# 5.2. CD8 Tcells
Tcells@meta.data[Tcells$Annotation_Level_4%in%c('Proliferative CD8 Tcells', 'Naive CD8 Tcells', 'Memory CD8 Tcells', 
                                                'CXCL13+ CD8 Tcells', 'Cytotoxic CD8 Tcells'),
                 'Annotation_Level_3'] = 'CD8 Tcells'
# 5.3. Order the levels:
Tcells$Annotation_Level_3 = factor(Tcells$Annotation_Level_3,
                                   levels=c('CD4 Tcells', 'CD8 Tcells', 'Double-Negative Tcells', 'gdTcells', 'NKT cells', 'CD8aa IELs',
                                            'NK cells', 'LTi cells'))


# 6. Annotate cells - Level 2
Tcells[['Annotation_Level_2']] = as.character(Tcells$Annotation_Level_3)
# 6.1. Unconventional Tcells
Tcells@meta.data[Tcells$Annotation_Level_3%in%c('gdTcells', 'NKT cells', 'CD8aa IELs'), 'Annotation_Level_2'] = 'Unconventional Tcells'
# 6.2. ILCs
Tcells@meta.data[Tcells$Annotation_Level_3%in%c('NK cells', 'LTi cells'), 'Annotation_Level_2'] = 'ILCs'
# 6.3. Order the levels:
Tcells$Annotation_Level_2 = factor(Tcells$Annotation_Level_2,
                                   levels=c('CD4 Tcells', 'CD8 Tcells', 'Double-Negative Tcells', 'Unconventional Tcells', 'ILCs'))


# 7. Annotate cells - Level 1
Tcells[['Annotation_Level_1']] = as.character(Tcells$Annotation_Level_2)
# 7.1. abTcells
Tcells@meta.data[Tcells$Annotation_Level_2%in%c('CD4 Tcells', 'CD8 Tcells', 'Double-Negative Tcells'), 'Annotation_Level_1'] = 'abTcells'
# 7.2. Order the levels:
Tcells$Annotation_Level_1 = factor(Tcells$Annotation_Level_1,
                                   levels=c('abTcells', 'Unconventional Tcells', 'ILCs'))


# 8. Annotate cells - Level 0
Tcells[['Annotation_Level_0']] = factor(rep('Tcells', dim(Tcells@meta.data)[1]))


# 9. Save Seurat object
SeuratDisk::SaveH5Seurat(Tcells, paste(project_dir, '2_annotation/results_Tcells/datasets/Tcells_finalAnnots.h5Seurat', sep='/'))


# 10. Visualize cells annotated in the different levels
annotation_plots = list()
for(annotation_level in grep('Annotation_Level', names(Tcells[[]]), value=T)[1:6]){
  annotation_plots[[annotation_level]] = Seurat::DimPlot(Tcells, reduction="umap", group.by=annotation_level,
                                                         label=TRUE, label.size=3, pt.size=.5, repel=T) + ggplot2::theme_minimal()
}
patchwork::wrap_plots(annotation_plots, ncol=3)


# 11. Separate cells in CD4 Tcells, run umap and visualize (annotations, feature plots, violin plots)
# 11.1. Subset
CD4Tcells = subset(Tcells, cells=rownames(Tcells@meta.data)[Tcells$Annotation_Level_2=='CD4 Tcells'])
# 11.2. PCA
CD4Tcells = Seurat::RunPCA(CD4Tcells, assay='integrated')
# 11.3. Choose number of PCs to use
pct = CD4Tcells[["pca"]]@stdev /sum(CD4Tcells[["pca"]]@stdev) * 100
cumu = cumsum(pct)
co1 = which(cumu > 90 & pct < 5)[1]
co2 = sort(which((pct[1:length(pct) - 1] - pct[2:length(pct)]) > 0.1), decreasing = T)[1] + 1
elbow = min(co1, co2) # 6
# 11.4. Determine the K-nearest neighbor graph
CD4Tcells = Seurat::FindNeighbors(CD4Tcells, dims=1:elbow)
# 11.5. UMAP:
CD4Tcells = Seurat::RunUMAP(CD4Tcells, dims=1:elbow, reduction='pca')
invisible(gc())
# 11.6. Visualize annotations:
Seurat::DimPlot(CD4Tcells, reduction="umap", group.by='Annotation_Level_4', label=TRUE, label.size=3, pt.size=.5, repel=T) +
  ggplot2::theme_minimal()
Seurat::DimPlot(CD4Tcells, reduction="umap", group.by='Annotation_Level_5', label=TRUE, label.size=3, pt.size=.5, repel=T) +
  ggplot2::theme_minimal()
Seurat::DimPlot(CD4Tcells, reduction="umap", group.by='Annotation_Level_6', label=TRUE, label.size=3, pt.size=.5, repel=T) +
  ggplot2::theme_minimal()
# 11.7. Feature Plots:
feature_plots(CD4Tcells, c('rna_CD3D', 'rna_CD3E', 'rna_CD3G', 'rna_CD8A', 'rna_CD8B', 'rna_CD4', 'rna_TRAC', 'rna_TRBC1', 'rna_TRBC2',
                          'rna_TRDC', 'rna_TRGC1'),
              ncol=3, with_dimplot=TRUE, dimplot_group='Annotation_Level_6', point_size = .5)
feature_plots(CD4Tcells, c('rna_CCR7', 'rna_SELL', 'rna_PASK', 'rna_ZNF683'),
              ncol=3, with_dimplot=TRUE, dimplot_group='Annotation_Level_6', point_size = .5) # Naive, CM
feature_plots(CD4Tcells, c('rna_KLRB1', 'rna_S100A4', 'rna_PTPRCAP', 'rna_PFN1', 'rna_GZMA'),
              ncol=2, with_dimplot=TRUE, dimplot_group='Annotation_Level_6', point_size = .5) # EM
feature_plots(CD4Tcells, c('rna_IL2RA', 'rna_IL2RB', 'rna_FOXP3'),
              ncol=2, with_dimplot=TRUE, dimplot_group='Annotation_Level_6', point_size = .5) # Regulatory CD4
feature_plots(CD4Tcells, c('rna_STAT4', 'rna_IL12RB2', 'rna_IFNG', 'rna_GATA3', 'rna_STAT6', 'rna_IL4',
                           'rna_CXCL13', 'rna_IL17A', 'rna_IL17F', 'rna_IL22'),
              ncol=3, with_dimplot=TRUE, dimplot_group='Annotation_Level_6', point_size = .5) # Other Thelper
# 11.8. Violin Plots:
violin_plots(CD4Tcells, c('rna_CD3D', 'rna_CD3E', 'rna_CD3G', 'rna_CD8A', 'rna_CD8B', 'rna_CD4', 'rna_TRAC', 'rna_TRBC1', 'rna_TRBC2',
                           'rna_TRDC', 'rna_TRGC1'), group.by='Annotation_Level_6', ncol=4, with_dimplot=FALSE)
violin_plots(CD4Tcells, c('rna_CCR7', 'rna_SELL', 'rna_PASK', 'rna_ZNF683'),
             group.by='Annotation_Level_6', ncol=2, with_dimplot=FALSE) # Naive, CM
violin_plots(CD4Tcells, c('rna_KLRB1', 'rna_S100A4', 'rna_PTPRCAP', 'rna_PFN1', 'rna_GZMA'),
             group.by='Annotation_Level_6', ncol=3,with_dimplot=FALSE) # EM
violin_plots(CD4Tcells, c('rna_IL2RA', 'rna_IL2RB', 'rna_FOXP3'),  group.by='Annotation_Level_6', ncol=2, with_dimplot=FALSE) # Reg CD4
violin_plots(CD4Tcells, c('rna_STAT4', 'rna_IL12RB2', 'rna_IFNG', 'rna_GATA3', 'rna_STAT6', 'rna_IL4',
                          'rna_CXCL13', 'rna_IL17A', 'rna_IL17F', 'rna_IL22'),
             group.by='Annotation_Level_6', ncol=4, with_dimplot=FALSE) # Other Thelper


# 12. Separate cells in CD8 Tcells, run umap and visualize (annotations, feature plots, violin plots)
# 12.1. Subset
CD8Tcells = subset(Tcells, cells=rownames(Tcells@meta.data)[Tcells$Annotation_Level_2=='CD8 Tcells'])
# 12.2. PCA
CD8Tcells = Seurat::RunPCA(CD8Tcells, assay='integrated')
# 12.3. Choose number of PCs to use
pct = CD8Tcells[["pca"]]@stdev /sum(CD8Tcells[["pca"]]@stdev) * 100
cumu = cumsum(pct)
co1 = which(cumu > 90 & pct < 5)[1]
co2 = sort(which((pct[1:length(pct) - 1] - pct[2:length(pct)]) > 0.1), decreasing = T)[1] + 1
elbow = min(co1, co2) # 8
# 12.4. Determine the K-nearest neighbor graph
CD8Tcells = Seurat::FindNeighbors(CD8Tcells, dims=1:elbow)
# 12.5. UMAP:
CD8Tcells = Seurat::RunUMAP(CD8Tcells, dims=1:elbow, reduction='pca')
invisible(gc())
# 12.6. Visualize annotations:
Seurat::DimPlot(CD8Tcells, reduction="umap", group.by='Annotation_Level_4', label=TRUE, label.size=3, pt.size=.5, repel=T) +
  ggplot2::theme_minimal()
Seurat::DimPlot(CD8Tcells, reduction="umap", group.by='Annotation_Level_5', label=TRUE, label.size=3, pt.size=.5, repel=T) +
  ggplot2::theme_minimal()
Seurat::DimPlot(CD8Tcells, reduction="umap", group.by='Annotation_Level_6', label=TRUE, label.size=3, pt.size=.5, repel=T) +
  ggplot2::theme_minimal()
# 12.7. Feature Plots:
feature_plots(CD8Tcells, c('rna_CD3D', 'rna_CD3E', 'rna_CD3G', 'rna_CD8A', 'rna_CD8B', 'rna_CD4', 'rna_TRAC', 'rna_TRBC1', 'rna_TRBC2',
                           'rna_TRDC', 'rna_TRGC1'),
              ncol=3, with_dimplot=TRUE, dimplot_group='Annotation_Level_6', point_size = .5)
feature_plots(CD8Tcells, c('rna_CCR7', 'rna_SELL', 'rna_PASK', 'rna_ZNF683'),
              ncol=3, with_dimplot=TRUE, dimplot_group='Annotation_Level_6', point_size = .5) # Naive, CM
feature_plots(CD8Tcells, c('rna_KLRB1', 'rna_S100A4', 'rna_PTPRCAP', 'rna_PFN1', 'rna_GZMA'),
              ncol=2, with_dimplot=TRUE, dimplot_group='Annotation_Level_6', point_size = .5) # EM
feature_plots(CD8Tcells, c('rna_GZMA', 'rna_GZMB', 'rna_GZMH', 'rna_GZMK', 'rna_GZMM',
                           'rna_PRF1', 'rna_KLRB1', 'rna_GNLY', 'rna_NKG7', 'rna_IFNG'),
              ncol=3, with_dimplot=TRUE, dimplot_group='Annotation_Level_6', point_size = .5) # Cytotoxic
# 12.7. Violin Plots:
violin_plots(CD8Tcells, c('rna_CD3D', 'rna_CD3E', 'rna_CD3G', 'rna_CD8A', 'rna_CD8B', 'rna_CD4', 'rna_TRAC', 'rna_TRBC1', 'rna_TRBC2',
                          'rna_TRDC', 'rna_TRGC1'), ncol=4, with_dimplot=FALSE, group.by='Annotation_Level_6')
violin_plots(CD8Tcells, c('rna_CCR7', 'rna_SELL', 'rna_PASK', 'rna_ZNF683'),
             ncol=2, with_dimplot=FALSE, group.by='Annotation_Level_6') # Naive, CM
violin_plots(CD8Tcells, c('rna_KLRB1', 'rna_S100A4', 'rna_PTPRCAP', 'rna_PFN1', 'rna_GZMA'),
             ncol=3, with_dimplot=FALSE, group.by='Annotation_Level_6') # EM
violin_plots(CD8Tcells, c('rna_GZMA', 'rna_GZMB', 'rna_GZMH', 'rna_GZMK', 'rna_GZMM', 'rna_PRF1', 'rna_KLRB1',
                          'rna_GNLY', 'rna_NKG7', 'rna_IFNG', 'rna_MKI67', 'rna_PCNA'),
             ncol=4, with_dimplot=FALSE, group.by='Annotation_Level_6') # Cytotoxic


# 13. Unconventional Tcells + ILCs + DNs, run umap and visualize (annotations, feature plots, violin plots)
# 13.1. Subset
uncT_ILCs_cells = subset(Tcells, cells=rownames(Tcells@meta.data)[Tcells$Annotation_Level_2%in%c('Unconventional Tcells', 'ILCs',
                                                                                                 'Double-Negative Tcells')])
# 13.2. PCA
uncT_ILCs_cells = Seurat::RunPCA(uncT_ILCs_cells, assay='integrated')
# 13.3. Choose number of PCs to use
pct = uncT_ILCs_cells[["pca"]]@stdev /sum(uncT_ILCs_cells[["pca"]]@stdev) * 100
cumu = cumsum(pct)
co1 = which(cumu > 90 & pct < 5)[1]
co2 = sort(which((pct[1:length(pct) - 1] - pct[2:length(pct)]) > 0.1), decreasing = T)[1] + 1
elbow = min(co1, co2) # 12
# 13.4. Determine the K-nearest neighbor graph
uncT_ILCs_cells = Seurat::FindNeighbors(uncT_ILCs_cells, dims=1:elbow)
# 13.5. UMAP:
uncT_ILCs_cells = Seurat::RunUMAP(uncT_ILCs_cells, dims=1:elbow, reduction='pca')
invisible(gc())
# 13.6. Visualize annotations:
Seurat::DimPlot(uncT_ILCs_cells, reduction="umap", group.by='Annotation_Level_4', label=TRUE, label.size=3, pt.size=.5, repel=T) +
  ggplot2::theme_minimal()
Seurat::DimPlot(uncT_ILCs_cells, reduction="umap", group.by='Annotation_Level_5', label=TRUE, label.size=3, pt.size=.5, repel=T) +
  ggplot2::theme_minimal()
Seurat::DimPlot(uncT_ILCs_cells, reduction="umap", group.by='Annotation_Level_6', label=TRUE, label.size=3, pt.size=.5, repel=T) +
  ggplot2::theme_minimal()
# 13.7. Feature Plots:
feature_plots(uncT_ILCs_cells, c('rna_CD3D', 'rna_CD3E', 'rna_CD3G', 'rna_CD8A', 'rna_CD8B', 'rna_CD4', 'rna_TRAC', 'rna_TRBC2',
                           'rna_TRDC', 'rna_TRGC1', 'rna_TRGC2'),
              ncol=3, with_dimplot=TRUE, dimplot_group='Annotation_Level_6', point_size = .5)
feature_plots(uncT_ILCs_cells, c('rna_GZMA', 'rna_GZMB', 'rna_GZMH', 'rna_GZMK', 'rna_GZMM', 'rna_PRF1', 'rna_KLRB1',
                                 'rna_GNLY', 'rna_NKG7'),
              ncol=3, with_dimplot=TRUE, dimplot_group='Annotation_Level_6', point_size = .5) # Cytotoxic
feature_plots(uncT_ILCs_cells, c('rna_NCR1', 'rna_NCAM1', 'rna_TYROBP', 'rna_FGFBP2', 'rna_KLRD1', 'rna_KLRF1', 'rna_KLRB1', 'rna_CX3CR1',
                                 'rna_FCGR3A', 'rna_XCL1', 'rna_XCL2', 'rna_GNLY', 'rna_NKG7', 'rna_PRF1'),
              ncol=3, with_dimplot=TRUE, dimplot_group='Annotation_Level_6', point_size = .5) # NK like
feature_plots(uncT_ILCs_cells, c('rna_RORC', 'rna_LTA', 'rna_LTB'),
              ncol=2, with_dimplot=TRUE, dimplot_group='Annotation_Level_6', point_size = .5) # ILC (LTi)
# 13.7. Violin Plots:
violin_plots(uncT_ILCs_cells, c('rna_CD3D', 'rna_CD3E', 'rna_CD3G', 'rna_CD8A', 'rna_CD8B', 'rna_CD4', 'rna_TRAC', 'rna_TRBC2',
                                 'rna_TRDC', 'rna_TRGC1', 'rna_TRGC2'),
              ncol=4, with_dimplot=FALSE, group.by='Annotation_Level_6')
violin_plots(uncT_ILCs_cells, c('rna_GZMA', 'rna_GZMB', 'rna_GZMH', 'rna_GZMK', 'rna_GZMM', 'rna_PRF1', 'rna_KLRB1',
                                 'rna_GNLY', 'rna_NKG7'),
              ncol=3, with_dimplot=FALSE, group.by='Annotation_Level_6') # Cytotoxic
violin_plots(uncT_ILCs_cells, c('rna_NCR1', 'rna_NCAM1', 'rna_TYROBP', 'rna_FGFBP2', 'rna_KLRD1', 'rna_KLRF1', 'rna_KLRB1', 'rna_CX3CR1',
                                 'rna_FCGR3A', 'rna_XCL1', 'rna_XCL2', 'rna_GNLY', 'rna_NKG7', 'rna_PRF1'),
              ncol=4, with_dimplot=FALSE, group.by='Annotation_Level_6') # NK like
violin_plots(uncT_ILCs_cells, c('rna_RORC', 'rna_LTA', 'rna_LTB'),
              ncol=2, with_dimplot=FALSE, group.by='Annotation_Level_6') # ILC (LTi)


# 14. Violin plot with only NK, NKT and CD8 subsets:
violin_plots(subset(Tcells, cells=rownames(Tcells@meta.data)[Tcells$Annotation_Level_3%in%c('CD8 Tcells', 'NK cells', 'NKT cells')]),
             c('rna_NCR1', 'rna_NCAM1', 'rna_TYROBP', 'rna_FGFBP2', 'rna_KLRD1', 'rna_KLRF1', 'rna_KLRB1', 'rna_CX3CR1', 'rna_FCGR3A',
               'rna_XCL1', 'rna_XCL2', 'rna_GNLY', 'rna_NKG7', 'rna_PRF1'),
             ncol=4, with_dimplot=FALSE, group.by='Annotation_Level_6') # NK like


# 15. Get Markers between the final annotations:
# 15.1. Annotation Level 6
Seurat::Idents(Tcells) = 'Annotation_Level_6'
Tcells_Level_6_markers = Seurat::FindAllMarkers(Tcells, assay='RNA', slot='data', logfc.threshold=.8, min.pct = 0.3, only.pos=TRUE)
View(Tcells_Level_6_markers)
write.csv(Tcells_Level_6_markers, paste(project_dir, '2_annotation/results_Tcells/markers/markers_Level_6.csv', sep='/'))
invisible(gc())

# 15.2. Annotation Level 5
Seurat::Idents(Tcells) = 'Annotation_Level_5'
Tcells_Level_5_markers = Seurat::FindAllMarkers(Tcells, assay='RNA', slot='data', logfc.threshold=.8, min.pct = 0.3, only.pos=TRUE)
View(Tcells_Level_5_markers)
write.csv(Tcells_Level_5_markers, paste(project_dir, '2_annotation/results_Tcells/markers/markers_Level_5.csv', sep='/'))
invisible(gc())

# 15.3. Annotation Level 4
Seurat::Idents(Tcells) = 'Annotation_Level_4'
Tcells_Level_4_markers = Seurat::FindAllMarkers(Tcells, assay='RNA', slot='data', logfc.threshold=.8, min.pct = 0.3, only.pos=TRUE)
View(Tcells_Level_4_markers)
write.csv(Tcells_Level_4_markers, paste(project_dir, '2_annotation/results_Tcells/markers/markers_Level_4.csv', sep='/'))
invisible(gc())

# 15.4. Annotation Level 3
Seurat::Idents(Tcells) = 'Annotation_Level_3'
Tcells_Level_3_markers = Seurat::FindAllMarkers(Tcells, assay='RNA', slot='data', logfc.threshold=.8, min.pct = 0.3, only.pos=TRUE)
View(Tcells_Level_3_markers)
write.csv(Tcells_Level_3_markers, paste(project_dir, '2_annotation/results_Tcells/markers/markers_Level_3.csv', sep='/'))
invisible(gc())

# 15.5. Annotation Level 2
Seurat::Idents(Tcells) = 'Annotation_Level_2'
Tcells_Level_2_markers = Seurat::FindAllMarkers(Tcells, assay='RNA', slot='data', logfc.threshold=.8, min.pct = 0.3, only.pos=TRUE)
View(Tcells_Level_2_markers)
write.csv(Tcells_Level_2_markers, paste(project_dir, '2_annotation/results_Tcells/markers/markers_Level_2.csv', sep='/'))
invisible(gc())

# 15.6. Annotation Level 1
Seurat::Idents(Tcells) = 'Annotation_Level_1'
Tcells_Level_1_markers = Seurat::FindAllMarkers(Tcells, assay='RNA', slot='data', logfc.threshold=.8, min.pct = 0.3, only.pos=TRUE)
View(Tcells_Level_1_markers)
write.csv(Tcells_Level_1_markers, paste(project_dir, '2_annotation/results_Tcells/markers/markers_Level_1.csv', sep='/'))
invisible(gc())






























