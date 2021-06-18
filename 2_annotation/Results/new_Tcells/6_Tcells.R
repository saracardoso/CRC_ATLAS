project_dir = '/home/scardoso/Documents/PhD/CRC_ATLAS'
source(paste(project_dir, 'utils/modified_plots.R', sep='/'))

Tcells = SeuratDisk::LoadH5Seurat('/home/scardoso/Documents/PhD/CRC_ATLAS/2_annotation/Results/new_Tcells/Tcells.h5Seurat')





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
SeuratDisk::SaveH5Seurat(Tcells, paste(project_dir, '2_annotation/Results/new_Tcells/Tcells_clustered.h5Seurat', sep='/'))



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
saveRDS(SIGMA_all_res02, paste(project_dir, '2_annotation/Results/new_Tcells/SIGMA_all_res02.Rdata', sep='/'))

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
#Tcells = SeuratDisk::LoadH5Seurat(paste(project_dir, '2_annotation/Results/new_Tcells/Tcells_clustered.h5Seurat', sep='/'))


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
SeuratDisk::SaveH5Seurat(Tcells_0, paste(project_dir, '2_annotation/Results/new_Tcells/Tcells_0.h5Seurat', sep='/'))


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
write.csv(Tcells_0_markers, paste(project_dir, '2_annotation/Results/new_Tcells/markers_C0_res02.csv', sep='/'))
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
write.csv(Tcells_0_markers_top20, paste(project_dir, '2_annotation/Results/new_Tcells/markers_C0_res02_top20.csv', sep='/'))


# 4. Annotation:
# Sub-cluster 0        --> Effector Memory (EM) CD4
# Sub-cluster 1        --> Naive CD4
# Sub-cluster 2        --> EM CD4
# Sub-cluster 3        --> Central Memory (CM) CD4


# 5. Check previous SIGMA result and colour cluster 0 by new sub-clusters
SIGMA_all_res02 = readRDS(paste(project_dir, '2_annotation/Results/new_Tcells/SIGMA_all_res02.Rdata', sep='/'))
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
SeuratDisk::SaveH5Seurat(Tcells_0, paste(project_dir, '2_annotation/Results/new_Tcells/Tcells_0.h5Seurat', sep='/'), overwrite=TRUE)





# -----
# - Sub-cluster cluster 1 from resolution 0.2 separately
# -----
#Tcells = SeuratDisk::LoadH5Seurat(paste(project_dir, '2_annotation/Results/new_Tcells/Tcells_integratedClusters.h5Seurat', sep='/'))


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
SeuratDisk::SaveH5Seurat(Tcells_1, paste(project_dir, '2_annotation/Results/new_Tcells/Tcells_1.h5Seurat', sep='/'))


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
write.csv(Tcells_1_markers, paste(project_dir, '2_annotation/Results/new_Tcells/markers_C1_res09.csv', sep='/'))
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
write.csv(Tcells_1_markers_top20, paste(project_dir, '2_annotation/Results/new_Tcells/markers_C1_res09_top20.csv', sep='/'))


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
SIGMA_all_res02 = readRDS(paste(project_dir, '2_annotation/Results/new_Tcells/SIGMA_all_res02.Rdata', sep='/'))
cluster_1_cells = rownames(Tcells_1@meta.data)
subclusters_1 = Tcells_1$Final_Annotation[intersect(colnames(SIGMA_all_res02$input_parameters$expr), cluster_1_cells)]
SIGMA::plot_singular_vectors(SIGMA_all_res02, '1', colour=as.character(subclusters_1))
# 4.3. Visualize annotations:
Seurat::DimPlot(Tcells_1, group.by='Final_Annotation', pt.size=0.5, label=F) +
  ggplot2::ggtitle('') + ggplot2::theme_minimal()# + Seurat::NoLegend()
# 4.4. Save changes:
SeuratDisk::SaveH5Seurat(Tcells_1, paste(project_dir, '2_annotation/Results/new_Tcells/Tcells_1.h5Seurat', sep='/'), overwrite=TRUE)





# -----
# - Sub-cluster cluster 2 from resolution 0.1 separately
# -----
#Tcells = SeuratDisk::LoadH5Seurat(paste(project_dir, '2_annotation/Results/new_Tcells/Tcells_integratedClusters.h5Seurat', sep='/'))


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
SeuratDisk::SaveH5Seurat(Tcells_2, paste(project_dir, '2_annotation/Results/new_Tcells/Tcells_2.h5Seurat', sep='/'), overwrite=TRUE)


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
write.csv(Tcells_2_markers, paste(project_dir, '2_annotation/Results/new_Tcells/markers_C2_res02.csv', sep='/'))
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
write.csv(Tcells_2_markers_top20, paste(project_dir, '2_annotation/Results/new_Tcells/markers_C2_res02_top20.csv', sep='/'))


# 4. Annotate sub-clusters:
Tcells_2[['Final_Annotation']] = rep('Regulatory Tcells', length(Tcells_2$integrated_snn_res.0.2))
# 4.1. Visualize annotations:
Seurat::DimPlot(Tcells_2, group.by='Final_Annotation', pt.size=0.5, label=T, label.size=4) +
  ggplot2::ggtitle('') + ggplot2::theme_minimal() + Seurat::NoLegend()
# 4.2. Save changes:
SeuratDisk::SaveH5Seurat(Tcells_2, paste(project_dir, '2_annotation/Results/new_Tcells/Tcells_2.h5Seurat', sep='/'), overwrite=TRUE)





# -----
# - Sub-cluster cluster 3 from resolution 0.1 separately
# -----
#Tcells = SeuratDisk::LoadH5Seurat(paste(project_dir, '2_annotation/Results/new_Tcells/Tcells_integratedClusters.h5Seurat', sep='/'))


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
SeuratDisk::SaveH5Seurat(Tcells_3, paste(project_dir, '2_annotation/Results/new_Tcells/Tcells_3.h5Seurat', sep='/'), overwrite=TRUE)


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
write.csv(Tcells_3_markers, paste(project_dir, '2_annotation/Results/new_Tcells/markers_C3_res07.csv', sep='/'))
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
write.csv(Tcells_3_markers_top20, paste(project_dir, '2_annotation/Results/new_Tcells/markers_C3_res08_top20.csv', sep='/'))

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
SIGMA_all_res02 = readRDS(paste(project_dir, '2_annotation/Results/new_Tcells/SIGMA_all_res02.Rdata', sep='/'))
cluster_3_cells = rownames(Tcells_3@meta.data)
subclusters_3 = Tcells_3$Final_Annotation[intersect(colnames(SIGMA_all_res02$input_parameters$expr), cluster_3_cells)]
SIGMA::plot_singular_vectors(SIGMA_all_res02, '3', colour=as.character(subclusters_3))
# 4.2. Visualize annotations:
Seurat::DimPlot(Tcells_3, group.by='Final_Annotation', pt.size=0.8, label=F) +#T, label.size=4) +
  ggplot2::ggtitle('') + ggplot2::theme_minimal()# + Seurat::NoLegend()
# 4.3. Save changes:
SeuratDisk::SaveH5Seurat(Tcells_3, paste(project_dir, '2_annotation/Results/new_Tcells/Tcells_3.h5Seurat', sep='/'), overwrite=TRUE)





















# -----
# - Sub-cluster cluster 4 from resolution 0.1 separately
# -----
#Tcells = SeuratDisk::LoadH5Seurat(paste(project_dir, '2_annotation/Results/new_Tcells/Tcells_integratedClusters.h5Seurat', sep='/'))


# 1. Sub-cluster:
Tcells_4 = subset(Tcells, cells=rownames(Tcells@meta.data)[Tcells$integrated_snn_res.0.1=='4'])
Seurat::DimPlot(Tcells, group.by = 'integrated_snn_res.0.1') | Seurat::DimPlot(Tcells_4, group.by = 'integrated_snn_res.0.1')
# 1.1. PCA
Tcells_4 = Seurat::RunPCA(Tcells_4, assay='integrated')
# 1.2. Determine the K-nearest neighbor graph
Tcells_4 = Seurat::FindNeighbors(Tcells_4, dims=1:50)
# 1.3. Find clusters:
Tcells_4 = Seurat::FindClusters(Tcells_4, resolution=seq(0.1, 1, by=.1))
# 1.4. UMAP:
Tcells_4 = Seurat::RunUMAP(Tcells_4, dims=1:50, reduction='pca')
invisible(gc())
# 1.5. Save object
SeuratDisk::SaveH5Seurat(Tcells_4, paste(project_dir, '2_annotation/Results/Tcells/Tcells_4.h5Seurat', sep='/'))


# 2. Choose best resolution using clustree:
library(ggraph)
clust_tree = clustree::clustree(Tcells_4, prefix='integrated_snn_res.')
clust_tree
# 2.1. Visualize different resolutions (first resolution was chosen to be the best for now):
clusters_01_0 = Seurat::DimPlot(Tcells_4, reduction="umap", group.by='integrated_snn_res.0.1', label=TRUE, label.size=6, pt.size=.5)
clusters_03_0 = Seurat::DimPlot(Tcells_4, reduction="umap", group.by='integrated_snn_res.0.3', label=TRUE, label.size=6, pt.size=.5)
clusters_06_0 = Seurat::DimPlot(Tcells_4, reduction="umap", group.by='integrated_snn_res.0.6', label=TRUE, label.size=6, pt.size=.5)
clusters_07_0 = Seurat::DimPlot(Tcells_4, reduction="umap", group.by='integrated_snn_res.0.7', label=TRUE, label.size=6, pt.size=.5)
((clusters_01_0 | clusters_03_0) / (clusters_06_0 | clusters_07_0)) | clust_tree
# 2.2. Feature Plots to assess known gene markers:
feature_plots(Tcells_4, c('rna_CCR7', 'rna_SELL', 'rna_PASK', 'rna_CXCL13'), point_size = 1,
              ncol=3, with_dimplot=TRUE, dimplot_group='integrated_snn_res.0.6')
# 2.3. Violin Plots to assess known markers:
violin_plots(Tcells_4, c('rna_CCR7', 'rna_SELL', 'rna_PASK', 'rna_CXCL13'),
             ncol=3, with_dimplot=TRUE, group.by='integrated_snn_res.0.6')
# 2.4. Cluster 4 will not be sub-clustered and will be considered as follicular CD4 Tcells

# Check markers between cluster 4 and all others, just to be sure it is follicular CD4 Tcells
# 3.1. Calculate all markers
Seurat::Idents(Tcells) = 'integrated_snn_res.0.1'
Tcells_4_markers = Seurat::FindMarkers(Tcells, ident.1='4', assay='RNA', slot='data', logfc.threshold=.5, only.pos=TRUE)
View(Tcells_4_markers)
# 3.2. Get top 20 markers for each sub-cluster
Tcells_4_markers_top20 = c()
n_markers = dim(Tcells_4_markers)[1]
order_by_avg_log2FC = order(Tcells_4_markers$avg_log2FC, decreasing=T)
if(n_markers > 20) top20_clust_markers = Tcells_4_markers[order_by_avg_log2FC[1:20],]
else top20_clust_markers = Tcells_4_markers[order_by_avg_log2FC,]
Tcells_4_markers_top20 = rbind(Tcells_4_markers_top20, top20_clust_markers)
View(Tcells_4_markers_top20)

# 4. Annotate cells:
# 4.1. Insert the annotations in the metadata:
Tcells_4[['Final_Annotation']] = rep('Follicular CD4 Tcells', length(Tcells_4$integrated_snn_res.0.6))
# 4.2. Visualize annotations:
Seurat::DimPlot(Tcells_4, group.by='Final_Annotation', pt.size=0.5, label=T, label.size=4) +
  ggplot2::ggtitle('') + ggplot2::theme_minimal() + Seurat::NoLegend()
# 4.3. Save changes:
SeuratDisk::SaveH5Seurat(Tcells_4, paste(project_dir, '2_annotation/Results/Tcells/Tcells_4.h5Seurat', sep='/'), overwrite=TRUE)





# -----
# - Sub-cluster cluster 6 from resolution 0.1 separately
# -----
Tcells = SeuratDisk::LoadH5Seurat(paste(project_dir, '2_annotation/data/Tcells_integratedClusters.h5Seurat', sep='/'))


# 1. Sub-cluster:
Tcells_6 = subset(Tcells, cells=rownames(Tcells@meta.data)[Tcells$integrated_snn_res.0.1=='6'])
Seurat::DimPlot(Tcells, group.by = 'integrated_snn_res.0.1') | Seurat::DimPlot(Tcells_6, group.by = 'integrated_snn_res.0.1')
# 1.1. PCA
Tcells_6 = Seurat::RunPCA(Tcells_6, assay='integrated')
# 1.2. Determine the K-nearest neighbor graph
Tcells_6 = Seurat::FindNeighbors(Tcells_6, dims=1:50)
# 1.3. Find clusters:
Tcells_6 = Seurat::FindClusters(Tcells_6, resolution=seq(0.1, 1, by=.1))
# 1.4. UMAP:
Tcells_6 = Seurat::RunUMAP(Tcells_6, dims=1:50, reduction='pca')
invisible(gc())
# 1.5. Save object
SeuratDisk::SaveH5Seurat(Tcells_6, paste(project_dir, '2_annotation/Results/Tcells/Tcells_6.h5Seurat', sep='/'))


# 2. Choose best resolution using clustree:
library(ggraph)
clust_tree = clustree::clustree(Tcells_6, prefix='integrated_snn_res.')
clust_tree
# 2.1. Visualize different resolutions (first resolution was chosen to be the best for now):
clusters_01_0 = Seurat::DimPlot(Tcells_6, reduction="umap", group.by='integrated_snn_res.0.1', label=TRUE, label.size=6, pt.size=.5)
clusters_07_0 = Seurat::DimPlot(Tcells_6, reduction="umap", group.by='integrated_snn_res.0.7', label=TRUE, label.size=6, pt.size=.5)
clusters_08_0 = Seurat::DimPlot(Tcells_6, reduction="umap", group.by='integrated_snn_res.0.8', label=TRUE, label.size=6, pt.size=.5)
clusters_09_0 = Seurat::DimPlot(Tcells_6, reduction="umap", group.by='integrated_snn_res.0.9', label=TRUE, label.size=6, pt.size=.5)
((clusters_01_0 | clusters_07_0) / (clusters_08_0 | clusters_09_0)) | clust_tree
# 2.2. Feature Plots to assess known gene markers:
feature_plots(Tcells_6, c('rna_CD3E', 'rna_CD3D', 'rna_CD3G', 'rna_CD3G', 'rna_KLRB1', 'rna_GNLY', 'rna_NKG7', 'rna_KIT'),
              ncol=3, with_dimplot=TRUE, dimplot_group='integrated_snn_res.0.9', point_size = 1)
feature_plots(Tcells_6, c('rna_TNFSF10', 'rna_FASLG', 'rna_TNF', 'rna_CSF2', 'rna_TBX21', 'rna_RORC'),
              ncol=3, with_dimplot=TRUE, dimplot_group='integrated_snn_res.0.9', point_size = 1)
feature_plots(Tcells_6, c('rna_IL5', 'rna_IL13', 'rna_IL17A', 'rna_IL22'),
              ncol=3, with_dimplot=TRUE, dimplot_group='integrated_snn_res.0.9', point_size = 1)
# 2.3. Violin Plots to assess known markers:
violin_plots(Tcells_6, c('rna_CD3E', 'rna_CD3D', 'rna_CD3G', 'rna_CD3G', 'rna_KLRB1'),
             ncol=3, with_dimplot=TRUE, group.by='integrated_snn_res.0.4')
# 2.4. Cluster 6 will not be sub-clustered and might be composed of NK cells

# Check markers between cluster 6 and all others, just to understand what this cells might be overall
# 3.1. Calculate all markers
Seurat::Idents(Tcells) = 'integrated_snn_res.0.1'
Tcells_6_markers_overall = Seurat::FindMarkers(Tcells, ident.1='6', assay='RNA', slot='data', logfc.threshold=.8, only.pos=TRUE)
View(Tcells_6_markers_overall)
# 3.2. Get top 20 markers for each sub-cluster
Tcells_6_markers_overall_top20 = c()
n_markers = dim(Tcells_6_markers_overall)[1]
order_by_avg_log2FC = order(Tcells_6_markers_overall$avg_log2FC, decreasing=T)
if(n_markers > 20) top20_clust_markers = Tcells_6_markers_overall[order_by_avg_log2FC[1:20],]
else top20_clust_markers = Tcells_6_markers_overall[order_by_avg_log2FC,]
Tcells_6_markers_overall_top20 = rbind(Tcells_6_markers_overall_top20, top20_clust_markers)
View(Tcells_6_markers_overall_top20)

# 4. Annotate cells:
# 4.1. Insert the annotations in the metadata:
Tcells_6[['Final_Annotation']] = rep('NK cells', length(Tcells_6$integrated_snn_res.0.1))
# 4.2. Visualize annotations:
Seurat::DimPlot(Tcells_6, group.by='Final_Annotation', pt.size=0.5, label=T, label.size=4) +
  ggplot2::ggtitle('') + ggplot2::theme_minimal() + Seurat::NoLegend()
# 4.3. Save changes:
SeuratDisk::SaveH5Seurat(Tcells_6, paste(project_dir, '2_annotation/Results/Tcells/Tcells_6.h5Seurat', sep='/'), overwrite=TRUE)





# -----
# - Sub-cluster cluster 5 from resolution 0.1 separately
# -----
Tcells = SeuratDisk::LoadH5Seurat(paste(project_dir, '2_annotation/data/Tcells_integratedClusters.h5Seurat', sep='/'))


# 1. Sub-cluster:
Tcells_5 = subset(Tcells, cells=rownames(Tcells@meta.data)[Tcells$integrated_snn_res.0.1=='5'])
Seurat::DimPlot(Tcells, group.by = 'integrated_snn_res.0.1') | Seurat::DimPlot(Tcells_5, group.by = 'integrated_snn_res.0.1')
# 1.1. PCA
Tcells_5 = Seurat::RunPCA(Tcells_5, assay='integrated')
# 1.2. Determine the K-nearest neighbor graph
Tcells_5 = Seurat::FindNeighbors(Tcells_5, dims=1:50)
# 1.3. Find clusters:
Tcells_5 = Seurat::FindClusters(Tcells_5, resolution=seq(0.1, 1, by=.1))
# 1.4. UMAP:
Tcells_5 = Seurat::RunUMAP(Tcells_5, dims=1:50, reduction='pca')
invisible(gc())
# 1.5. Save object
SeuratDisk::SaveH5Seurat(Tcells_5, paste(project_dir, '2_annotation/Results/Tcells/Tcells_5.h5Seurat', sep='/'))


# 2. Choose best resolution using clustree:
library(ggraph)
clust_tree = clustree::clustree(Tcells_5, prefix='integrated_snn_res.')
clust_tree
# 2.1. Visualize different resolutions (first resolution was chosen to be the best for now):
clusters_02_0 = Seurat::DimPlot(Tcells_5, reduction="umap", group.by='integrated_snn_res.0.2', label=TRUE, label.size=6, pt.size=.5)
clusters_06_0 = Seurat::DimPlot(Tcells_5, reduction="umap", group.by='integrated_snn_res.0.6', label=TRUE, label.size=6, pt.size=.5)
clusters_07_0 = Seurat::DimPlot(Tcells_5, reduction="umap", group.by='integrated_snn_res.0.7', label=TRUE, label.size=6, pt.size=.5)
clusters_08_0 = Seurat::DimPlot(Tcells_5, reduction="umap", group.by='integrated_snn_res.0.8', label=TRUE, label.size=6, pt.size=.5)
((clusters_02_0 | clusters_06_0) / (clusters_07_0 | clusters_08_0)) | clust_tree
# 2.2. Feature Plots to assess known gene markers:
feature_plots(Tcells_5, c('rna_CD3E', 'rna_CD3D', 'rna_CD3G', 'rna_CD3G', 'rna_CD8A', 'rna_CD4', 'rna_TRAC', 'rna_TRBC1', 'rna_TRBC2'),
              ncol=3, with_dimplot=TRUE, dimplot_group='integrated_snn_res.0.8', point_size = 1)
feature_plots(Tcells_5, c('rna_KLRB1', 'rna_GNLY', 'rna_NKG7', 'rna_GZMA', 'rna_GZMB', 'rna_GZMH', 'rna_GZMK', 'rna_GZMM'),
              ncol=3, with_dimplot=TRUE, dimplot_group='integrated_snn_res.0.8', point_size = 1)
feature_plots(Tcells_5, c('rna_CCR7', 'rna_SELL', 'rna_ZBTB16', 'rna_TBX21', 'rna_RORC', 'rna_GZMB', 'rna_PRF1', 'rna_IFNG',
                          'rna_TNF', 'rna_IL2', 'rna_IL17A'),
              ncol=3, with_dimplot=TRUE, dimplot_group='integrated_snn_res.0.8', point_size = 1)
# 2.3. Violin Plots to assess known markers:
violin_plots(Tcells_5, c('rna_CD3E', 'rna_CD3D', 'rna_CD3G', 'rna_CD3G', 'rna_KLRB1', 'rna_GNLY', 'rna_NKG7', 'rna_CD8A', 'rna_CD4',
                         'rna_TRDC', 'rna_TRGC1', 'rna_TRGC2', 'rna_TRAC', 'rna_TRBC1', 'rna_TRBC2', 'rna_CCR7', 'rna_CD62L'),
             ncol=3, with_dimplot=TRUE, group.by='integrated_snn_res.0.8')
# Resolution 0.8 will be used. Could most (except cluster 4) be MAIT cells?

# 3. Check markers between cluster 5 and all others, just to understand what this cells might be overall
# 3.1. Calculate all markers
Seurat::Idents(Tcells) = 'integrated_snn_res.0.1'
Tcells_5_markers_overall = Seurat::FindMarkers(Tcells, ident.1='5', assay='RNA', slot='data', logfc.threshold=.5, only.pos=TRUE)
View(Tcells_5_markers_overall)
# 3.2. Get top 20 markers for each sub-cluster
Tcells_5_markers_overall_top20 = c()
n_markers = dim(Tcells_5_markers_overall)[1]
order_by_avg_log2FC = order(Tcells_5_markers_overall$avg_log2FC, decreasing=T)
if(n_markers > 20) top20_clust_markers = Tcells_5_markers_overall[order_by_avg_log2FC[1:20],]
else top20_clust_markers = Tcells_5_markers_overall[order_by_avg_log2FC,]
Tcells_5_markers_overall_top20 = rbind(Tcells_5_markers_overall_top20, top20_clust_markers)
View(Tcells_5_markers_overall_top20)
Tcells_5_markers_overall[c(grep('^GZM', rownames(Tcells_5_markers_overall), value=T), 'KLRB1'),]

# 4. Find markers between the sub-clusters
# 4.1. Calculate all markers
Seurat::Idents(Tcells_5) = 'integrated_snn_res.0.8'
Tcells_5_markers = Seurat::FindAllMarkers(Tcells_5, assay='RNA', slot='data', logfc.threshold=.5, only.pos=TRUE)
View(Tcells_5_markers)
write.csv(Tcells_5_markers, paste(project_dir, '2_annotation/Results/Tcells/markers_C5_res08.csv', sep='/'))
# 4.2. Get top 20 markers for each sub-cluster
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
write.csv(Tcells_5_markers_top20, paste(project_dir, '2_annotation/Results/Tcells/markers_C5_res08_top20.csv', sep='/'))

# 7. Annotate sub-clusters:
# Sub-cluster 0           --> MAIT? NKT?
# Sub-cluster 1           --> MAIT? NKT?
# Sub-cluster 2           --> MAIT? NKT?
# Sub-cluster 3           --> MAIT? NKT?
# Sub-cluster 4           --> ???
# 7.1. Insert the annotations in the metadata:
Tcells_5[['Final_Annotation']] = rep('', length(Tcells_5$integrated_snn_res.0.8))
MAIT_NKT_cells = rownames(Tcells_5@meta.data)[Tcells_5$integrated_snn_res.0.8!='4']
Tcells_5@meta.data[MAIT_NKT_cells, 'Final_Annotation'] = 'MAIT? NKT?'
unkown_cells = rownames(Tcells_5@meta.data)[Tcells_5$integrated_snn_res.0.8=='4']
Tcells_5@meta.data[unkown_cells, 'Final_Annotation'] = 'Unkown'
# 7.2. Visualize annotations:
Seurat::DimPlot(Tcells_5, group.by='Final_Annotation', pt.size=0.5, label=T, label.size=4) +
  ggplot2::ggtitle('') + ggplot2::theme_minimal() + Seurat::NoLegend()
# 7.3. Save changes:
SeuratDisk::SaveH5Seurat(Tcells_5, paste(project_dir, '2_annotation/Results/Tcells/Tcells_5.h5Seurat', sep='/'), overwrite=TRUE)









# -----
# - Analyse clusters 1, 2, 3, 5, 6 together
# -----
Tcells = SeuratDisk::LoadH5Seurat(paste(project_dir, '2_annotation/data/Tcells_integratedClusters.h5Seurat', sep='/'))


# 1. Sub-cluster:
Tcells_wo04 = subset(Tcells, cells=rownames(Tcells@meta.data)[Tcells$integrated_snn_res.0.1%in%c('1', '2', '3', '5', '6')])
remove(Tcells)
invisible(gc())


# 2. Change metadata according to clusters performed individually:
Tcells_wo04[['Temporary_Annotation']] = rep('', dim(Tcells_wo04@meta.data)[1])
Tcells_wo04[['big_Tcell_clust']] = rep('', dim(Tcells_wo04@meta.data)[1])
Tcells_wo04[['orig_clusts']] = rep('', dim(Tcells_wo04@meta.data)[1])
# 2.1. C1
Tcells_1 = SeuratDisk::LoadH5Seurat(paste(project_dir, '2_annotation/Results/Tcells/Tcells_1.h5Seurat', sep='/'), overwrite=TRUE)
Tcells_wo04@meta.data[rownames(Tcells_1@meta.data), 'Temporary_Annotation'] = paste('C1', Tcells_1$Final_Annotation, sep='_')
Tcells_wo04@meta.data[rownames(Tcells_1@meta.data), 'big_Tcell_clust'] = rep('1', dim(Tcells_1@meta.data)[1])
Tcells_wo04@meta.data[rownames(Tcells_1@meta.data), 'orig_clusts'] = paste('C1', Tcells_1$integrated_snn_res.0.5, sep='_')
# 2.2. C2
Tcells_2 = SeuratDisk::LoadH5Seurat(paste(project_dir, '2_annotation/Results/Tcells/Tcells_2.h5Seurat', sep='/'), overwrite=TRUE)
Tcells_wo04@meta.data[rownames(Tcells_2@meta.data), 'Temporary_Annotation'] = paste('C2', Tcells_2$Final_Annotation, sep='_')
Tcells_wo04@meta.data[rownames(Tcells_2@meta.data), 'big_Tcell_clust'] = rep('2', dim(Tcells_2@meta.data)[1])
Tcells_wo04@meta.data[rownames(Tcells_2@meta.data), 'orig_clusts'] =  paste('C2', Tcells_2$integrated_snn_res.0.8, sep='_')
# 2.3. C3
Tcells_3 = SeuratDisk::LoadH5Seurat(paste(project_dir, '2_annotation/Results/Tcells/Tcells_3.h5Seurat', sep='/'), overwrite=TRUE)
Tcells_wo04@meta.data[rownames(Tcells_3@meta.data), 'Temporary_Annotation'] = paste('C3', Tcells_3$Final_Annotation, sep='_')
Tcells_wo04@meta.data[rownames(Tcells_3@meta.data), 'big_Tcell_clust'] = rep('3', dim(Tcells_3@meta.data)[1])
Tcells_wo04@meta.data[rownames(Tcells_3@meta.data), 'orig_clusts'] = paste('C3', Tcells_3$integrated_snn_res.0.8, sep='_')
# 2.4. C5
Tcells_5 = SeuratDisk::LoadH5Seurat(paste(project_dir, '2_annotation/Results/Tcells/Tcells_5.h5Seurat', sep='/'), overwrite=TRUE)
Tcells_wo04@meta.data[rownames(Tcells_5@meta.data), 'Temporary_Annotation'] = paste('C5', Tcells_5$Final_Annotation, sep='_')
Tcells_wo04@meta.data[rownames(Tcells_5@meta.data), 'big_Tcell_clust'] = rep('5', dim(Tcells_5@meta.data)[1])
Tcells_wo04@meta.data[rownames(Tcells_5@meta.data), 'orig_clusts'] = paste('C5', Tcells_5$integrated_snn_res.0.8, sep='_')
# 2.5. C6
Tcells_6 = SeuratDisk::LoadH5Seurat(paste(project_dir, '2_annotation/Results/Tcells/Tcells_6.h5Seurat', sep='/'), overwrite=TRUE)
Tcells_wo04@meta.data[rownames(Tcells_6@meta.data), 'Temporary_Annotation'] = paste('C6', Tcells_6$Final_Annotation, sep='_')
Tcells_wo04@meta.data[rownames(Tcells_6@meta.data), 'big_Tcell_clust'] = rep('6', dim(Tcells_6@meta.data)[1])
Tcells_wo04@meta.data[rownames(Tcells_6@meta.data), 'orig_clusts'] = paste('C6', rep('6', dim(Tcells_6@meta.data)[1]), sep='_')
# 2.6. Remove [integrated...] metadata variables
Tcells_wo04@meta.data = Tcells_wo04@meta.data[, !colnames(Tcells_wo04@meta.data) %in%
                                                c(grep('^in', colnames(Tcells_wo04@meta.data), value=T), 'seurat_clusters')]
invisible(gc())
# 2.7. Remove individual datasets:
remove(Tcells_1, Tcells_2, Tcells_3, Tcells_5, Tcells_6)
invisible(gc())


# 3. Re-group them (without clustering):
# 3.1. PCA
Tcells_wo04 = Seurat::RunPCA(Tcells_wo04, assay='integrated')
# 3.2. UMAP:
Tcells_wo04 = Seurat::RunUMAP(Tcells_wo04, dims=1:50, reduction='pca')
invisible(gc())
# 3.3. Save object
SeuratDisk::SaveH5Seurat(Tcells_wo04, paste(project_dir, '2_annotation/Results/Tcells/Tcells_wo04.h5Seurat', sep='/'), overwrite = TRUE)


# 4. Visualize the clusters:
big_Tcell_clust = Seurat::DimPlot(Tcells_wo04, group.by='big_Tcell_clust', label=T, label.size=4, pt.size=.5) +
  ggplot2::theme_minimal() + Seurat::NoLegend()
orig_clusts = Seurat::DimPlot(Tcells_wo04, group.by='orig_clusts', label=T, label.size=4, pt.size=.5) +
  ggplot2::theme_minimal() + Seurat::NoLegend()
Temporary_Annotation = Seurat::DimPlot(Tcells_wo04, group.by='Temporary_Annotation', label=T, label.size=4, pt.size=.5) +
  ggplot2::theme_minimal() + Seurat::NoLegend()
orig_annotation = Seurat::DimPlot(subset(Tcells_wo04, cells=rownames(Tcells_wo04@meta.data)[Tcells_wo04$original_annotation_2!='NA']),
                                  group.by='original_annotation_2', label=T, label.size=4, pt.size=.5) +
  ggplot2::theme_minimal() + Seurat::NoLegend()
(big_Tcell_clust | orig_clusts) / (Temporary_Annotation | orig_annotation)
orig_clusts | Temporary_Annotation

# 5. Assess known gene markers
# 5.1. Feature Plots:
feature_plots(Tcells_wo04, c('rna_CD3E', 'rna_CD3D', 'rna_CD3G'),
              ncol=2, with_dimplot=TRUE, dimplot_group='Temporary_Annotation', point_size = 1)
feature_plots(Tcells_wo04, c('rna_CD8A', 'rna_CD4', 'rna_TRDC', 'rna_TRGC1', 'rna_TRGC2'),
              ncol=3, with_dimplot=TRUE, dimplot_group='Temporary_Annotation', point_size = 1)
feature_plots(Tcells_wo04, c('rna_KLRB1', 'rna_GNLY', 'rna_NKG7', 'rna_GZMA', 'rna_GZMB', 'rna_GZMH', 'rna_GZMK', 'rna_GZMM'),
              ncol=3, with_dimplot=TRUE, dimplot_group='Temporary_Annotation', point_size = 1)
feature_plots(Tcells_wo04, c('rna_CCR7', 'rna_SELL', 'rna_PASK', 'rna_CREM', 'rna_CD69',
                          'rna_KLRB1', 'rna_S100A4', 'rna_PTPRCAP', 'rna_PRF1', 'rna_GZMA'),
              ncol=3, with_dimplot=TRUE, dimplot_group='Temporary_Annotation', point_size = 1)
feature_plots(Tcells_wo04, c('rna_TRAV1-2', 'rna_KLRB1', 'rna_ZBTB16', 'rna_TBX21', 'rna_RORC', 'rna_GZMB', 'rna_PRF1', 'rna_IFNG',
                          'rna_TNF', 'rna_IL2', 'rna_IL17A'),
              ncol=3, with_dimplot=TRUE, dimplot_group='Temporary_Annotation', point_size = 1)
# 2.3. Violin Plots:
violin_plots(Tcells_wo04, c('rna_CD3E', 'rna_CD3D', 'rna_CD3G', 'rna_CD8A', 'rna_CD4', 'rna_TRDC', 'rna_TRGC1', 'rna_TRGC2',
                         'rna_TRAC', 'rna_TRBC1', 'rna_TRBC2', 'rna_CCR7', 'rna_SELL', 'rna_PASK', 'rna_CREM', 'rna_CD69'),
             ncol=4, with_dimplot=TRUE, group.by='Temporary_Annotation')
violin_plots(Tcells_wo04, c('rna_KLRB1', 'rna_GNLY', 'rna_NKG7', 'rna_GZMA', 'rna_GZMB', 'rna_GZMH', 'rna_GZMK', 'rna_GZMM',
                            'rna_S100A4', 'rna_PTPRCAP', 'rna_PRF1', 'rna_GZMA'),
             ncol=4, with_dimplot=TRUE, group.by='Temporary_Annotation')

violin_plots(Tcells_wo04, c('rna_CD3E', 'rna_CD3D', 'rna_CD3G'),
             ncol=2, with_dimplot=TRUE, group.by='orig_clusts')
violin_plots(Tcells_wo04, c('rna_CD8A', 'rna_CD4', 'rna_TRDC', 'rna_TRGC1', 'rna_TRGC2',
                            'rna_TRAC', 'rna_TRBC1', 'rna_TRBC2'),
             ncol=3, with_dimplot=TRUE, group.by='orig_clusts')
violin_plots(Tcells_wo04, c('rna_CCR7', 'rna_SELL', 'rna_PASK', 'rna_CREM', 'rna_CD69'),
             ncol=3, with_dimplot=TRUE, group.by='orig_clusts')
violin_plots(Tcells_wo04, c('rna_KLRB1', 'rna_GNLY', 'rna_NKG7', 'rna_GZMA', 'rna_GZMB', 'rna_GZMH', 'rna_GZMK', 'rna_GZMM'),
             ncol=3, with_dimplot=TRUE, group.by='orig_clusts')

#'rna_S100A4', 'rna_PTPRCAP', 'rna_PRF1'


























