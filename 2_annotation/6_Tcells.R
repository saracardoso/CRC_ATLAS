project_dir = '/home/scardoso/Documents/PhD/CRC_ATLAS'
source(paste(project_dir, 'utils/modified_plots.R', sep='/'))



# -----
# - Integrate T cells
# -----

# Load subset of T cells:
Tcells = SeuratDisk::LoadH5Seurat('/home/scardoso/Documents/PhD/CRC_ATLAS/2_annotation/data/Tcells.h5Seurat')

# Separate T cells by dataset:
Tcells = Seurat::SplitObject(Tcells, split.by='dataset')
# SCTransform data:
for(dts_name in names(Tcells)) {
  Tcells[[dts_name]] <- Seurat::SCTransform(Tcells[[dts_name]], vars.to.regress = c('CC.Difference', 'percent.mitochondrial_RNA'))
  gc()
}
gc()
# Select integration features and prepare for integration:
features = Seurat::SelectIntegrationFeatures(Tcells, nfeatures = 3000)
Tcells = Seurat::PrepSCTIntegration(Tcells, anchor.features=features)
# Find integration achors and integrate the data:
anchors = Seurat::FindIntegrationAnchors(Tcells, normalization.method='SCT', anchor.features=features)
gc()
Tcells = Seurat::IntegrateData(anchorset=anchors, normalization.method="SCT")
gc()
# Save object:
SeuratDisk::SaveH5Seurat(Tcells, paste(project_dir, '2_annotation/data/Tcells_integratedClusters.h5Seurat', sep='/'))



# -----
# - Find First Clusters and Annotate Them
# -----

# 1. Run PCA:
Tcells = Seurat::RunPCA(Tcells, assay='integrated')

# 1. Find first clusters
# 2.1. Determine the K-nearest neighbor graph
Tcells = Seurat::FindNeighbors(Tcells, dims=1:50)
# 2.2. Find clusters:
Tcells = Seurat::FindClusters(Tcells, resolution=seq(0.1, 1, by=.1))
# 2.3. UMAP:
Tcells = Seurat::RunUMAP(Tcells, dims=1:50, reduction='pca')
invisible(gc())

# 3. Choose best resolution using clustree:
library(ggraph)
clust_tree = clustree::clustree(Tcells, prefix='integrated_snn_res.')
clust_tree

# 4. Visualize different resolutions (first resolution was chosen to be the best for now):
clusters_01 = Seurat::DimPlot(Tcells, reduction="umap", group.by='integrated_snn_res.0.1', label=TRUE, label.size=6, pt.size=.2)
clusters_02 = Seurat::DimPlot(Tcells, reduction="umap", group.by='integrated_snn_res.0.2', label=TRUE, label.size=6, pt.size=.2)
clusters_03 = Seurat::DimPlot(Tcells, reduction="umap", group.by='integrated_snn_res.0.3', label=TRUE, label.size=6, pt.size=.2)
clusters_05 = Seurat::DimPlot(Tcells, reduction="umap", group.by='integrated_snn_res.0.5', label=TRUE, label.size=6, pt.size=.2)
((clusters_01 | clusters_02) / (clusters_03 | clusters_05)) | clust_tree

# 5. Assess Known gene markers:
# 5.1. Feature Plots:
feature_plots(Tcells, c('rna_CD3E', 'rna_CD3D', 'rna_CD8A', 'rna_CD4', 'rna_TRDC', 'rna_TRGC1', 'rna_TRGC2'),
              ncol=2, with_dimplot=TRUE, dimplot_group='integrated_snn_res.0.1')
feature_plots(Tcells, c('rna_TRAC', 'rna_TRBC1', 'rna_TRBC2', 'rna_PDCD1', 'rna_GZMB', 'rna_PRF1', 'rna_CCR7'),
              ncol=2, with_dimplot=TRUE, dimplot_group='integrated_snn_res.0.1')
feature_plots(Tcells, c('rna_CXCL13', 'rna_IL2RB', 'rna_KLRB1', 'rna_CD69'),
              ncol=2, with_dimplot=TRUE, dimplot_group='integrated_snn_res.0.1')
#Seurat::FeaturePlot(Tcells_subset, c('rna_TRAC', 'rna_TRDC'), blend = TRUE) / clusters_01
# 5.2. Violin Plots:
violin_plots(Tcells, c('rna_CD3E', 'rna_CD3D', 'rna_CD8A', 'rna_CD4', 'rna_TRDC', 'rna_TRGC1',
                       'rna_TRGC2', 'rna_TRAC', 'rna_TRBC1', 'rna_TRBC2', 'rna_PDCD1'),
              ncol=3, with_dimplot=TRUE, group.by='integrated_snn_res.0.1')
violin_plots(Tcells, c('rna_GZMB', 'rna_PRF1', 'rna_CCR7', 'rna_CXCL13', 'rna_IL2RB', 'rna_KLRB1', 'rna_CD69'),
              ncol=3, with_dimplot=TRUE, group.by='integrated_snn_res.0.1')

# 6. Save object:
SeuratDisk::SaveH5Seurat(Tcells, paste(project_dir, '2_annotation/data/Tcells_integratedClusters.h5Seurat', sep='/'))



# ---
# - Check SIGMA for further clusterability (res.0.1)
# ---

# 1. Get data processed for SIGMA
# 1.1. Subset cells to half of the original
n_cells = round(dim(Tcells@assays$integrated@data)[2] / 3)
sampled_cells = colnames(Tcells@assays$integrated@data)[sample(1:dim(Tcells@assays$integrated@data)[2], n_cells)]
# 1.2. Check if subsetting still resembles unsubsetted dataset
clusters_01 | Seurat::DimPlot(subset(Tcells, cells=sampled_cells), reduction="umap", group.by='integrated_snn_res.0.1',
                              label=TRUE, label.size=6, pt.size=.2)
# 1.3. Subset dataset
Tcells_subset = subset(Tcells, cells=sampled_cells)
invisible(gc())
remove(Tcells)
invisible(gc())
# 1.4. Get integrated dataset (because of the integration between datasets)
sigma_data = as.matrix(SeuratObject::GetAssayData(Tcells_subset, assay='integrated', slot='data'))
invisible(gc())
# 1.5. Get vector with resolution 0.1 clusters and dataset
clusters_01_vector = Tcells_subset@meta.data[colnames(sigma_data), 'integrated_snn_res.0.1']
dataset_vector = Tcells_subset@meta.data[colnames(sigma_data), 'dataset']

# 2. Run SIGMA
SIGMA_all_res01 = SIGMA::sigma_funct(sigma_data, clusters_01_vector)
# 2.1. Save SIGMA result
saveRDS(SIGMA_all_res01, paste(project_dir, '2_annotation/Results/Tcells/SIGMA_all_res01.Rdata', sep='/'))

# 3. Check results
# 3.1. Clusterability of all clusters
SIGMA::plot_sigma(SIGMA_all_res01) # All clusters have a sigma value greater than 0.90.
SIGMA_all_res01$maximum_measure # min (Cluster 4) = 0.9424960; max (Cluster 0) = 0.9761377
# 3.2. Check that clusterability is not due to datasets
SIGMA::plot_singular_vectors(SIGMA_all_res01, '0', colour=dataset_vector[Tcells_subset@meta.data$integrated_snn_res.0.1=='0']) +
  ggplot2::ggtitle('Cluster 0')
SIGMA::plot_singular_vectors(SIGMA_all_res01, '1', colour=dataset_vector[Tcells_subset@meta.data$integrated_snn_res.0.1=='1']) +
  ggplot2::ggtitle('Cluster 1')
SIGMA::plot_singular_vectors(SIGMA_all_res01, '2', colour=dataset_vector[Tcells_subset@meta.data$integrated_snn_res.0.1=='2']) +
  ggplot2::ggtitle('Cluster 2')
SIGMA::plot_singular_vectors(SIGMA_all_res01, '3', colour=dataset_vector[Tcells_subset@meta.data$integrated_snn_res.0.1=='3']) +
  ggplot2::ggtitle('Cluster 3')
SIGMA::plot_singular_vectors(SIGMA_all_res01, '4', colour=dataset_vector[Tcells_subset@meta.data$integrated_snn_res.0.1=='4']) +
  ggplot2::ggtitle('Cluster 4')
SIGMA::plot_singular_vectors(SIGMA_all_res01, '5', colour=dataset_vector[Tcells_subset@meta.data$integrated_snn_res.0.1=='5']) +
  ggplot2::ggtitle('Cluster 5')
SIGMA::plot_singular_vectors(SIGMA_all_res01, '6', colour=dataset_vector[Tcells_subset@meta.data$integrated_snn_res.0.1=='6']) +
  ggplot2::ggtitle('Cluster 6')



# -----
# - Sub-cluster cluster 0 from resolution 0.1 separately
# -----
Tcells = SeuratDisk::LoadH5Seurat(paste(project_dir, '2_annotation/data/Tcells_integratedClusters.h5Seurat', sep='/'))


# 1. Sub-cluster:
Tcells_0 = subset(Tcells, cells=rownames(Tcells@meta.data)[Tcells$integrated_snn_res.0.1=='0'])
Seurat::DimPlot(Tcells, group.by = 'integrated_snn_res.0.1') | Seurat::DimPlot(Tcells_0, group.by = 'integrated_snn_res.0.1')
# 1.1. PCA
Tcells_0 = Seurat::RunPCA(Tcells_0, assay='integrated')
# 1.2. Determine the K-nearest neighbor graph
Tcells_0 = Seurat::FindNeighbors(Tcells_0, dims=1:50)
# 1.3. Find clusters:
Tcells_0 = Seurat::FindClusters(Tcells_0, resolution=seq(0.1, 1, by=.1))
# 1.4. UMAP:
Tcells_0 = Seurat::RunUMAP(Tcells_0, dims=1:50, reduction='pca')
invisible(gc())
# 1.5. Save object
SeuratDisk::SaveH5Seurat(Tcells_0, paste(project_dir, '2_annotation/Results/Tcells/Tcells_0.h5Seurat', sep='/'))


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
feature_plots(Tcells_0, c('rna_CD3E', 'rna_CD3D', 'rna_CD8A', 'rna_CD4', 'rna_TRAC', 'rna_TRBC1', 'rna_TRBC2'),
              ncol=2, with_dimplot=TRUE, dimplot_group='integrated_snn_res.0.3')
feature_plots(Tcells_0, c('rna_CCR7', 'rna_PDCD1', 'rna_IL2RB', 'rna_KLRB1', 'rna_CD69'),
              ncol=2, with_dimplot=TRUE, dimplot_group='integrated_snn_res.0.3')
feature_plots(Tcells_0, c('rna_CCR7', 'rna_SELL', 'rna_IL2RB', 'rna_IL2RA', 'rna_CD69', 'rna_CREM'),
              ncol=2, with_dimplot=TRUE, dimplot_group='integrated_snn_res.0.4')
feature_plots(Tcells_0, c('rna_CCR7', 'rna_SELL', 'rna_NPM1', 'rna_LRRN3', 'rna_PASK', 'rna_IL7R'),
              ncol=2, with_dimplot=TRUE, dimplot_group='integrated_snn_res.0.4')
# 2.3. Violin Plots to assess known markers:
violin_plots(Tcells_0, c('rna_CD3E', 'rna_CD3D', 'rna_CD8A', 'rna_CD4', 'rna_TRAC', 'rna_TRBC1', 'rna_TRBC2'),
             ncol=3, with_dimplot=TRUE, group.by='integrated_snn_res.0.3')
violin_plots(Tcells_0, c('rna_CCR7', 'rna_PDCD1', 'rna_IL2RB', 'rna_KLRB1', 'rna_CD69'),
             ncol=3, with_dimplot=TRUE, group.by='integrated_snn_res.0.3')
violin_plots(Tcells_0, c('rna_CCR7', 'rna_SELL', 'rna_IL2RB', 'rna_IL2RA', 'rna_CD69', 'rna_CREM'),
             ncol=3, with_dimplot=TRUE, group.by='integrated_snn_res.0.4')
# 2.4. Resolution 0.4 will be chosen


# 3. Find markers
# 3.1. Calculate all markers
Seurat::Idents(Tcells_0) = 'integrated_snn_res.0.4'
Tcells_0_markers = Seurat::FindAllMarkers(Tcells_0, assay='RNA', slot='data', logfc.threshold=.5, only.pos=TRUE)
View(Tcells_0_markers)
write.csv(Tcells_0_markers, paste(project_dir, '2_annotation/Results/Tcells/markers_C0_res04.csv', sep='/'))
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
write.csv(Tcells_0_markers_top20, paste(project_dir, '2_annotation/Results/Tcells/markers_C0_res04_top20.csv', sep='/'))


# 4. Annotate sub-clusters:
# Sub-cluster 0        --> CD4 central memory (CM)
# Sub-cluster 1        --> naive CD4
# Sub-cluster 2        --> CD4, TEMRA?
# Sub-clusters 3, 5, 6 --> CD4, TEM? TEMRA?
# Sub-cluster 7        --> Th17
# Sub-cluster8         --> IFN response related CD4


# 5. Perform trajectory analysis to help in this annotation:
# 5.1. Convert Seurat dataset to monocle and use monocle's UMAP:
monocle_Tcells_0 = SeuratWrappers::as.cell_data_set(Tcells_0, assay='integrated')
monocle_Tcells_0 = monocle3::reduce_dimension(monocle_Tcells_0, reduction_method='UMAP', preprocess_method='PCA')
seurat_clusters = monocle3::plot_cells(monocle_Tcells_0, show_trajectory_graph=FALSE, color_cells_by="integrated_snn_res.0.4",
                                       group_label_size=5, cell_size=.8) + ggplot2::ggtitle("Colored by Seurat's clusters")
# 5.2. Cluster cells to get different partitions (cells belonging to different trajectories, i.e., different ancestors):
monocle_Tcells_0 = monocle3::cluster_cells(monocle_Tcells_0)
p2 = monocle3::plot_cells(monocle_Tcells_0, color_cells_by="cluster", show_trajectory_graph=FALSE,
                          cell_size=.8, group_label_size=5) + ggplot2::ggtitle("Colored by Monocle's clusters")
p3 = monocle3::plot_cells(monocle_Tcells_0, color_cells_by="partition", show_trajectory_graph=FALSE,
                          cell_size=.8, group_label_size=5) + ggplot2::ggtitle("Monocle's partitions")
seurat_clusters | p2 | p3
# 5.3. Learn the trajectory graph:
monocle_Tcells_0 = monocle3::learn_graph(monocle_Tcells_0)
tj1 = monocle3::plot_cells(monocle_Tcells_0, color_cells_by = "integrated_snn_res.0.4", label_cell_groups=FALSE, labels_per_group=2,
                           label_leaves=TRUE, label_branch_points=TRUE, graph_label_size=2, group_label_size=6, cell_size=.8) +
  ggplot2::ggtitle("Trajectory") + Seurat::NoLegend()
(seurat_clusters / p2) | tj1
# 5.4. Order cells in pseudo-time
monocle_Tcells_0 = monocle3::order_cells(monocle_Tcells_0)
tjtime =  monocle3::plot_cells(monocle_Tcells_0, color_cells_by = "pseudotime", label_cell_groups=FALSE,
                               label_leaves=FALSE, label_branch_points=FALSE, graph_label_size=3, group_label_size=6, cell_size=.8) +
  ggplot2::ggtitle("Trajectory coloured by pseudo-time") + Seurat::NoLegend()
tj1 | tjtime
# 5.5. Visualize cells ordered by pseudotime
pdata_monocle_Tcells_0 = monocle3::pData(monocle_Tcells_0)
pdata_monocle_Tcells_0$pseudotime_monocle3 = monocle3::pseudotime(monocle_Tcells_0)
pdata_monocle_Tcells_0$integrated_snn_res.0.4 = factor(pdata_monocle_Tcells_0$integrated_snn_res.0.4,
                                                       levels=c('4', '0', '6', '3', '5', '2', '1', '7', '8'))
ggplot2::ggplot(as.data.frame(pdata_monocle_Tcells_0), 
                ggplot2::aes(x = pseudotime_monocle3, 
                             y = integrated_snn_res.0.4, colour = integrated_snn_res.0.4)) +
  ggbeeswarm::geom_quasirandom(groupOnX = FALSE) + ggplot2::theme_classic() +
  ggplot2::xlab("monocle3 pseudotime") + ggplot2::ylab("Seurat Clusters") +
  ggplot2::ggtitle("Cells ordered by monocle3 pseudotime") + Seurat::NoLegend()
# 5.6. Visualize ordered pseudotime without Treg branch
monocle_Tcells_0_noTreg = monocle3::choose_cells(monocle_Tcells_0)
pdata_monocle_Tcells_0_noTreg = monocle3::pData(monocle_Tcells_0_noTreg)
pdata_monocle_Tcells_0_noTreg$pseudotime_monocle3 = monocle3::pseudotime(monocle_Tcells_0_noTreg)
pdata_monocle_Tcells_0_noTreg$integrated_snn_res.0.4 = factor(pdata_monocle_Tcells_0_noTreg$integrated_snn_res.0.4,
                                                       levels=c('4', '0', '6', '3', '5', '2', '1', '7', '8'))
ggplot2::ggplot(as.data.frame(pdata_monocle_Tcells_0_noTreg), 
                ggplot2::aes(x = pseudotime_monocle3, 
                             y = integrated_snn_res.0.4, colour = integrated_snn_res.0.4)) +
  ggbeeswarm::geom_quasirandom(groupOnX = FALSE) + ggplot2::theme_classic() +
  ggplot2::xlab("monocle3 pseudotime") + ggplot2::ylab("Seurat Clusters") +
  ggplot2::ggtitle("Cells ordered by monocle3 pseudotime without Treg branch") + Seurat::NoLegend()


# 6. Annotate sub-clusters:
# Sub-cluster 0        --> CD4 central memory (CM)
# Sub-cluster 1        --> naive CD4
# Sub-cluster 2, 3, 5  --> CD4, TEM (TEMRA?)
# Sub-clusters 6       --> CD4, TEM?
# Sub-cluster 7        --> Th17
# Sub-cluster8         --> IFN response related CD4


# 7. Check previous SIGMA result and colour cluster 0 by new sub-clusters
SIGMA_all_res01 = readRDS(paste(project_dir, '2_annotation/Results/Tcells/SIGMA_all_res01.Rdata', sep='/'))
cluster_0_cells = rownames(Tcells_0@meta.data)
subclusters_0 = Tcells_0$integrated_snn_res.0.4[intersect(colnames(SIGMA_all_res01$input_parameters$expr), cluster_0_cells)]
subclusters_firstAnnotations = as.character(subclusters_0)
subclusters_firstAnnotations[subclusters_firstAnnotations%in%c('2','3','5')] = 'TEMRA?'
SIGMA::plot_singular_vectors(SIGMA_all_res01, '0', colour=subclusters_0)
SIGMA::plot_singular_vectors(SIGMA_all_res01, '0', colour=subclusters_firstAnnotations)


# 7. Check SIGMA for further clusterability
# 7.1. Subset cells to half of the original
n_cells = round(dim(Tcells_0@assays$integrated@data)[2] / 3)
sampled_cells = colnames(Tcells_0@assays$integrated@data)[sample(1:dim(Tcells_0@assays$integrated@data)[2], n_cells)]
# 7.2. Check if subsetting still resembles unsubsetted dataset
clusters_04_0 | Seurat::DimPlot(subset(Tcells_0, cells=sampled_cells), reduction="umap", group.by='integrated_snn_res.0.4',
                                label=TRUE, label.size=6, pt.size=.2)
# 7.3. Subset dataset
Tcells_0_subset = subset(Tcells_0, cells=sampled_cells)
invisible(gc())
remove(Tcells_0)
invisible(gc())
# 7.4. Get integrated dataset (because of the integration between datasets)
sigma_data = as.matrix(SeuratObject::GetAssayData(Tcells_0_subset, assay='integrated', slot='data'))
invisible(gc())
# 7.5. Get vector with resolution 0.4 clusters
clusters_04_vector = Tcells_0_subset@meta.data[colnames(sigma_data), 'integrated_snn_res.0.4']
dataset_vector = Tcells_0_subset@meta.data[colnames(sigma_data), 'dataset']
# 7.6. Run SIGMA
SIGMA_C0_res04 = SIGMA::sigma_funct(sigma_data, clusters_04_vector)
invisible(gc())
# 7.7. Save SIGMA result
saveRDS(SIGMA_C0_res04, paste(project_dir, '2_annotation/Results/Tcells/SIGMA_C0_res04.Rdata', sep='/'))
# 7.8. Check results
# 7.8.1. Clusterability of all clusters
SIGMA::plot_sigma(SIGMA_C0_res04) # All clusters have a sigma value greater than 0.8
SIGMA_C0_res04$maximum_measure # min (Cluster 4) = 0.8873414; max (Cluster 0) = 0.9594508
# 7.8.2. Check that clusterability is not due to datasets
SIGMA::plot_singular_vectors(SIGMA_C0_res04, '0', colour=dataset_vector[Tcells_0_subset@meta.data$integrated_snn_res.0.4=='0']) +
  ggplot2::ggtitle('Cluster 0')
SIGMA::plot_singular_vectors(SIGMA_C0_res04, '1', colour=dataset_vector[Tcells_0_subset@meta.data$integrated_snn_res.0.4=='1']) +
  ggplot2::ggtitle('Cluster 1')
SIGMA::plot_singular_vectors(SIGMA_C0_res04, '2', colour=dataset_vector[Tcells_0_subset@meta.data$integrated_snn_res.0.4=='2']) +
  ggplot2::ggtitle('Cluster 2')
SIGMA::plot_singular_vectors(SIGMA_C0_res04, '3', colour=dataset_vector[Tcells_0_subset@meta.data$integrated_snn_res.0.4=='3']) +
  ggplot2::ggtitle('Cluster 3')
SIGMA::plot_singular_vectors(SIGMA_C0_res04, '4', colour=dataset_vector[Tcells_0_subset@meta.data$integrated_snn_res.0.4=='4']) +
  ggplot2::ggtitle('Cluster 4')
SIGMA::plot_singular_vectors(SIGMA_C0_res04, '5', colour=dataset_vector[Tcells_0_subset@meta.data$integrated_snn_res.0.4=='5']) +
  ggplot2::ggtitle('Cluster 5')
SIGMA::plot_singular_vectors(SIGMA_C0_res04, '6', colour=dataset_vector[Tcells_0_subset@meta.data$integrated_snn_res.0.4=='6']) +
  ggplot2::ggtitle('Cluster 6')
SIGMA::plot_singular_vectors(SIGMA_C0_res04, '7', colour=dataset_vector[Tcells_0_subset@meta.data$integrated_snn_res.0.4=='7']) +
  ggplot2::ggtitle('Cluster 7')
SIGMA::plot_singular_vectors(SIGMA_C0_res04, '8', colour=dataset_vector[Tcells_0_subset@meta.data$integrated_snn_res.0.4=='8']) +
  ggplot2::ggtitle('Cluster 8')


# 8. For now, these cluster will not be further clustered
























# -----
# - Visualize Clusters
# -----

# Cluster from the three resolutions:
clusters_04 = Seurat::DimPlot(Tcells, reduction="umap", group.by='integrated_snn_res.0.4', label=TRUE, label.size=6)
clusters_06 = Seurat::DimPlot(Tcells, reduction="umap", group.by='integrated_snn_res.0.6', label=TRUE, label.size=6)
clusters_1 = Seurat::DimPlot(Tcells, reduction="umap", group.by='integrated_snn_res.1', label=TRUE, label.size=6)
(clusters_04 | clusters_06) / (clusters_1 | ggplot2::ggplot())

# Choose resolution 0.6
res = 'integrated_snn_res.1'

# Metadata:
state = Seurat::DimPlot(Tcells, reduction="umap", group.by='state', label=FALSE, label.size=6)
patients = Seurat::DimPlot(Tcells, reduction="umap", group.by='patient', label=FALSE, label.size=6, pt.size=.1)
datasets = Seurat::DimPlot(Tcells, reduction="umap", group.by='dataset', label=FALSE, label.size=6, pt.size=.1)
(clusters_1 | patients) / (state | datasets)

# Normal/ tumour/ healthy distribution across clusters:
clusters = c()
percs = c()
fill = c()
for(cluster in unique(as.character(Tcells@meta.data[,res]))){
  x = table(Tcells$state[Tcells@meta.data[,res]==cluster])[c('Tumor', 'Normal', 'Healthy')]
  x[is.na(x)] = 0
  y = sum(x)
  percs = c(percs, x/y)
  fill = c(fill, c('Tumor', 'Normal', 'Healthy'))
  clusters = c(clusters, rep(cluster, 3))
}
dist_nt = data.frame(percs=percs, fill=fill, clusters=clusters)
dist_nt$clusters = factor(dist_nt$clusters,
                          levels=as.character(0:(length(unique(Tcells@meta.data[,res]))-1)))
dist_state = ggplot2::ggplot(dist_nt, mapping = ggplot2::aes(clusters, percs, fill=fill)) +
  ggplot2::geom_bar(position='stack', stat='identity')
dist_state

# Gene expression:
Seurat::VlnPlot(Tcells, features=c('rna_CD3E', 'rna_CD3D', 'rna_CD8A', 'rna_CD4', 'rna_TRDC', 'rna_TRGC1'),
                pt.size = 0.001, ncol=2, group.by = res)# + ggmin::theme_min()
Seurat::VlnPlot(Tcells, features=c('rna_TRGC2', 'rna_TRAC', 'rna_TRBC1', 'rna_TRBC2', 'rna_PDCD1', 'rna_GZMB'),
                pt.size = 0.001, ncol=2, group.by = res)# + ggmin::theme_min()
Seurat::VlnPlot(Tcells, features=c('rna_PRF1', 'rna_CXCL13', 'rna_CCR7', 'rna_IL2RB', 'rna_KLRB1', 'rna_CD69'),
                pt.size = 0.001, ncol=2, group.by = res)# + ggmin::theme_min()

Seurat::FeaturePlot(Tcells, features=c('rna_CD3E', 'rna_CD3D', 'rna_CD8A', 'rna_CD4', 'rna_TRDC', 'rna_TRGC1'))
Seurat::FeaturePlot(Tcells, features=c('rna_TRGC2', 'rna_TRAC', 'rna_TRBC1', 'rna_TRBC2', 'rna_PDCD1', 'rna_GZMB'))
Seurat::FeaturePlot(Tcells, features=c('rna_PRF1', 'rna_CXCL13', 'rna_CCR7', 'rna_IL2RB', 'rna_KLRB1', 'rna_CD69'))

Seurat::DotPlot(Tcells, features=c('rna_CD3E', 'rna_CD3D', 'rna_CD8A', 'rna_CD4', 'rna_TRDC', 'rna_TRGC1',
                                   'rna_TRGC2', 'rna_TRAC', 'rna_TRBC1', 'rna_TRBC2', 'rna_PDCD1', 'rna_GZMB',
                                   'rna_PRF1', 'rna_CXCL13', 'rna_CCR7', 'rna_IL2RB', 'rna_KLRB1', 'rna_CD69'), assay='RNA',
                group.by=res) + ggmin::theme_min()





# -----
# - Find markers
# -----

Seurat::Idents(Tcells) = res
Tcells_markers = Seurat::FindAllMarkers(Tcells, assay='RNA', slot='data', logfc.threshold=.8, only.pos=TRUE)
View(Tcells_markers)
write.csv(Tcells_markers, paste(project_dir, '2_annotation/markers/Tcells/Tcells_res1.csv', sep='/'))

# Get top 20 markers for each cluster:
Tcells_markers_top20 = dplyr::top_n(dplyr::group_by(Tcells_markers[Tcells_markers$p_val_adj<0.05,], cluster), n=20, wt=avg_log2FC)
View(Tcells_markers_top20)
write.csv(Tcells_markers_top20, paste(project_dir, '2_annotation/markers/Tcells/Tcells_res1_top20.csv', sep='/'))




# -----
# - Store Tcells data for cellxgene
# -----


for(column in colnames(Tcells@meta.data)){
  if(class(Tcells@meta.data[[column]])=='factor') print(column)
}
loomR::create('/home/scardoso/Documents/PhD/CRC_ATLAS/cellxgene/TcellsLoom.loom', data=Seurat::GetAssayData(Tcells, assay='RNA', slot='data'),
              cell.attrs=as.list(Tcells@meta.data))

# Get necessary python modules:
reticulate::py_config()
anndata_module = reticulate::import('anndata')

# Initialize the h5ad object:
scData = anndata_module$AnnData(Seurat::GetAssayData(Tcells, assay='RNA', slot='data'))
invisible(gc())
scData = scData$transpose()
scData$obs = Tcells@meta.data
scData$write_h5ad('/home/scardoso/Documents/PhD/CRC_ATLAS/cellxgene/Tcells.h5ad')
