project_dir = '/home/scardoso/Documents/PhD/CRC_ATLAS'
source(paste(project_dir, 'utils/modified_plots.R', sep='/'))

# Load subset of myeloid cells:
Myeloid = SeuratDisk::LoadH5Seurat(paste(project_dir, '/2_annotation/results_Myeloid/datasets/myeloid_cells.h5Seurat', sep='/'))





# -----
# - Find clusters
# -----


# 1. Run PCA:
Myeloid = Seurat::RunPCA(Myeloid, assay='integrated')
invisible(gc())
Seurat::DefaultAssay(Myeloid)='integrated'

# 2. Choose number of PCs to use (where the elbow occurs):
pct = Myeloid[["pca"]]@stdev /sum(Myeloid[["pca"]]@stdev) * 100
cumu = cumsum(pct)
co1 = which(cumu > 90 & pct < 5)[1]
co2 = sort(which((pct[1:length(pct) - 1] - pct[2:length(pct)]) > 0.1), decreasing = T)[1] + 1
# 2.1 Number of PCs to use
elbow = min(co1, co2) # 12
# 2.2. Visual representation of the PCs to use
plot_df = data.frame(pct=pct, cumu=cumu, rank=1:length(pct))
ggplot2::ggplot(plot_df, ggplot2::aes(cumu, pct, label = rank, color = rank > elbow)) + 
  ggplot2::geom_text() + 
  ggplot2::geom_vline(xintercept = 90, color = "grey") + 
  ggplot2::geom_hline(yintercept = min(pct[pct > 5]), color = "grey") +
  ggplot2::theme_bw()
# 2.3. Heatmap of the first 12 PCs
Seurat::DimHeatmap(Myeloid, dims=1:12, cells=500, balanced=TRUE)


# 3. Find clusters
# 3.1. Determine the K-nearest neighbor graph
Myeloid = Seurat::FindNeighbors(Myeloid, dims=1:elbow)
# 3.2. Find clusters:
Myeloid = Seurat::FindClusters(Myeloid, resolution=seq(0.1, 1, by=.1))
# 3.3. UMAP:
Myeloid = Seurat::RunUMAP(Myeloid, dims=1:elbow, reduction='pca')
invisible(gc())


# 4. Choose best resolution:
library(ggraph)
clust_tree = clustree::clustree(Myeloid, prefix='integrated_snn_res.')
clust_tree
# 4.1. Visualize different resolutions (first resolution was chosen to be the best for now):
clusters_01 = Seurat::DimPlot(Myeloid, reduction="umap", group.by='integrated_snn_res.0.1', label=TRUE, label.size=6, pt.size=.2)
clusters_04 = Seurat::DimPlot(Myeloid, reduction="umap", group.by='integrated_snn_res.0.4', label=TRUE, label.size=6, pt.size=.2)
clusters_07 = Seurat::DimPlot(Myeloid, reduction="umap", group.by='integrated_snn_res.0.7', label=TRUE, label.size=6, pt.size=.2)
clusters_08 = Seurat::DimPlot(Myeloid, reduction="umap", group.by='integrated_snn_res.0.8', label=TRUE, label.size=6, pt.size=.2)
((clusters_01 | clusters_04) / (clusters_07 | clusters_08)) | clust_tree
# 4.2. Plot Metrics
metrics =  c("nUMI", "nGene", "S.Score", "G2M.Score", "percent.mitochondrial_RNA")
feature_plots(Myeloid, metrics, ncol=3, with_dimplot=TRUE, dimplot_group='integrated_snn_res.0.7')
# 4.3. Plot metadata
state = Seurat::DimPlot(Myeloid, reduction="umap", group.by='state', label=FALSE, label.size=6) + Seurat::NoLegend()
patients = Seurat::DimPlot(Myeloid, reduction="umap", group.by='patient', label=FALSE, label.size=6, pt.size=.1) + Seurat::NoLegend()
datasets = Seurat::DimPlot(Myeloid, reduction="umap", group.by='dataset', label=FALSE, label.size=6, pt.size=.1) + Seurat::NoLegend()
(clusters_08 | patients) / (state | datasets)
# 4.4. Normal/ tumour/ healthy distribution across clusters:
clusters = c()
percs = c()
fill = c()
for(cluster in unique(as.character(Myeloid@meta.data[,'integrated_snn_res.0.8']))){
  x = table(Myeloid$state[Myeloid@meta.data[,'integrated_snn_res.0.8']==cluster])[c('Tumor', 'Normal', 'Healthy')]
  x[is.na(x)] = 0
  y = sum(x)
  percs = c(percs, x/y)
  fill = c(fill, c('Tumor', 'Normal', 'Healthy'))
  clusters = c(clusters, rep(cluster, 3))
}
dist_nt = data.frame(percs=percs, fill=fill, clusters=clusters)
dist_nt$clusters = factor(dist_nt$clusters,
                          levels=as.character(0:(length(unique(Myeloid@meta.data[,'integrated_snn_res.0.8']))-1)))
dist_state = ggplot2::ggplot(dist_nt, mapping = ggplot2::aes(clusters, percs, fill=fill)) +
  ggplot2::geom_bar(position='stack', stat='identity')
dist_state
# 4.5. Feature Plots:
feature_plots(Myeloid, c('rna_KIT', 'rna_TPSB2', 'rna_TPSAB1'),
              ncol=2, with_dimplot=TRUE, dimplot_group='integrated_snn_res.0.8') # Mast cells
feature_plots(Myeloid, c('rna_PPBP', 'rna_PF4'),
              ncol=2, with_dimplot=TRUE, dimplot_group='integrated_snn_res.0.8') # Megakaryocytes / Platelets
feature_plots(Myeloid, c('rna_FCER1A', 'rna_CST3'),
              ncol=2, with_dimplot=TRUE, dimplot_group='integrated_snn_res.0.8') # cDCs
feature_plots(Myeloid, c('rna_IL3RA', 'rna_SERPINF1', 'rna_GZMB', 'rna_ITM2C'),
              ncol=3, with_dimplot=TRUE, dimplot_group='integrated_snn_res.0.8') # pDCs
feature_plots(Myeloid, c('rna_CD14', 'rna_FCGR3A', 'rna_MARCO', 'rna_ADGRE1', 'rna_ITGAM'),
              ncol=2, with_dimplot=TRUE, dimplot_group='integrated_snn_res.0.8') # Monocyte/Macrophage lineage
# 4.6. Resolution 0.8 is chosen
# 4.7. Known clusters for now:
# Cluster 0           --> cDCs
# Clusters 1-4, 6-12  --> Monocyte/Macrophage lineage
# Cluster 5           --> Mast cells
# Cluster 14          --> pDCs
# Cluster 13          --> Unknown

# 5. Save object:
SeuratDisk::SaveH5Seurat(Myeloid, paste(project_dir, '2_annotation/results_Myeloid/datasets/Myeloid_clustered.h5Seurat', sep='/'))




# -----
# - Cluster only monocyte/macropphage clusters
# -----

# 1. Sub-cluster:
monocytes_macrophages = subset(Myeloid, cells=rownames(Myeloid@meta.data)[!Myeloid$integrated_snn_res.0.8%in%c('0', '5', '13', '14')])
(Seurat::DimPlot(Myeloid, group.by = 'integrated_snn_res.0.8', label=T, label.size=4) + ggplot2::theme_minimal() + Seurat::NoLegend()) |
  (Seurat::DimPlot(monocytes_macrophages, group.by = 'integrated_snn_res.0.8', label=T, label.size=4) + ggplot2::theme_minimal() + Seurat::NoLegend())
# 1.1. PCA
monocytes_macrophages = Seurat::RunPCA(monocytes_macrophages, assay='integrated')
# 1.2. Choose number of PCs to use
pct = monocytes_macrophages[["pca"]]@stdev /sum(monocytes_macrophages[["pca"]]@stdev) * 100
cumu = cumsum(pct)
co1 = which(cumu > 90 & pct < 5)[1]
co2 = sort(which((pct[1:length(pct) - 1] - pct[2:length(pct)]) > 0.1), decreasing = T)[1] + 1
elbow = min(co1, co2) # 10
# 1.3. Visual representation of the PCs to use
plot_df = data.frame(pct=pct, cumu=cumu, rank=1:length(pct))
ggplot2::ggplot(plot_df, ggplot2::aes(cumu, pct, label = rank, color = rank > elbow)) + 
  ggplot2::geom_text() + 
  ggplot2::geom_vline(xintercept = 90, color = "grey") + 
  ggplot2::geom_hline(yintercept = min(pct[pct > 5]), color = "grey") +
  ggplot2::theme_bw()
# 1.4. Heatmap of the first 12 PCs
Seurat::DimHeatmap(monocytes_macrophages, dims=1:12, cells=500, balanced=TRUE)
# 1.5. Determine the K-nearest neighbor graph
monocytes_macrophages = Seurat::FindNeighbors(monocytes_macrophages, dims=1:elbow)
# 1.6. Find clusters:
monocytes_macrophages = Seurat::FindClusters(monocytes_macrophages, resolution=seq(0.1, 1, by=.1))
# 1.7. UMAP:
monocytes_macrophages = Seurat::RunUMAP(monocytes_macrophages, dims=1:elbow, reduction='pca')
invisible(gc())
# 1.8. Save object
SeuratDisk::SaveH5Seurat(monocytes_macrophages, paste(project_dir, '2_annotation/results_Myeloid/datasets/monocytes_macrophages.h5Seurat', sep='/'))


# 2. Choose best resolution using clustree:
library(ggraph)
clust_tree = clustree::clustree(monocytes_macrophages, prefix='integrated_snn_res.')
clust_tree
# 2.1. Visualize different resolutions (first resolution was chosen to be the best for now):
clusters_01_0 = Seurat::DimPlot(monocytes_macrophages, reduction="umap", group.by='integrated_snn_res.0.1', label=TRUE, label.size=6, pt.size=.2)
clusters_02_0 = Seurat::DimPlot(monocytes_macrophages, reduction="umap", group.by='integrated_snn_res.0.2', label=TRUE, label.size=6, pt.size=.2)
clusters_03_0 = Seurat::DimPlot(monocytes_macrophages, reduction="umap", group.by='integrated_snn_res.0.3', label=TRUE, label.size=6, pt.size=.2)
clusters_04_0 = Seurat::DimPlot(monocytes_macrophages, reduction="umap", group.by='integrated_snn_res.0.4', label=TRUE, label.size=6, pt.size=.2)
((clusters_01_0 | clusters_02_0) / (clusters_03_0 | clusters_04_0)) | clust_tree
# 2.2. Feature Plots to assess known gene markers:
feature_plots(monocytes_macrophages, c('rna_IL1B', 'rna_IL6', 'rna_S100A8', 'rna_S100A9',
                                       'rna_CD163', 'rna_SEPP1', 'rna_APOE', 'rna_MAF',
                                       'rna_MKI67', 'rna_KIAA0101', 'rna_STMN1',
                                       'rna_SPP1', 'rna_FCN1'),
              ncol=3, with_dimplot=TRUE, dimplot_group='integrated_snn_res.0.3') # Pro-, anti- inflammatory, proliferative, SPP1+, FCN1+
violin_plots(monocytes_macrophages, c('rna_IL1B', 'rna_IL6', 'rna_S100A8', 'rna_S100A9',
                                      'rna_CD163', 'rna_SEPP1', 'rna_APOE', 'rna_MAF',
                                      'rna_MKI67', 'rna_KIAA0101', 'rna_STMN1',
                                      'rna_SPP1', 'rna_FCN1'),
             ncol=5, with_dimplot=TRUE, group.by='integrated_snn_res.0.3') # Pro-, anti- inflammatory, proliferative, SPP1+, FCN1+
# Resolution 0.3

# 3. Find markers
# 3.1. Calculate all markers
Seurat::Idents(monocytes_macrophages) = 'integrated_snn_res.0.3'
monocytes_macrophages_markers = Seurat::FindAllMarkers(monocytes_macrophages, assay='RNA', slot='data', logfc.threshold=.5, only.pos=TRUE)
View(monocytes_macrophages_markers)
write.csv(monocytes_macrophages_markers, paste(project_dir, '2_annotation/results_Myeloid/markers/markers_monoMacroLineage_res03.csv', sep='/'))
# 3.2. Get top 20 markers for each sub-cluster
monocytes_macrophages_markers_top20 = c()
for(clust in as.character(unique(monocytes_macrophages_markers$cluster))){
  clust_markers = monocytes_macrophages_markers[monocytes_macrophages_markers$cluster==clust,]
  n_markers = dim(clust_markers)[1]
  order_by_avg_log2FC = order(clust_markers$avg_log2FC, decreasing=T)
  if(n_markers > 20) top20_clust_markers = clust_markers[order_by_avg_log2FC[1:20],]
  else top20_clust_markers = clust_markers[order_by_avg_log2FC,]
  monocytes_macrophages_markers_top20 = rbind(monocytes_macrophages_markers_top20, top20_clust_markers)
}
View(monocytes_macrophages_markers_top20)
write.csv(monocytes_macrophages_markers_top20, paste(project_dir, '2_annotation/results_Myeloid/markers/markers_monoMacroLineage_res03_top20.csv', sep='/'))


# 4. Annotation for the monocyte/macrophage lineage:
# Sub-clusters 0 + 6      --> SPP1+ anti-inflammatory macro/mono
# Sub-cluster  1          --> pro-inflammatory macro/mono
# Sub-clusters 2 + 3 + 5  --> anti-inflammatory macro/mono
# Sub-cluster  4          --> FCN1+ pro-inflammatory macro/mono
# 4.1 Get cells in each of these annotations and annotate directly the big Myeloid dataset
spp1_cells = rownames(monocytes_macrophages@meta.data)[monocytes_macrophages$integrated_snn_res.0.3%in%c('0','6')]
pro_cells = rownames(monocytes_macrophages@meta.data)[monocytes_macrophages$integrated_snn_res.0.3=='1']
anti_cells = rownames(monocytes_macrophages@meta.data)[monocytes_macrophages$integrated_snn_res.0.3%in%c('2','3','5')]
fcn1_cells = rownames(monocytes_macrophages@meta.data)[monocytes_macrophages$integrated_snn_res.0.3=='4']




# -----
# - Final Annotations
# -----

# 1. Annotation level 2:
# 1.1. Annotate of the cells not in the monocyte/macrophage lineage:
# Sub-cluster  0   --> cDCs
# Sub-cluster  5   --> Mast cells
# Sub-cluster  14  --> pDCs
# Sub-cluster  13  --> Unknown
Myeloid[['Annotation_Level_2']] = rep('', length(Myeloid$integrated_snn_res.0.8))
cDCs_cells = rownames(Myeloid@meta.data)[Myeloid$integrated_snn_res.0.8=='0']
Myeloid@meta.data[cDCs_cells, 'Annotation_Level_2'] = 'cDCs'
mast_cells = rownames(Myeloid@meta.data)[Myeloid$integrated_snn_res.0.8=='5']
Myeloid@meta.data[mast_cells, 'Annotation_Level_2'] = 'Mast cells'
pDCs_cells = rownames(Myeloid@meta.data)[Myeloid$integrated_snn_res.0.8=='14']
Myeloid@meta.data[pDCs_cells, 'Annotation_Level_2'] = 'pDCs'
unknown_cells = rownames(Myeloid@meta.data)[Myeloid$integrated_snn_res.0.8=='13']
Myeloid@meta.data[unknown_cells, 'Annotation_Level_2'] = 'Unknown'

# 1.2. Annotate the cells in the monocyte/macrophage lineage:
Myeloid@meta.data[spp1_cells, 'Annotation_Level_2'] = 'SPP1+ Anti-inflammatory macro/mono'
Myeloid@meta.data[pro_cells, 'Annotation_Level_2'] = 'Other pro-inflammatory macro/mono'
Myeloid@meta.data[anti_cells, 'Annotation_Level_2'] = 'Other anti-inflammatory macro/mono'
Myeloid@meta.data[fcn1_cells, 'Annotation_Level_2'] = 'FCN1+ Pro-inflammatory macro/mono'


# 2. Annotation level 1:
Myeloid[['Annotation_Level_1']] = as.character(Myeloid$Annotation_Level_2)
Myeloid@meta.data[Myeloid$Annotation_Level_2%in%c('SPP1+ Anti-inflammatory macro/mono', 'Other anti-inflammatory macro/mono'),
                  'Annotation_Level_1'] = 'Anti-inflammatory macro/mono'
Myeloid@meta.data[Myeloid$Annotation_Level_2%in%c('FCN1+ Pro-inflammatory macro/mono', 'Other pro-inflammatory macro/mono'),
                  'Annotation_Level_1'] = 'Pro-inflammatory macro/mono'


# 3. Annotation level 0:
Myeloid[['Annotation_Level_0']] = factor(rep('Myeloid cells', dim(Myeloid@meta.data)[1]))


# 4. Remove [integrated...] metadata variables
Myeloid@meta.data = Myeloid@meta.data[, !colnames(Myeloid@meta.data) %in%
                                        c(grep('^in', colnames(Myeloid@meta.data), value=T), 'seurat_clusters')]


# 5. Visualize annotations:
Seurat::DimPlot(Myeloid, reduction="umap", group.by='Annotation_Level_2', label=TRUE, label.size=3, pt.size=.5, repel=T) +
  ggplot2::theme_minimal()
Seurat::DimPlot(Myeloid, reduction="umap", group.by='Annotation_Level_1', label=TRUE, label.size=3, pt.size=.5, repel=T) +
  ggplot2::theme_minimal()

# 6. Save Seurat object
SeuratDisk::SaveH5Seurat(Myeloid, paste(project_dir, '2_annotation/results_Myeloid/datasets/Myeloid_finalAnnots.h5Seurat', sep='/'))


# 7. Get Markers between the final annotations:
# 7.1. Annotation Level 2
Seurat::Idents(Myeloid) = 'Annotation_Level_2'
Myeloid_Level_2_markers = Seurat::FindAllMarkers(Myeloid, assay='RNA', slot='data', logfc.threshold=.8, min.pct = 0.3, only.pos=TRUE)
View(Myeloid_Level_2_markers)
write.csv(Myeloid_Level_2_markers, paste(project_dir, '2_annotation/results_Myeloid/markers/markers_Level_2.csv', sep='/'))
invisible(gc())
# 7.2. Annotation Level 1
Seurat::Idents(Myeloid) = 'Annotation_Level_1'
Myeloid_Level_1_markers = Seurat::FindAllMarkers(Myeloid, assay='RNA', slot='data', logfc.threshold=.8, min.pct = 0.3, only.pos=TRUE)
View(Myeloid_Level_1_markers)
write.csv(Myeloid_Level_1_markers, paste(project_dir, '2_annotation/results_Myeloid/markers/markers_Level_1.csv', sep='/'))
invisible(gc())

