project_dir = '/home/scardoso/Documents/PhD/CRC_ATLAS'
source(paste(project_dir, 'utils/modified_plots.R', sep='/'))

# Load subset of epithelial cells:
Bcells = SeuratDisk::LoadH5Seurat(paste(project_dir, '/2_annotation/results_Bcells/datasets/Bcells.h5Seurat', sep='/'))





# -----
# - Find clusters
# -----


# 1. Run PCA:
Bcells = Seurat::RunPCA(Bcells, assay='integrated')
invisible(gc())
Seurat::DefaultAssay(Bcells)='integrated'

# 2. Choose number of PCs to use (where the elbow occurs):
pct = Bcells[["pca"]]@stdev /sum(Bcells[["pca"]]@stdev) * 100
cumu = cumsum(pct)
co1 = which(cumu > 90 & pct < 5)[1]
co2 = sort(which((pct[1:length(pct) - 1] - pct[2:length(pct)]) > 0.1), decreasing = T)[1] + 1
# 2.1 Number of PCs to use
elbow = min(co1, co2) # 10
# 2.2. Visual representation of the PCs to use
plot_df = data.frame(pct=pct, cumu=cumu, rank=1:length(pct))
ggplot2::ggplot(plot_df, ggplot2::aes(cumu, pct, label = rank, color = rank > elbow)) + 
  ggplot2::geom_text() + 
  ggplot2::geom_vline(xintercept = 90, color = "grey") + 
  ggplot2::geom_hline(yintercept = min(pct[pct > 5]), color = "grey") +
  ggplot2::theme_bw()
# 2.3. Heatmap of the first 12 PCs
Seurat::DimHeatmap(Bcells, dims=1:12, cells=500, balanced=TRUE)


# 3. Find clusters
# 3.1. Determine the K-nearest neighbor graph
Bcells = Seurat::FindNeighbors(Bcells, dims=1:elbow)
# 3.2. Find clusters:
Bcells = Seurat::FindClusters(Bcells, resolution=seq(0.1, 1, by=.1))
# 3.3. UMAP:
Bcells = Seurat::RunUMAP(Bcells, dims=1:elbow, reduction='pca')
invisible(gc())


# 4. Choose best resolution:
library(ggraph)
clust_tree = clustree::clustree(Bcells, prefix='integrated_snn_res.')
clust_tree
# 4.1. Visualize different resolutions (first resolution was chosen to be the best for now):
clusters_01 = Seurat::DimPlot(Bcells, reduction="umap", group.by='integrated_snn_res.0.1', label=TRUE, label.size=6, pt.size=.2)
clusters_04 = Seurat::DimPlot(Bcells, reduction="umap", group.by='integrated_snn_res.0.4', label=TRUE, label.size=6, pt.size=.2)
clusters_05 = Seurat::DimPlot(Bcells, reduction="umap", group.by='integrated_snn_res.0.5', label=TRUE, label.size=6, pt.size=.2)
clusters_07 = Seurat::DimPlot(Bcells, reduction="umap", group.by='integrated_snn_res.0.7', label=TRUE, label.size=6, pt.size=.2)
((clusters_01 | clusters_04) / (clusters_05 | clusters_07)) | clust_tree
# 4.2. Plot Metrics
metrics =  c("nUMI", "nGene", "S.Score", "G2M.Score", "percent.mitochondrial_RNA")
feature_plots(Bcells, metrics, ncol=3, with_dimplot=TRUE, dimplot_group='integrated_snn_res.0.7')
# 4.3. Plot metadata
state = Seurat::DimPlot(Bcells, reduction="umap", group.by='state', label=FALSE, label.size=6) + Seurat::NoLegend()
patients = Seurat::DimPlot(Bcells, reduction="umap", group.by='patient', label=FALSE, label.size=6, pt.size=.1) + Seurat::NoLegend()
datasets = Seurat::DimPlot(Bcells, reduction="umap", group.by='dataset', label=FALSE, label.size=6, pt.size=.1) + Seurat::NoLegend()
(clusters_07 | patients) / (state | datasets)
# 4.2. Feature Plots:
feature_plots(Bcells, c('rna_MZB1', 'rna_IL2RA', 'rna_MS4A1', 'rna_CD27', 'rna_IGHD'),
              ncol=3, with_dimplot=TRUE, dimplot_group='integrated_snn_res.0.7')
# 4.4. Resolution 0.7 is chosen

# 5. Save object:
SeuratDisk::SaveH5Seurat(Bcells, paste(project_dir, '2_annotation/results_Bcells/datasets/Bcells_clustered.h5Seurat', sep='/'))





# -----
# - Annotate Clusters
# -----


# 1. IGs expression:
ighm = 'IGHM'
ighd = 'IGHD'
igha = grep('^IGHA', rownames(Bcells@assays$RNA@data), value=T)
ighg = grep('^IGHG', rownames(Bcells@assays$RNA@data), value=T)
feature_plots(Bcells, features=paste('rna', c(igha, ighd, ighg, ighm), sep='_'),
              ncol=3, with_dimplot=TRUE, dimplot_group='integrated_snn_res.0.7')
violin_plots(Bcells, features=paste('rna', c(igha, ighd, ighg, ighm), sep='_'),
             ncol=3, with_dimplot=TRUE, group.by='integrated_snn_res.0.7')


# 2. Kappa vs lamda expression:
kappa = grep('^IGK', rownames(Bcells@assays$RNA@data), value=T)
kappa_expression = Matrix::colSums(Seurat::GetAssayData(Bcells, assay='RNA', slot='data')[kappa,])
Bcells[['kappa_expression']] = kappa_expression
lambda = grep('^IGL', rownames(Bcells@assays$RNA@data), value=T)
lambda_expression = Matrix::colSums(Seurat::GetAssayData(Bcells, assay='RNA', slot='data')[lambda[2:44],])
Bcells[['lambda_expression']] = lambda_expression
feature_plots(Bcells, features=c('kappa_expression', 'lambda_expression'), ncol=2, with_dimplot=TRUE, dimplot_group='integrated_snn_res.0.7')
violin_plots(Bcells, features=c('kappa_expression', 'lambda_expression'), ncol=2, with_dimplot=TRUE, group.by='integrated_snn_res.0.7')


# 3. Memory vs activated vs plasma vs naive
feature_plots(Bcells, features=c('rna_MZB1', 'rna_MS4A1', 'rna_CD27', 'rna_CXCR4', 'rna_NR4A2', 'rna_RGS13'),
              ncol=4, with_dimplot=TRUE, dimplot_group='integrated_snn_res.0.7')
violin_plots(Bcells, features=c('rna_MZB1', 'rna_MS4A1', 'rna_CD27', 'rna_CXCR4', 'rna_NR4A2', 'rna_RGS13'),
             ncol=4, with_dimplot=TRUE, group.by='integrated_snn_res.0.7')


# 4. Calculate similarity between clusters
# 4.1. Calculate average expression profiles of each cluster:
Bcells_average_profiles = Seurat::AverageExpression(Bcells, group.by = 'integrated_snn_res.0.7', assays='RNA', slor='data')$RNA
# 4.2. Calculate distance matrices between clusters:
Bcells_dist_average_clusts = dist(t(Bcells_average_profiles), method='euclidean')
# 4.3. Visualize distance matrices:
dist_average_clusts_matrix = as.matrix(Bcells_dist_average_clusts)
#dist_average_clusts_matrix[upper.tri(dist_average_clusts_matrix)] <- NA
pheatmap::pheatmap(dist_average_clusts_matrix, cluster_rows=F, cluster_cols=F, na_col="white", display_numbers=TRUE,
                   main='Clusters distance based on average profiles')
# 4.4. Create dendogram:
hclust_average_clusters = hclust(Bcells_dist_average_clusts, method='complete')
plot(hclust_average_clusters, main='Clustering based on average profiles')


# 5. Find markers
# 5.1. Calculate all markers
Seurat::Idents(Bcells) = 'integrated_snn_res.0.7'
Bcells_markers = Seurat::FindAllMarkers(Bcells, assay='RNA', slot='data', logfc.threshold=.5, only.pos=TRUE)
View(Bcells_markers)
write.csv(Bcells_markers, paste(project_dir, '2_annotation/results_Bcells/markers/markers_res07.csv', sep='/'))
# 5.2. Get top 20 markers for each sub-cluster
Bcells_markers_top20 = c()
for(clust in as.character(unique(Bcells_markers$cluster))){
  clust_markers = Bcells_markers[Bcells_markers$cluster==clust,]
  n_markers = dim(clust_markers)[1]
  order_by_avg_log2FC = order(clust_markers$avg_log2FC, decreasing=T)
  if(n_markers > 20) top20_clust_markers = clust_markers[order_by_avg_log2FC[1:20],]
  else top20_clust_markers = clust_markers[order_by_avg_log2FC,]
  Bcells_markers_top20 = rbind(Bcells_markers_top20, top20_clust_markers)
}
View(Bcells_markers_top20)
write.csv(Bcells_markers_top20, paste(project_dir, '2_annotation/results_Bcells/markers/markers_res07_top20.csv', sep='/'))


# 6. The following clusters will be merged, based on previous markers and similarity measures, and new markers calculated:
# 10 + 11
# 3 + 4 + 6 + 8 + 15
# 2 + 12 + 13
# 0 + 1 + 5 + 7 + 14
# 6.1. Merge Clusters:
Bcells[['temp_clusters']] = as.character(Bcells@meta.data[,'integrated_snn_res.0.7'])
Bcells$temp_clusters[Bcells$temp_clusters%in%c('0', '1', '5', '7', '14')] = '0_1_5_7_14'
Bcells$temp_clusters[Bcells$temp_clusters%in%c('2', '12', '13')] = '2_12_14'
Bcells$temp_clusters[Bcells$temp_clusters%in%c('3', '4', '6', '8', '15')] = '3_4_6_8_15'
Bcells$temp_clusters[Bcells$temp_clusters%in%c('10', '11')] = '10_11'
Seurat::DimPlot(Bcells, reduction="umap", group.by='integrated_snn_res.0.7', label=TRUE, label.size=6) |
  Seurat::DimPlot(Bcells, reduction="umap", group.by='temp_clusters', label=TRUE, label.size=6)
# 6.3. Calculate markers for new clusters:
Seurat::Idents(Bcells) = 'temp_clusters'
Bcells_markers = Seurat::FindAllMarkers(Bcells, assay='RNA', slot='data', logfc.threshold=.8, only.pos=TRUE)
View(Bcells_markers)
write.csv(Bcells_markers, paste(project_dir, '2_annotation/markers/Bcells/markers_tempclusters.csv', sep='/'))
# 6.4. Get top 20 markers for each cluster:
Bcells_markers_top20 = dplyr::top_n(dplyr::group_by(Bcells_markers[Bcells_markers$p_val_adj<0.05,], cluster), n=20, wt=avg_log2FC)
View(Bcells_markers_top20)
write.csv(Bcells_markers_top20, paste(project_dir, '2_annotation/results_Bcells/markers/markers_tempclusters_top20.csv', sep='/'))


# 7. Final feature plots
# 7.1. Ig expression
feature_plots(Bcells, features=paste('rna', c(igha, ighd, ighg, ighm), sep='_'),
              ncol=3, with_dimplot=TRUE, dimplot_group='temp_clusters')
# 7.2.  Memory vs plasma vs naive
feature_plots(Bcells, features=c('rna_MZB1', 'rna_MS4A1', 'rna_CD27', 'rna_CXCR4', 'rna_NR4A2'),
              ncol=3, with_dimplot=TRUE, dimplot_group='temp_clusters')
# 7.3. Proliferative B cells
feature_plots(Bcells, c('rna_STMN1', 'rna_ACTB', 'rna_RGS13', 'rna_MKI67', 'rna_PCNA', 'rna_MARCKSL1', 'rna_HMGN1', 'rna_HMGN2'),
              ncol=3, with_dimplot=TRUE, dimplot_group='temp_clusters')


# 8. Final violin plots
violin_plots(Bcells, features=c(paste('rna', c(igha, ighd, ighg, ighm), sep='_'),
                                'rna_MZB1', 'rna_MS4A1', 'rna_CD27', 'rna_CXCR4', 'rna_NR4A2',
                                'rna_STMN1', 'rna_ACTB', 'rna_RGS13', 'rna_MKI67', 'rna_PCNA', 'rna_MARCKSL1', 'rna_HMGN1', 'rna_HMGN2'),
             ncol=5, with_dimplot=TRUE, group.by='temp_clusters')


# 9. Annotate cells:
# Sub-clusters 0 + 1 + 5 + 7 + 14  --> Memory Bcells
# Sub-clusters 2 + 12 + 14         --> Naive Bcells
# Sub-clusters 3 + 4 + 6 + 8 + 15  --> IgA+ Plasma cells 1
# Sub-cluster  9                   --> IgG+ Plasma cells
# Sub-clusters 10 + 11             --> Proliferative Bcells
# Sub-cluster  16                  --> IgA+ Plasma cells 2 
# 9.1. Annotation level 2:
Bcells[['Annotation_Level_2']] = rep('', length(Bcells$temp_clusters))
memory_cells = rownames(Bcells@meta.data)[Bcells$temp_clusters=='0_1_5_7_14']
Bcells@meta.data[memory_cells, 'Annotation_Level_2'] = 'Memory Bcells'
naive_cells = rownames(Bcells@meta.data)[Bcells$temp_clusters=='2_12_14']
Bcells@meta.data[naive_cells, 'Annotation_Level_2'] = 'Naive Bcells'
iga1_cells = rownames(Bcells@meta.data)[Bcells$temp_clusters=='3_4_6_8_15']
Bcells@meta.data[iga1_cells, 'Annotation_Level_2'] = 'IgA+ Plasma cells 1'
igg_cells = rownames(Bcells@meta.data)[Bcells$temp_clusters=='9']
Bcells@meta.data[igg_cells, 'Annotation_Level_2'] = 'IgG+ Plasma cells'
prol_cells = rownames(Bcells@meta.data)[Bcells$temp_clusters=='10_11']
Bcells@meta.data[prol_cells, 'Annotation_Level_2'] = 'Proliferative Bcells'
iga2_cells = rownames(Bcells@meta.data)[Bcells$temp_clusters=='16']
Bcells@meta.data[iga2_cells, 'Annotation_Level_2'] = 'IgA+ Plasma cells 2'
# 9.2. Annotation level 1:
Bcells[['Annotation_Level_1']] = as.character(Bcells$Annotation_Level_2)
Bcells@meta.data[Bcells$Annotation_Level_2%in%c('IgA+ Plasma cells 1', 'IgA+ Plasma cells 2'), 'Annotation_Level_1'] = 'IgA+ Plasma cells'
# 9.3. Annotation level 0:
Bcells[['Annotation_Level_0']] = factor(rep('Bcells', dim(Bcells@meta.data)[1]))
# 9.4. Remove [integrated...] metadata variables
Bcells@meta.data = Bcells@meta.data[, !colnames(Bcells@meta.data) %in%
                                      c(grep('^in', colnames(Bcells@meta.data), value=T), 'seurat_clusters')]
# 9.5. Visualize annotations:
Seurat::DimPlot(Bcells, reduction="umap", group.by='Annotation_Level_2', label=TRUE, label.size=3, pt.size=.5, repel=T) +
  ggplot2::theme_minimal()
Seurat::DimPlot(Bcells, reduction="umap", group.by='Annotation_Level_1', label=TRUE, label.size=3, pt.size=.5, repel=T) +
  ggplot2::theme_minimal()
# 9.6. Save Seurat object
SeuratDisk::SaveH5Seurat(Bcells, paste(project_dir, '2_annotation/results_Bcells/datasets/Bcells_finalAnnots.h5Seurat', sep='/'))


# 10. Get Markers between the final annotations:
# 10.1. Annotation Level 2
Seurat::Idents(Bcells) = 'Annotation_Level_2'
Bcells_Level_2_markers = Seurat::FindAllMarkers(Bcells, assay='RNA', slot='data', logfc.threshold=.8, min.pct = 0.3, only.pos=TRUE)
View(Bcells_Level_2_markers)
write.csv(Bcells_Level_2_markers, paste(project_dir, '2_annotation/results_Bcells/markers/markers_Level_2.csv', sep='/'))
invisible(gc())
# 10.1. Annotation Level 1
Seurat::Idents(Bcells) = 'Annotation_Level_1'
Bcells_Level_1_markers = Seurat::FindAllMarkers(Bcells, assay='RNA', slot='data', logfc.threshold=.8, min.pct = 0.3, only.pos=TRUE)
View(Bcells_Level_1_markers)
write.csv(Bcells_Level_1_markers, paste(project_dir, '2_annotation/results_Bcells/markers/markers_Level_1.csv', sep='/'))
invisible(gc())


