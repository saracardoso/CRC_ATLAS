project_dir = '/home/scardoso/Documents/PhD/CRC_ATLAS'

# Load subset of epithelial cells:
Stromal = SeuratDisk::LoadH5Seurat('/home/scardoso/Documents/PhD/CRC_ATLAS/2_annotation/data/stromal_cells.h5Seurat')



# -----
# - Integrate epithelial cells
# -----

# Separate epithelial cells by dataset:
Stromal = Seurat::SplitObject(Stromal, split.by='dataset')
# SCTransform data:
for(dts_name in names(Stromal)) {
  Stromal[[dts_name]] <- Seurat::SCTransform(Stromal[[dts_name]], vars.to.regress = c('CC.Difference', 'percent.mitochondrial_RNA'))
  gc()
}
gc()
# Select integration features and prepare for integration:
features = Seurat::SelectIntegrationFeatures(Stromal, nfeatures = 3000)
Stromal = Seurat::PrepSCTIntegration(Stromal, anchor.features=features)
# Find integration achors and integrate the data:
anchors = Seurat::FindIntegrationAnchors(Stromal, normalization.method='SCT', anchor.features=features)
gc()
Stromal = Seurat::IntegrateData(anchorset=anchors, normalization.method="SCT")
gc()



# -----
# - Prepare data to find clusters
# -----

# Run PCA:
Stromal = Seurat::RunPCA(Stromal, assay='integrated')
# Find the number of PCs to use:
#pct = Stromal[["pca"]]@stdev /sum(Stromal[["pca"]]@stdev) * 100
#cumu = cumsum(pct)
#co1 = which(cumu > 90 & pct < 5)[1]
#co2 = sort(which((pct[1:length(pct) - 1] - pct[2:length(pct)]) > 0.1), decreasing = T)[1] + 1
#elbow = min(co1, co2) # 20

#plot_df = data.frame(pct=pct, cumu=cumu, rank=1:length(pct))
#ggplot2::ggplot(plot_df, ggplot2::aes(cumu, pct, label = rank, color = rank > elbow)) + 
#  ggplot2::geom_text() + 
#  ggplot2::geom_vline(xintercept = 90, color = "grey") + 
#  ggplot2::geom_hline(yintercept = min(pct[pct > 5]), color = "grey") +
#  ggplot2::theme_bw()

#Seurat::DimHeatmap(Stromal, dims=1:21, cells=500, balanced=TRUE)



# -----
# - Find Clusters
# -----

# Determine the K-nearest neighbor graph
Stromal = Seurat::FindNeighbors(Stromal, dims=1:50)#dims=1:elbow)

# Find clusters:
Stromal = Seurat::FindClusters(Stromal, resolution=c(0.4, 0.6, 1))

# UMAP:
Stromal = Seurat::RunUMAP(Stromal, dims=1:50, reduction='pca')#dims=1:elbow, reduction='pca')

# Save object with clusters:
SeuratDisk::SaveH5Seurat(Stromal, paste(project_dir, '2_annotation/data/Stromal_integratedClusters.h5Seurat', sep='/'))



# -----
# - Visualize Clusters
# -----

# Cluster from the three resolutions:
clusters_04 = Seurat::DimPlot(Stromal, reduction="umap", group.by='integrated_snn_res.0.4', label=TRUE, label.size=6)
clusters_06 = Seurat::DimPlot(Stromal, reduction="umap", group.by='integrated_snn_res.0.6', label=TRUE, label.size=6)
clusters_1 = Seurat::DimPlot(Stromal, reduction="umap", group.by='integrated_snn_res.1', label=TRUE, label.size=6)
(clusters_04 | clusters_06) / (clusters_1 | ggplot2::ggplot())

# Choose resolution 0.6
res = 'integrated_snn_res.0.6'

# Metadata:
state = Seurat::DimPlot(Stromal, reduction="umap", group.by='state', label=FALSE, label.size=6)
patients = Seurat::DimPlot(Stromal, reduction="umap", group.by='patient', label=FALSE, label.size=6, pt.size=.1)
datasets = Seurat::DimPlot(Stromal, reduction="umap", group.by='dataset', label=FALSE, label.size=6, pt.size=.1)
(clusters_06 | patients) / (state | datasets)

# Normal/ tumour/ healthy distribution across clusters:
clusters = c()
percs = c()
fill = c()
for(cluster in unique(as.character(Stromal@meta.data[,res]))){
  x = table(Stromal$state[Stromal@meta.data[,res]==cluster])[c('Tumor', 'Normal', 'Healthy')]
  x[is.na(x)] = 0
  y = sum(x)
  percs = c(percs, x/y)
  fill = c(fill, c('Tumor', 'Normal', 'Healthy'))
  clusters = c(clusters, rep(cluster, 3))
}
dist_nt = data.frame(percs=percs, fill=fill, clusters=clusters)
dist_nt$clusters = factor(dist_nt$clusters,
                          levels=as.character(0:(length(unique(Stromal@meta.data[,res]))-1)))
dist_state = ggplot2::ggplot(dist_nt, mapping = ggplot2::aes(clusters, percs, fill=fill)) +
  ggplot2::geom_bar(position='stack', stat='identity')
dist_state

# Expression markers - endothelial cells:
featurePlot_ec = Seurat::FeaturePlot(Stromal, features=c('rna_VWF', 'rna_PLVAP', 'rna_CDH5', 'rna_ENG'), min.cutoff = 'q1')
vlnPlot_ec = Seurat::VlnPlot(Stromal, features=c('rna_VWF', 'rna_PLVAP', 'rna_CDH5', 'rna_ENG'),
                                 pt.size = 0, ncol=2, group.by = res)
clusters_06 | ( featurePlot_ec / vlnPlot_ec )

# Expression markers - lymphatic endothelial cells:
featurePlot_ec = Seurat::FeaturePlot(Stromal, features=c('rna_LYVE1', 'rna_PROX1'), min.cutoff = 'q1')
vlnPlot_ec = Seurat::VlnPlot(Stromal, features=c('rna_LYVE1', 'rna_PROX1'),
                             pt.size = 0, ncol=2, group.by = res)
clusters_06 | ( featurePlot_ec / vlnPlot_ec )

# Expression markers - enteric glia:
featurePlot_eg = Seurat::FeaturePlot(Stromal, features=c('rna_SOX10', 'rna_S100B'), min.cutoff = 'q1')
vlnPlot_eg = Seurat::VlnPlot(Stromal, features=c('rna_SOX10', 'rna_S100B'),
                             pt.size = 0, ncol=2, group.by = res)
clusters_06 | ( featurePlot_eg / vlnPlot_eg )

# Expression markers - VSMC (vascular smooth muscle cells):
featurePlot_vsmc = Seurat::FeaturePlot(Stromal, features=c('rna_SYNPO2', 'rna_CRYAB', 'rna_CNN1', 'rna_DES'), min.cutoff = 'q1')
vlnPlot_vsmc = Seurat::VlnPlot(Stromal, features=c('rna_SYNPO2', 'rna_CRYAB', 'rna_CNN1', 'rna_DES'),
                             pt.size = 0, ncol=2, group.by = res)
clusters_06 | ( featurePlot_vsmc / vlnPlot_vsmc )

# Expression markers - pericytes
featurePlot_p = Seurat::FeaturePlot(Stromal, features=c('rna_RGS5', 'rna_CSPG4', 'rna_ABCC9', 'rna_KCNJ8'), min.cutoff = 'q1')
vlnPlot_p = Seurat::VlnPlot(Stromal, features=c('rna_RGS5', 'rna_CSPG4', 'rna_ABCC9', 'rna_KCNJ8'),
                               pt.size = 0, ncol=2, group.by = res)
clusters_06 | ( featurePlot_p / vlnPlot_p )

# Expression markers - fibroblasts
featurePlot_f = Seurat::FeaturePlot(Stromal, features=c('rna_COL1A1', 'rna_COL1A2', 'rna_COL6A1', 'rna_COL6A2', 'rna_COL3A1', 'rna_DCN'),
                                    min.cutoff = 'q1')
vlnPlot_f = Seurat::VlnPlot(Stromal, features=c('rna_COL1A1', 'rna_COL1A2', 'rna_COL6A1', 'rna_COL6A2', 'rna_COL3A1', 'rna_DCN'),
                            pt.size = 0, ncol=2, group.by = res)
clusters_06 | ( featurePlot_f / vlnPlot_f )

# Expression markers - myofibroblasts
featurePlot_mf = Seurat::FeaturePlot(Stromal, features=c('rna_TAGLN', 'rna_ACTA2', 'rna_ACTG2', 'rna_MYH11', 'rna_MYLK'),
                                    min.cutoff = 'q1')
vlnPlot_mf = Seurat::VlnPlot(Stromal, features=c('rna_TAGLN', 'rna_ACTA2', 'rna_ACTG2', 'rna_MYH11', 'rna_MYLK'),
                            pt.size = 0, ncol=2, group.by = res)
clusters_06 | ( featurePlot_mf / vlnPlot_mf )

# Expression markers - CAFs
featurePlot_caf = Seurat::FeaturePlot(Stromal, features=c('rna_THY1', 'rna_FAP'), ncol=2)
vlnPlot_caf = Seurat::VlnPlot(Stromal, features=c('rna_THY1', 'rna_FAP'), pt.size = 0, ncol=2, group.by = res)
clusters_06 | ( featurePlot_caf / vlnPlot_caf )




# -----
# - Calculate similarity between clusters
# -----

# Calculate average and median expression profiles of each cluster:
n_genes = dim(Stromal@assays$RNA@data)[1]
n_clusters = length(unique(Stromal@meta.data[,res]))
average_profiles = matrix(rep(0, n_genes * n_clusters), nrow = n_genes)
colnames(average_profiles) = as.character(c(0:(n_clusters-1)))
rownames(average_profiles) = rownames(Stromal@assays$RNA@data)
for(clust in as.character(c(0:(n_clusters-1)))){
  message('Cluster:', clust)
  clust_cells = rownames(Stromal@meta.data)[Stromal@meta.data[,res] == clust]
  clust_matrix = as.matrix(Stromal@assays$RNA@data[,clust_cells])
  average_profiles[,clust] = rowMeans(clust_matrix)
  invisible(gc())
}

# Calculate distance matrices between clusters:
dist_average_clusts = dist(t(average_profiles), method='euclidean')

# Visualize distance matrices:
dist_average_clusts_matrix = as.matrix(dist_average_clusts)
#dist_average_clusts_matrix[upper.tri(dist_average_clusts_matrix)] <- NA
pheatmap::pheatmap(dist_average_clusts_matrix, cluster_rows=F, cluster_cols=F, na_col="white", display_numbers=TRUE,
                   main='Clusters distance based on average profiles')

# Create dendogram:
hclust_average_clusters = hclust(dist_average_clusts, method='complete')
plot(hclust_average_clusters, main='Clustering based on average profiles')




# -----
# - CAFs vs Fibroblasts
# -----

Stromal[['state2']] = as.character(Stromal@meta.data[,'state'])
Stromal$state2[Stromal$state2=='Healthy'] = 'Normal'

Stromal[['temp_clusters0']] = as.character(Stromal@meta.data[,res])
Stromal$temp_clusters0[Stromal$temp_clusters0=='1'] = paste(Stromal$temp_clusters0[Stromal$temp_clusters0=='1'],
                                                            Stromal$state2[Stromal$temp_clusters0=='1'], sep='_')
Stromal$temp_clusters0[Stromal$temp_clusters0=='14'] = paste(Stromal$temp_clusters0[Stromal$temp_clusters0=='14'],
                                                             Stromal$state2[Stromal$temp_clusters0=='14'], sep='_')

Seurat::VlnPlot(subset(Stromal, cells=rownames(Stromal@meta.data)[!Stromal$integrated_snn_res.0.6%in%c('5', '6', '7', '8', '9', '13', '16',
                                                                                                       '17', '18', '21')]),
                features=c('rna_THY1', 'rna_FAP'), pt.size = 0, ncol=1, group.by = 'temp_clusters0')

# Calculate average expression profiles of each cluster:
clusts = unique(Stromal$temp_clusters0)
n_genes = dim(Stromal@assays$RNA@data)[1]
n_clusters = length(clusts)
average_profiles = matrix(rep(0, n_genes * n_clusters), nrow = n_genes)
colnames(average_profiles) = clusts
rownames(average_profiles) = rownames(Stromal@assays$RNA@data)
for(clust in clusts){
  clust_cells = rownames(Stromal@meta.data)[Stromal@meta.data[,'temp_clusters0'] == clust]
  clust_matrix = as.matrix(Stromal@assays$RNA@data[,clust_cells])
  average_profiles[,clust] = rowMeans(clust_matrix)
  invisible(gc())
}

# Calculate distance matrices between clusters:
dist_average_clusts = dist(t(average_profiles), method='euclidean')

# Visualize distance matrices:
dist_average_clusts_matrix = as.matrix(dist_average_clusts)
pheatmap::pheatmap(dist_average_clusts_matrix, cluster_rows=F, cluster_cols=F, na_col="white", display_numbers=TRUE,
                   main='Clusters distance based on average profiles')




# -----
# - Find markers
# -----

Seurat::Idents(Stromal) = res
Stromal_markers = Seurat::FindAllMarkers(Stromal, assay='RNA', slot='data', logfc.threshold=1, only.pos=TRUE)
View(Stromal_markers)
write.csv(Stromal_markers, paste(project_dir, '2_annotation/markers/Stromal/Stromal_res06.csv', sep='/'))

# Get top 10 markers for each cluster:
Stromal_markers_top20 = dplyr::top_n(dplyr::group_by(Stromal_markers, cluster), n=40, wt=avg_log2FC)
View(Stromal_markers_top20)
write.csv(Stromal_markers_top20, paste(project_dir, '2_annotation/markers/Stromal/Stromal_res06_top20.csv', sep='/'))

# Merge clusters together, based on previous markers and similarity measures:
Stromal[['temp_clusters']] = as.character(Stromal@meta.data[,'temp_clusters0'])
Stromal$temp_clusters[Stromal$temp_clusters%in%c('5', '6', '7', '17', '21')] = '5_6_7_17_21'
Stromal$temp_clusters[Stromal$temp_clusters%in%c('1_Tumor', '14_Tumor')] = '1Tumor_14Tumor'
Stromal$temp_clusters[Stromal$temp_clusters%in%c('1_Normal', '14_Normal', '0', '2', '3', '4', '10', '11', '12', '15', '19', '20')] =
  '1Normal_14Normal_0_2_3_4_10_11_12_15_19_20'
Seurat::DimPlot(Stromal, reduction="umap", group.by=res, label=TRUE, label.size=4) + Seurat::NoLegend() |
  Seurat::DimPlot(Stromal, reduction="umap", group.by='temp_clusters', label=TRUE, label.size=4) + Seurat::NoLegend()

# Calculate markers for new clusters:
Seurat::Idents(Stromal) = 'temp_clusters'
Stromal_markers = Seurat::FindAllMarkers(Stromal, assay='RNA', slot='data', logfc.threshold=.8, only.pos=TRUE)
View(Stromal_markers)
write.csv(Stromal_markers, paste(project_dir, '2_annotation/markers/Stromal/Stromal_tempclusters.csv', sep='/'))

# Get top 20 markers for each cluster:
Stromal_markers_top20 = dplyr::top_n(dplyr::group_by(Stromal_markers[Stromal_markers$p_val_adj<0.05,], cluster), n=20, wt=avg_log2FC)
View(Stromal_markers_top20)
write.csv(Stromal_markers_top20, paste(project_dir, '2_annotation/markers/Stromal/Stromal_tempclusters_top20.csv', sep='/'))
