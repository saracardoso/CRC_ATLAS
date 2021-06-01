project_dir = '/home/scardoso/Documents/PhD/CRC_ATLAS'

# Load subset of epithelial cells:
Bcells = SeuratDisk::LoadH5Seurat('/home/scardoso/Documents/PhD/CRC_ATLAS/2_annotation/data/Bcells.h5Seurat')



# -----
# - Integrate epithelial cells
# -----

# Separate epithelial cells by dataset:
Bcells = Seurat::SplitObject(Bcells, split.by='dataset')
# SCTransform data:
for(dts_name in names(Bcells)) {
  Bcells[[dts_name]] <- Seurat::SCTransform(Bcells[[dts_name]], vars.to.regress = c('CC.Difference', 'percent.mitochondrial_RNA'))
  gc()
}
gc()
# Select integration features and prepare for integration:
features = Seurat::SelectIntegrationFeatures(Bcells, nfeatures = 3000)
Bcells = Seurat::PrepSCTIntegration(Bcells, anchor.features=features)
# Find integration achors and integrate the data:
anchors = Seurat::FindIntegrationAnchors(Bcells, normalization.method='SCT', anchor.features=features)
gc()
Bcells = Seurat::IntegrateData(anchorset=anchors, normalization.method="SCT")
gc()



# -----
# - Prepare data to find clusters
# -----

# Run PCA:
Bcells = Seurat::RunPCA(Bcells, assay='integrated')
# Find the number of PCs to use:
#pct = Bcells[["pca"]]@stdev /sum(Bcells[["pca"]]@stdev) * 100
#cumu = cumsum(pct)
#co1 = which(cumu > 90 & pct < 5)[1]
#co2 = sort(which((pct[1:length(pct) - 1] - pct[2:length(pct)]) > 0.1), decreasing = T)[1] + 1
#elbow = min(co1, co2) # 12

#plot_df = data.frame(pct=pct, cumu=cumu, rank=1:length(pct))
#ggplot2::ggplot(plot_df, ggplot2::aes(cumu, pct, label = rank, color = rank > elbow)) + 
#  ggplot2::geom_text() + 
#  ggplot2::geom_vline(xintercept = 90, color = "grey") + 
#  ggplot2::geom_hline(yintercept = min(pct[pct > 5]), color = "grey") +
#  ggplot2::theme_bw()

#Seurat::DimHeatmap(Bcells, dims=1:15, cells=500, balanced=TRUE)



# -----
# - Find Clusters
# -----

# Determine the K-nearest neighbor graph
Bcells = Seurat::FindNeighbors(Bcells, dims=1:50)#dims=1:elbow)

# Find clusters:
Bcells = Seurat::FindClusters(Bcells, resolution=c(0.4, 0.6, 1))

# UMAP:
Bcells = Seurat::RunUMAP(Bcells, dims=1:50, reduction='pca')#dims=1:elbow, reduction='pca')

# Save object with clusters:
SeuratDisk::SaveH5Seurat(Bcells, paste(project_dir, '2_annotation/data/Bcells_integratedClusters.h5Seurat', sep='/'))



# -----
# - Visualize Clusters
# -----

# Cluster from the three resolutions:
clusters_04 = Seurat::DimPlot(Bcells, reduction="umap", group.by='integrated_snn_res.0.4', label=TRUE, label.size=6)
clusters_06 = Seurat::DimPlot(Bcells, reduction="umap", group.by='integrated_snn_res.0.6', label=TRUE, label.size=6)
clusters_1 = Seurat::DimPlot(Bcells, reduction="umap", group.by='integrated_snn_res.1', label=TRUE, label.size=6)
(clusters_04 | clusters_06) / clusters_1

# Choose resolution 0.6
res = 'integrated_snn_res.0.6'

# Metadata:
state = Seurat::DimPlot(Bcells, reduction="umap", group.by='state', label=FALSE, label.size=6)
patients = Seurat::DimPlot(Bcells, reduction="umap", group.by='patient', label=FALSE, label.size=6, pt.size=.1)
datasets = Seurat::DimPlot(Bcells, reduction="umap", group.by='dataset', label=FALSE, label.size=6, pt.size=.1)
(clusters_06 | patients) / (state | datasets)

# Normal/ tumour/ healthy distribution across clusters:
clusters = c()
percs = c()
fill = c()
for(cluster in unique(as.character(Bcells@meta.data[,res]))){
  x = table(Bcells$state[Bcells@meta.data[,res]==cluster])[c('Tumor', 'Normal', 'Healthy')]
  x[is.na(x)] = 0
  y = sum(x)
  percs = c(percs, x/y)
  fill = c(fill, c('Tumor', 'Normal', 'Healthy'))
  clusters = c(clusters, rep(cluster, 3))
}
dist_nt = data.frame(percs=percs, fill=fill, clusters=clusters)
dist_nt$clusters = factor(dist_nt$clusters,
                          levels=as.character(0:(length(unique(Bcells@meta.data[,res]))-1)))
dist_state = ggplot2::ggplot(dist_nt, mapping = ggplot2::aes(clusters, percs, fill=fill)) +
  ggplot2::geom_bar(position='stack', stat='identity')
dist_state

# IGs expression:
ighm = 'IGHM'
ighd = 'IGHD'
igha = grep('^IGHA', rownames(Bcells@assays$RNA@data), value=T)
ighg = grep('^IGHG', rownames(Bcells@assays$RNA@data), value=T)
clusters_06 | Seurat::FeaturePlot(Bcells, features=paste('rna', c(igha, ighd, ighg, ighm), sep='_'), min.cutoff = 'q1')
clusters_06 / Seurat::VlnPlot(Bcells, features=paste('rna', c(igha, ighd, ighg, ighm), sep='_'), ncol=5, pt.size=0, y.max=10, group.by=res)

# Kappa vs lamda expression:
kappa = grep('^IGK', rownames(Bcells@assays$RNA@data), value=T)
kappa_expression = colSums(Seurat::GetAssayData(Bcells, assay='RNA', slot='data')[kappa,])
Bcells[['kappa_expression']] = kappa_expression
lambda = grep('^IGL', rownames(Bcells@assays$RNA@data), value=T)
lambda_expression = colSums(Seurat::GetAssayData(Bcells, assay='RNA', slot='data')[lambda[2:44],])
Bcells[['lambda_expression']] = lambda_expression
clusters_06 | Seurat::FeaturePlot(Bcells, features=c('kappa_expression', 'lambda_expression'), ncol=1)
clusters_06 | Seurat::VlnPlot(Bcells, features=c('kappa_expression', 'lambda_expression'), ncol=1, pt.size=0, y.max=40, group.by=res)

# Memory vs activated vs plasma vs naive
clusters_06 / ( Seurat::FeaturePlot(Bcells, features=c('rna_MZB1', 'rna_IL2RA', 'rna_MS4A1', 'rna_CD27'), ncol=2) |
                  Seurat::VlnPlot(Bcells, features=c('rna_MZB1', 'rna_IL2RA', 'rna_MS4A1', 'rna_CD27'), ncol=2, pt.size=0, group.by=res))

clusters_06 / ( Seurat::FeaturePlot(Bcells, features=c('rna_MZB1', 'rna_CXCR4', 'rna_SDC1', 'rna_NR4A2'), ncol=2) |
                  Seurat::VlnPlot(Bcells, features=c('rna_MZB1', 'rna_CXCR4', 'rna_SDC1', 'rna_NR4A2'), ncol=2, pt.size=0, group.by=res))




# -----
# - Calculate similarity between clusters
# -----

# Calculate average and median expression profiles of each cluster:
n_genes = dim(Bcells@assays$RNA@data)[1]
n_clusters = length(unique(Bcells@meta.data[,res]))
average_profiles = matrix(rep(0, n_genes * n_clusters), nrow = n_genes)
colnames(average_profiles) = as.character(c(0:(n_clusters-1)))
rownames(average_profiles) = rownames(Bcells@assays$RNA@data)
for(clust in as.character(c(0:(n_clusters-1)))){
  message('Cluster:', clust)
  clust_cells = rownames(Bcells@meta.data)[Bcells@meta.data[,res] == clust]
  clust_matrix = as.matrix(Bcells@assays$RNA@data[,clust_cells])
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
# - Find markers
# -----

Seurat::Idents(Bcells) = res
Bcells_markers = Seurat::FindAllMarkers(Bcells, assay='RNA', slot='data', logfc.threshold=.8, only.pos=TRUE)
View(Bcells_markers)
write.csv(Bcells_markers, paste(project_dir, '2_annotation/markers/Bcells/Bcells_res06.csv', sep='/'))

# Get top 20 markers for each cluster:
Bcells_markers_top20 = dplyr::top_n(dplyr::group_by(Bcells_markers[Bcells_markers$p_val_adj<0.05,], cluster), n=20, wt=avg_log2FC)
View(Bcells_markers_top20)
write.csv(Bcells_markers_top20, paste(project_dir, '2_annotation/markers/Bcells/Bcells_res06_top20.csv', sep='/'))

# Merge clusters together, based on previous markers and similarity measures:
Bcells[['temp_clusters']] = as.character(Bcells@meta.data[,res])
Bcells$temp_clusters[Bcells$temp_clusters%in%c('0', '2')] = '0_2'
Bcells$temp_clusters[Bcells$temp_clusters%in%c('1', '16')] = '1_16'
Bcells$temp_clusters[Bcells$temp_clusters%in%c('4', '13')] = '4_13'
Bcells$temp_clusters[Bcells$temp_clusters%in%c('3', '5', '6', '8')] = '3_5_6_8'
Bcells$temp_clusters[Bcells$temp_clusters%in%c('9', '10')] = '9_10'
Bcells$temp_clusters[Bcells$temp_clusters%in%c('11', '17')] = '11_17'
Seurat::DimPlot(Bcells, reduction="umap", group.by=res, label=TRUE, label.size=6) |
  Seurat::DimPlot(Bcells, reduction="umap", group.by='temp_clusters', label=TRUE, label.size=6)

# Calculate markers for new clusters:
Seurat::Idents(Bcells) = 'temp_clusters'
Bcells_markers = Seurat::FindAllMarkers(Bcells, assay='RNA', slot='data', logfc.threshold=.8, only.pos=TRUE)
View(Bcells_markers)
write.csv(Bcells_markers, paste(project_dir, '2_annotation/markers/Bcells/Bcells_tempclusters.csv', sep='/'))

# Get top 20 markers for each cluster:
Bcells_markers_top20 = dplyr::top_n(dplyr::group_by(Bcells_markers[Bcells_markers$p_val_adj<0.05,], cluster), n=20, wt=avg_log2FC)
View(Bcells_markers_top20)
write.csv(Bcells_markers_top20, paste(project_dir, '2_annotation/markers/Bcells/Bcells_tempclusters_top20.csv', sep='/'))
