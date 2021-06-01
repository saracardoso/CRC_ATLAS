project_dir = '/home/scardoso/Documents/PhD/CRC_ATLAS'

# Load subset of epithelial cells:
Myeloid = SeuratDisk::LoadH5Seurat('/home/scardoso/Documents/PhD/CRC_ATLAS/2_annotation/data/myeloid_cells.h5Seurat')



# -----
# - Integrate epithelial cells
# -----

# Separate epithelial cells by dataset:
Myeloid = Seurat::SplitObject(Myeloid, split.by='dataset')
# SCTransform data:
for(dts_name in names(Myeloid)) {
  Myeloid[[dts_name]] <- Seurat::SCTransform(Myeloid[[dts_name]], vars.to.regress = c('CC.Difference', 'percent.mitochondrial_RNA'))
  gc()
}
gc()
# Select integration features and prepare for integration:
features = Seurat::SelectIntegrationFeatures(Myeloid, nfeatures = 3000)
Myeloid = Seurat::PrepSCTIntegration(Myeloid, anchor.features=features)
# Find integration achors and integrate the data:
anchors = Seurat::FindIntegrationAnchors(Myeloid, normalization.method='SCT', anchor.features=features)
gc()
Myeloid = Seurat::IntegrateData(anchorset=anchors, normalization.method="SCT")
gc()



# -----
# - Prepare data to find clusters
# -----

# Run PCA:
Myeloid = Seurat::RunPCA(Myeloid, assay='integrated')



# -----
# - Find First Clusters and Annotate Them
# -----

# 1. Find first clusters

# Determine the K-nearest neighbor graph
Myeloid = Seurat::FindNeighbors(Myeloid, dims=1:50)
# Find clusters:
Myeloid = Seurat::FindClusters(Myeloid, resolution=seq(0.1, 1, by=.1))
# UMAP:
Myeloid = Seurat::RunUMAP(Myeloid, dims=1:50, reduction='pca')
# Save object with first clusters:
SeuratDisk::SaveH5Seurat(Myeloid, paste(project_dir, '2_annotation/data/Myeloid_integratedClusters.h5Seurat', sep='/'))

# 2. Choose best resolution using clustree:
library(ggraph)
clust_tree = clustree::clustree(Myeloid, prefix='integrated_snn_res.')
clust_tree

# 3. Visualize the first two resolution (second resolution was chosen to be the best for now):
clusters_01 = Seurat::DimPlot(Myeloid, reduction="umap", group.by='integrated_snn_res.0.1', label=TRUE, label.size=6, pt.size=.8)
clusters_02 = Seurat::DimPlot(Myeloid, reduction="umap", group.by='integrated_snn_res.0.2', label=TRUE, label.size=6, pt.size=.8)
(clusters_01 | clusters_02)



# ---
# - Check SIGMA for further clusterability (res.0.2)
# ---

# 1. Get data processed for SIGMA
# 1.1. Subset cells to half of the original
n_cells = round(dim(Myeloid@assays$integrated@data)[2] / 2)
sampled_cells = colnames(Myeloid@assays$integrated@data)[sample(1:dim(Myeloid@assays$integrated@data)[2], n_cells)]
# 1.2. Check if subsetting still resembles unsubsetted dataset
clusters_02 | Seurat::DimPlot(subset(Myeloid, cells=sampled_cells), reduction="umap", group.by='integrated_snn_res.0.2',
                              label=TRUE, label.size=6, pt.size=.8)
# 1.3. Subset dataset
Myeloid_subset = subset(Myeloid, cells=sampled_cells)
# 1.4. Get raw counts data, perform normalization, log transformation and use only genes present in the integrated dataset (the one used for
#      clustering)
sigma_data = as.matrix(SeuratObject::GetAssayData(Myeloid_subset, assay='RNA', slot='counts'))
sigma_data = sigma_data[rowSums(sigma_data) > 0,]
invisible(gc())
sigma_data_norm = t(t(sigma_data)/colSums(sigma_data))*10000
invisible(gc())
sigma_data_norm = log(sigma_data_norm + 1)
invisible(gc())
sigma_data_norm = sigma_data_norm[rownames(Myeloid_subset@assays$integrated@data), ]
invisible(gc())
# 1.5. Get vector with resolution 0.2 clusters
clusters_02_vector = Myeloid_subset@meta.data[colnames(sigma_data), 'integrated_snn_res.0.2']

# 2. Factors to be excluded from considering while calculating the clusterability
data("ribosomal_genes", package='SIGMA')
data("stress_genes", package='SIGMA')
rb.genes =  intersect(rb.genes, rownames(sigma_data_norm)) 
stress.genes = intersect(stress.genes, rownames(sigma_data_norm))
exclude = data.frame(clsm = log(colSums(sigma_data[rownames(Myeloid_subset@assays$integrated@data),]) + 1),
                     cellcycle = Myeloid_subset$CC.Difference,
                     mt = colMeans(sigma_data_norm[grep("^MT-", rownames(sigma_data_norm)),]),
                     ribosomal = colMeans(sigma_data_norm[rb.genes,]), stress = colMeans(sigma_data_norm[stress.genes,]))
remove(sigma_data)
invisible(gc())

# 3. Run SIGMA
clusterability_res02_1 = SIGMA::sigma_funct(sigma_data_norm, clusters=clusters_02_vector, exclude=exclude)

# 4. Check results
# 4.1. Clusterability of all clusters
SIGMA::plot_sigma(clusterability_res02_1)
# 4.2. MP distribution for each cluster
for(i in as.character(0:9)) SIGMA::plot_MP(clusterability_res02_1, i)
# 4.3. Plot singular vectors for each cluster
for(i in as.character(0:9)){
  plot_svi = SIGMA::plot_singular_vectors(clusterability_res02_1, i) + ggplot2::ggtitle(i)
  print(plot_svi)
}
# 4.4. Analyse cluster 0
SIGMA::get_var_genes(clusterability_res02_1, '0')[,1:2]
SIGMA::plot_singular_vectors(clusterability_res02_1, '0', colour='SOD2')













# -----
# - Visualize Clusters
# -----

# Choose resolution 0.6
res = 'integrated_snn_res.0.6'

# Metadata:
state = Seurat::DimPlot(Myeloid, reduction="umap", group.by='state', label=FALSE, label.size=6)
patients = Seurat::DimPlot(Myeloid, reduction="umap", group.by='patient', label=FALSE, label.size=6, pt.size=.1)
datasets = Seurat::DimPlot(Myeloid, reduction="umap", group.by='dataset', label=FALSE, label.size=6, pt.size=.1)
(clusters_06 | patients) / (state | datasets)

# Normal/ tumour/ healthy distribution across clusters:
clusters = c()
percs = c()
fill = c()
for(cluster in unique(as.character(Myeloid@meta.data[,res]))){
  x = table(Myeloid$state[Myeloid@meta.data[,res]==cluster])[c('Tumor', 'Normal', 'Healthy')]
  x[is.na(x)] = 0
  y = sum(x)
  percs = c(percs, x/y)
  fill = c(fill, c('Tumor', 'Normal', 'Healthy'))
  clusters = c(clusters, rep(cluster, 3))
}
dist_nt = data.frame(percs=percs, fill=fill, clusters=clusters)
dist_nt$clusters = factor(dist_nt$clusters,
                          levels=as.character(0:(length(unique(Myeloid@meta.data[,res]))-1)))
dist_state = ggplot2::ggplot(dist_nt, mapping = ggplot2::aes(clusters, percs, fill=fill)) +
  ggplot2::geom_bar(position='stack', stat='identity')
dist_state

# Some gene expression:
vlngenePlot = Seurat::VlnPlot(Myeloid, features=c('rna_CD14', 'rna_FCGR3A', 'rna_MARCO', 'rna_TPSB2', 'rna_PPBP', 'rna_FCER1A',
                                                 'rna_GZMB'),
                              pt.size = 0, ncol=2, group.by = res)
clusters_06 | vlngenePlot




# -----
# - Calculate similarity between clusters
# -----

# Calculate average and median expression profiles of each cluster:
n_genes = dim(Myeloid@assays$RNA@data)[1]
n_clusters = length(unique(Myeloid@meta.data[,res]))
average_profiles = matrix(rep(0, n_genes * n_clusters), nrow = n_genes)
colnames(average_profiles) = as.character(c(0:(n_clusters-1)))
rownames(average_profiles) = rownames(Myeloid@assays$RNA@data)
for(clust in as.character(c(0:(n_clusters-1)))){
  message('Cluster:', clust)
  clust_cells = rownames(Myeloid@meta.data)[Myeloid@meta.data[,res] == clust]
  clust_matrix = as.matrix(Myeloid@assays$RNA@data[,clust_cells])
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

Seurat::Idents(Myeloid) = res
Myeloid_markers = Seurat::FindAllMarkers(Myeloid, assay='RNA', slot='data', logfc.threshold=.8, only.pos=TRUE)
View(Myeloid_markers)
write.csv(Myeloid_markers, paste(project_dir, '2_annotation/markers/Myeloid/Myeloid_res06.csv', sep='/'))

# Get top 10 markers for each cluster:
Myeloid_markers_top20 = dplyr::top_n(dplyr::group_by(Myeloid_markers[Myeloid_markers$p_val_adj<0.05,], cluster), n=20, wt=avg_log2FC)
View(Myeloid_markers_top20)
write.csv(Myeloid_markers_top20, paste(project_dir, '2_annotation/markers/Myeloid/Myeloid_res06_top20.csv', sep='/'))




# -----
# - Perform trajectory analysis for macro/mono + DCs + unknown clusters
# -----

# Convert Seurat dataset to monocle, after removing mast cells cluster (5):
monocle_Myeloid = SeuratWrappers::as.cell_data_set(subset(Myeloid, cells=rownames(Myeloid@meta.data)[Myeloid$integrated_snn_res.0.6!='5']),
                                                   assay='integrated')

# Pre-process and reduce dimensions using monocle??

# Cluster cells to get different partitions (cells belonging to different trajectories, i.e., different ancestors):
monocle_Myeloid = monocle3::cluster_cells(monocle_Myeloid)
p1 = monocle3::plot_cells(monocle_Myeloid, color_cells_by="integrated_snn_res.0.6", show_trajectory_graph=FALSE, label_cell_groups=FALSE, 
                          labels_per_group=1, cell_size=.8, group_label_size=5) + ggplot2::ggtitle("Colored by Seurat's clusters")
p2 = monocle3::plot_cells(monocle_Myeloid, color_cells_by="cluster", show_trajectory_graph=FALSE, label_cell_groups=FALSE, 
                          labels_per_group=1, cell_size=.8, group_label_size=5) + ggplot2::ggtitle("Colored by Monocle's clusters")
p3 = monocle3::plot_cells(monocle_Myeloid, color_cells_by="partition", show_trajectory_graph=FALSE, label_cell_groups=FALSE, 
                          labels_per_group=1, cell_size=.8, group_label_size=5) + ggplot2::ggtitle("Monocle's partitions")
p1 | p2 | p3

# Learn the trajectory graph:
monocle_Myeloid = monocle3::learn_graph(monocle_Myeloid)
tj1 = monocle3::plot_cells(monocle_Myeloid, color_cells_by = "integrated_snn_res.0.6", label_cell_groups=FALSE, labels_per_group=2,
                           label_leaves=TRUE, label_branch_points=TRUE, graph_label_size=2, group_label_size=6, cell_size=.8) +
  ggplot2::ggtitle("Trajectory with partitions") + Seurat::NoLegend()
monocle_Myeloid = monocle3::learn_graph(monocle_Myeloid, use_partition=FALSE)
tj2 = monocle3::plot_cells(monocle_Myeloid, color_cells_by = "integrated_snn_res.0.6", label_cell_groups=FALSE, labels_per_group=2,
                           label_leaves=TRUE, label_branch_points=TRUE, graph_label_size=2, group_label_size=6, cell_size=.8) +
  ggplot2::ggtitle("Trajectory without partitions") + Seurat::NoLegend()

seurat_clusters = Seurat::DimPlot(subset(Myeloid, cells=rownames(Myeloid@meta.data)[Myeloid$integrated_snn_res.0.6!='5']),
                                  reduction="umap", group.by='integrated_snn_res.0.6', label=TRUE, label.size=6, pt.size=.8)
(seurat_clusters / p2) | tj1 | tj2

