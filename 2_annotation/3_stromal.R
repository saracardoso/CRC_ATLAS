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
pct = Stromal[["pca"]]@stdev /sum(Stromal[["pca"]]@stdev) * 100
cumu = cumsum(pct)
co1 = which(cumu > 90 & pct < 5)[1]
co2 = sort(which((pct[1:length(pct) - 1] - pct[2:length(pct)]) > 0.1), decreasing = T)[1] + 1
elbow = min(co1, co2) # 20

plot_df = data.frame(pct=pct, cumu=cumu, rank=1:length(pct))
ggplot2::ggplot(plot_df, ggplot2::aes(cumu, pct, label = rank, color = rank > elbow)) + 
  ggplot2::geom_text() + 
  ggplot2::geom_vline(xintercept = 90, color = "grey") + 
  ggplot2::geom_hline(yintercept = min(pct[pct > 5]), color = "grey") +
  ggplot2::theme_bw()

Seurat::DimHeatmap(Stromal, dims=1:21, cells=500, balanced=TRUE)



# -----
# - Find Clusters
# -----

# Determine the K-nearest neighbor graph
Stromal = Seurat::FindNeighbors(Stromal, dims=1:elbow)

# Find clusters:
Stromal = Seurat::FindClusters(Stromal, resolution=c(0.4, 0.6, 1))

# UMAP:
Stromal = Seurat::RunUMAP(Stromal, dims=1:elbow, reduction='pca')

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

# Some gene expression:
vlngenePlot = Seurat::VlnPlot(Stromal, features=c('rna_COL1A1', 'rna_VWF', 'rna_SOX10', 'rna_SYNPO2', 'rna_RGS5', 'rna_TAGLN', 'rna_THY1',
                                                  'rna_FAP'),
                              pt.size = 0, ncol=2, group.by = res)
clusters_06 | vlngenePlot



# -----
# - Find markers
# -----

Seurat::Idents(Stromal) = res
Stromal_markers = Seurat::FindAllMarkers(Stromal, assay='RNA', slot='data', logfc.threshold=1, only.pos=TRUE)
View(Stromal_markers)
write.csv(Stromal_markers, paste(project_dir, '2_annotation/markers/Stromal/Stromal_res06.csv', sep='/'))

# Get top 10 markers for each cluster:
Stromal_markers_top20 = dplyr::top_n(dplyr::group_by(Stromal_markers[Stromal_markers$p_val_adj<0.05,], cluster), n=20, wt=avg_log2FC)
View(Stromal_markers_top20)
write.csv(Stromal_markers_top20, paste(project_dir, '2_annotation/markers/Stromal/Stromal_res06_top20.csv', sep='/'))

# Clusters highlighted in the full atlas:
CRCatlas_integrated = SeuratDisk::LoadH5Seurat('/home/scardoso/Documents/PhD/CRC_ATLAS/2_annotation/data/CRC_annotations.h5Seurat')
for(cluster in 0:(length(unique(Stromal_markers_top20$cluster))-1)){
  cluster_cells = rownames(Stromal@meta.data)[Stromal@meta.data$integrated_snn_res.0.6==cluster]
  grDevices::png(filename = paste(project_dir, '/2_annotation/figures/Stromal/Cluster',
                                  cluster, '.png', sep=''), width=805, height=461)
  print(Seurat::DimPlot(CRCatlas_integrated, cells.highlight = cluster_cells) +
          ggplot2::ggtitle(paste('Cluster', cluster)) + Seurat::NoLegend())
  dev.off()
}

# Heatmap with the gene markers:
sub_cells = subset(Stromal, downsample=600)
sub_cells = Seurat::NormalizeData(sub_cells, assay='RNA')
sub_cells = Seurat::ScaleData(sub_cells, assay='RNA')
for(cluster_i in unique(Stromal_markers_top20$cluster)){
  gene_markers = Stromal_markers_top20$gene[Stromal_markers_top20$cluster==cluster_i]
  grDevices::png(filename = paste(project_dir, '/2_annotation/figures/Stromal/Heatmap',
                                  cluster_i, '.png', sep=''), width=1500, height=800)
  print(Seurat::DoHeatmap(sub_cells, features=gene_markers, assay='RNA') + ggplot2::ggtitle(paste('Cluster', cluster_i, 'markers')))
  dev.off()
}

Seurat::DoHeatmap(sub_cells, features=c('COL1A1', 'COL1A2', 'COL6A1', 'COL6A2', 'SPARC', 'COL14A', 'COL3A1', 'DCN',
                                        'VWF', 'CDH5', 'PLVAP', 'ENG',
                                        'SOX10', 'S100B', 'GFAP',
                                        'SYNPO2', 'CRYAB', 'CNN1', 'DES',
                                        'RGS5', 'RGS5', 'CSPG4', 'ABCC9', 'KCNJ8',
                                        'TAGLN', 'ACTA2', 'ACTG2', 'MYH11', 'MYLK', 'MMP1', 'MMP3',
                                        'THY1', 'FAP'), assay='RNA')

