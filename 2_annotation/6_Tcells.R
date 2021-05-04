project_dir = '/home/scardoso/Documents/PhD/CRC_ATLAS'

# Load subset of epithelial cells:
Tcells = SeuratDisk::LoadH5Seurat('/home/scardoso/Documents/PhD/CRC_ATLAS/2_annotation/data/Tcells.h5Seurat')



# -----
# - Integrate epithelial cells
# -----

# Separate epithelial cells by dataset:
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



# -----
# - Prepare data to find clusters
# -----

# Run PCA:
Tcells = Seurat::RunPCA(Tcells, assay='integrated')
# Find the number of PCs to use:
pct = Tcells[["pca"]]@stdev /sum(Tcells[["pca"]]@stdev) * 100
cumu = cumsum(pct)
co1 = which(cumu > 90 & pct < 5)[1]
co2 = sort(which((pct[1:length(pct) - 1] - pct[2:length(pct)]) > 0.1), decreasing = T)[1] + 1
elbow = min(co1, co2) # 15

plot_df = data.frame(pct=pct, cumu=cumu, rank=1:length(pct))
ggplot2::ggplot(plot_df, ggplot2::aes(cumu, pct, label = rank, color = rank > elbow)) + 
  ggplot2::geom_text() + 
  ggplot2::geom_vline(xintercept = 90, color = "grey") + 
  ggplot2::geom_hline(yintercept = min(pct[pct > 5]), color = "grey") +
  ggplot2::theme_bw()

Seurat::DimHeatmap(Tcells, dims=1:21, cells=500, balanced=TRUE)



# -----
# - Find Clusters
# -----

# Determine the K-nearest neighbor graph
Tcells = Seurat::FindNeighbors(Tcells, dims=1:elbow)

# Find clusters:
Tcells = Seurat::FindClusters(Tcells, resolution=c(0.4, 0.6, 1))

# UMAP:
Tcells = Seurat::RunUMAP(Tcells, dims=1:elbow, reduction='pca')

# Save object with clusters:
SeuratDisk::SaveH5Seurat(Tcells, paste(project_dir, '2_annotation/data/Tcells_integratedClusters.h5Seurat', sep='/'))



# -----
# - Visualize Clusters
# -----

# Cluster from the three resolutions:
clusters_04 = Seurat::DimPlot(Tcells, reduction="umap", group.by='integrated_snn_res.0.4', label=TRUE, label.size=6)
clusters_06 = Seurat::DimPlot(Tcells, reduction="umap", group.by='integrated_snn_res.0.6', label=TRUE, label.size=6)
clusters_1 = Seurat::DimPlot(Tcells, reduction="umap", group.by='integrated_snn_res.1', label=TRUE, label.size=6)
(clusters_04 | clusters_06) / (clusters_1 | ggplot2::ggplot())

# Choose resolution 0.6
res = 'integrated_snn_res.0.6'

# Metadata:
state = Seurat::DimPlot(Tcells, reduction="umap", group.by='state', label=FALSE, label.size=6)
patients = Seurat::DimPlot(Tcells, reduction="umap", group.by='patient', label=FALSE, label.size=6, pt.size=.1)
datasets = Seurat::DimPlot(Tcells, reduction="umap", group.by='dataset', label=FALSE, label.size=6, pt.size=.1)
(clusters_06 | patients) / (state | datasets)

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
Tcells_markers = Seurat::FindAllMarkers(Tcells, assay='RNA', slot='data', logfc.threshold=1, only.pos=TRUE)
View(Tcells_markers)
write.csv(Tcells_markers, paste(project_dir, '2_annotation/markers/Tcells/Tcells_res06.csv', sep='/'))

# Get top 10 markers for each cluster:
Tcells_markers_top20 = dplyr::top_n(dplyr::group_by(Tcells_markers[Tcells_markers$p_val_adj<0.05,], cluster), n=20, wt=avg_log2FC)
View(Tcells_markers_top20)
write.csv(Tcells_markers_top20, paste(project_dir, '2_annotation/markers/Tcells/Tcells_res06_top20.csv', sep='/'))

# Clusters highlighted in the full atlas:
CRCatlas_integrated = SeuratDisk::LoadH5Seurat('/home/scardoso/Documents/PhD/CRC_ATLAS/2_annotation/data/CRC_annotations.h5Seurat')
for(cluster in as.character(0:15)){
  cluster_cells = rownames(Tcells@meta.data)[Tcells@meta.data$integrated_snn_res.0.6==cluster]
  grDevices::png(filename = paste(project_dir, '/2_annotation/figures/Tcells/Cluster',
                                  cluster, '.png', sep=''), width=805, height=461)
  print(Seurat::DimPlot(CRCatlas_integrated, cells.highlight = cluster_cells) +
          ggplot2::ggtitle(paste('Cluster', cluster)) + Seurat::NoLegend())
  dev.off()
}

# Heatmap with the gene markers:
Tcells = Seurat::NormalizeData(Tcells, assay='RNA')
Tcells = Seurat::ScaleData(Tcells, assay='RNA')
gene_markers = unique(Tcells_markers_top20$gene)
sub_cells = subset(Tcells, downsample=600)
Seurat::DoHeatmap(sub_cells, features=gene_markers[1:60], assay='RNA', group.by = res)
Seurat::DoHeatmap(sub_cells, features=gene_markers[61:120], assay='RNA', group.by = res)
Seurat::DoHeatmap(sub_cells, features=gene_markers[121:184], assay='RNA', group.by = res)

Seurat::DoHeatmap(sub_cells, features=gene_markers[1:60], assay='RNA', group.by = 'state')
Seurat::DoHeatmap(sub_cells, features=gene_markers[61:120], assay='RNA', group.by = 'state')
Seurat::DoHeatmap(sub_cells, features=gene_markers[121:184], assay='RNA', group.by = 'state')





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
