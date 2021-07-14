project_dir = '/home/scardoso/Documents/PhD/CRC_ATLAS'

# Load integrated CRC atlas:
CRCatlas_integrated = SeuratDisk::LoadH5Seurat(paste(project_dir, '1_merge/CRC_integrated.h5Seurat', sep='/'))



# -----
# - Prepare data to find clusters
# -----

# Run PCA:
CRCatlas_integrated = Seurat::ScaleData(CRCatlas_integrated, assay='integrated')
CRCatlas_integrated = Seurat::RunPCA(CRCatlas_integrated, assay='integrated')
gc()

# Choose number of PCs to use (where the elbow occurs):
pct = CRCatlas_integrated[["pca"]]@stdev /sum(CRCatlas_integrated[["pca"]]@stdev) * 100
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

Seurat::DimHeatmap(CRCatlas_integrated, dims=1:21, cells=500, balanced=TRUE)



# -----
# - Find Clusters
# -----

# Determine the K-nearest neighbor graph
CRCatlas_integrated = Seurat::FindNeighbors(CRCatlas_integrated, dims=1:elbow)

# Find clusters (the resolution is set to a small value because we first only want to separate the cells into bigger clusters):
CRCatlas_integrated = Seurat::FindClusters(CRCatlas_integrated, resolution = 0.1)

# Run UMAP visualization:
CRCatlas_integrated = Seurat::RunUMAP(CRCatlas_integrated, dims=1:elbow, reduction='pca')



# -----
# - Evaluate quality of clusters
# -----

# Color clusters according to metadata:
Seurat::DimPlot(CRCatlas_integrated, reduction="umap", group.by='dataset', label=FALSE)
Seurat::DimPlot(CRCatlas_integrated, reduction="umap", group.by='integrated_snn_res.0.1', label=TRUE, label.size=5)
Seurat::DimPlot(CRCatlas_integrated, reduction="umap", group.by='Phase', label=FALSE)
Seurat::DimPlot(CRCatlas_integrated, reduction="umap", group.by='integrated_snn_res.0.1', split.by='Phase', label=FALSE)
Seurat::DimPlot(CRCatlas_integrated, reduction="umap", group.by='state', label=FALSE)
Seurat::DimPlot(CRCatlas_integrated, reduction="umap", group.by='integrated_snn_res.0.1', split.by='state', label=T)

# Explore if certain metrics are source of variation between clusters:
metrics =  c("nUMI", "nGene", "S.Score", "G2M.Score", "", "percent.mitochondrial_RNA")
Seurat::FeaturePlot(CRCatlas_integrated, reduction = "umap", features = metrics,
                    pt.size = 0.4, order = TRUE, min.cutoff = 'q10', label = TRUE)

# Normal/ tumour/ healthy distribution across clusters:
clusters = c()
percs = c()
fill = c()
for(cluster in unique(as.character(CRCatlas_integrated@meta.data[,'integrated_snn_res.0.1']))){
  x = table(CRCatlas_integrated$state[CRCatlas_integrated@meta.data[,'integrated_snn_res.0.1']==cluster])[c('Tumor', 'Normal', 'Healthy')]
  x[is.na(x)] = 0
  y = sum(x)
  percs = c(percs, x/y)
  fill = c(fill, c('Tumor', 'Normal', 'Healthy'))
  clusters = c(clusters, rep(cluster, 3))
}
dist_nt = data.frame(percs=percs, fill=fill, clusters=clusters)
dist_nt$clusters = factor(dist_nt$clusters,
                          levels=as.character(0:(length(unique(CRCatlas_integrated@meta.data[,'integrated_snn_res.0.1']))-1)))
dist_state = ggplot2::ggplot(dist_nt, mapping = ggplot2::aes(clusters, percs, fill=fill)) +
  ggplot2::geom_bar(position='stack', stat='identity')
dist_state




# -----
# - Visualize Clusters
# -----

dimplot = Seurat::DimPlot(CRCatlas_integrated, reduction="umap", group.by='integrated_snn_res.0.1', label=TRUE, label.size=6)
dimplot

dimplot_normTum = Seurat::DimPlot(CRCatlas_integrated, reduction="umap", group.by='state', label=FALSE, label.size=6)
dimplot_normTum

# T cell markers:
tcells = Seurat::FeaturePlot(CRCatlas_integrated, features='rna_CD3D')
dimplot | tcells
# B cell markers:
bcells = Seurat::FeaturePlot(CRCatlas_integrated, features=c('rna_CD79A', 'rna_MZB1'), min.cutoff='q10', ncol=1)
dimplot | bcells
# Monocytes markers
monocytes = Seurat::FeaturePlot(CRCatlas_integrated, features=c('rna_CD14', 'rna_FCGR3A'), min.cutoff='q10', ncol=1)
dimplot | monocytes
# DCs markers
dcs = Seurat::FeaturePlot(CRCatlas_integrated, features=c('rna_FCER1A', 'rna_GZMB'), min.cutoff='q10', ncol=1)
dimplot | dcs
# Mast cells markers
mast = Seurat::FeaturePlot(CRCatlas_integrated, features=c('rna_TPSAB1', 'rna_TPSB2'), min.cutoff='q10', ncol=1)
dimplot | mast
# Macrophages markers
macro = Seurat::FeaturePlot(CRCatlas_integrated, features=c('rna_MARCO', 'rna_ITGAM'), ncol=1)
dimplot | macro
# Megakaryocytes/ platelets
mega = Seurat::FeaturePlot(CRCatlas_integrated, features='rna_PPBP', min.cutoff='q10')
dimplot | mega
# Erythrocytes
ery = Seurat::FeaturePlot(CRCatlas_integrated, features=c('rna_HBB', 'rna_HBA2'), min.cutoff='q10', ncol=1)
dimplot | ery
# Fibroblasts
fibro = Seurat::FeaturePlot(CRCatlas_integrated, features='rna_COL1A1', min.cutoff='q10')
dimplot | fibro
# CAFs
cafs = Seurat::FeaturePlot(CRCatlas_integrated, features=c('rna_THY1', 'rna_FAP'), min.cutoff='q10', ncol=1)
dimplot | cafs
# Endothelial cells
endo = Seurat::FeaturePlot(CRCatlas_integrated, features=c('rna_VWF', 'rna_ENG'), min.cutoff='q10', ncol=1)
dimplot | endo
# Enteric glia
enteric = Seurat::FeaturePlot(CRCatlas_integrated, features='rna_S100B', min.cutoff='q10')
dimplot | enteric
# Epithelial cells
epithelial = Seurat::FeaturePlot(CRCatlas_integrated, features='rna_EPCAM', min.cutoff='q10')
dimplot | epithelial
dimplot_normTum | epithelial
# Oncogenes
onco = Seurat::FeaturePlot(CRCatlas_integrated, features=c('rna_MYC', 'rna_AXIN2', 'rna_RNF43', 'rna_AREG'), min.cutoff='q10', ncol=2)
dimplot | onco



# -----
# - Annotate Clusters
# -----

"
(Resolution 0.1)
Tcells - 0
Bcells - 4, 8
Myeloid cells - 2, 12
Stromal cells - 3, 6, 7, 10, 11
Epithelial cells - 1, 5, 9, 13
"

CRCatlas_integrated[['Level_0']] = as.character(CRCatlas_integrated$integrated_snn_res.0.1)
CRCatlas_integrated$Level_0[CRCatlas_integrated$Level_0=='0'] = 'Tcells'
CRCatlas_integrated$Level_0[CRCatlas_integrated$Level_0%in%c('4', '8')] = 'Bcells'
CRCatlas_integrated$Level_0[CRCatlas_integrated$Level_0%in%c('2', '12')] = 'Myeloid cells'
CRCatlas_integrated$Level_0[CRCatlas_integrated$Level_0%in%c('3', '6', '7', '10', '11')] = 'Stromal cells'
CRCatlas_integrated$Level_0[CRCatlas_integrated$Level_0%in%c('1', '5', '9', '13')] = 'Epithelial cells'

Seurat::DimPlot(CRCatlas_integrated, reduction="umap", group.by='Level_0', label=FALSE, label.size=6)
Seurat::DimPlot(CRCatlas_integrated, reduction="umap", split.by='Level_0', group.by='Level_0', label=FALSE, pt.size=1)

# Store integrated CRCatlas for annotation purposes:
SeuratDisk::SaveH5Seurat(CRCatlas_integrated, paste(project_dir, '2_annotation/results_globalAnnotation/datasets/CRC_annotations.h5Seurat', sep='/'))



# -----
# - Subset dataset for further annotation
# -----

Seurat::DefaultAssay(CRCatlas_integrated) = 'RNA'
CRCatlas_integrated@meta.data = CRCatlas_integrated@meta.data[, !colnames(CRCatlas_integrated@meta.data) %in%
                                                                c('integrated_snn_res.0.1', 'seurat_clusters')]

# Epithelial cells: (these will not be re-integrated)
epithelial_cells = subset(CRCatlas_integrated,
                          cells=rownames(CRCatlas_integrated@meta.data)[CRCatlas_integrated$Level_0=='Epithelial cells'])
SeuratDisk::SaveH5Seurat(epithelial_cells, paste(project_dir, '/2_annotation/results_Epithelial/datasets/epithelial_cells.h5Seurat', sep='/'))

# Stromal cells:
stromal_cells = subset(CRCatlas_integrated,
                       cells=rownames(CRCatlas_integrated@meta.data)[CRCatlas_integrated$Level_0=='Stromal cells'])
SeuratDisk::SaveH5Seurat(stromal_cells, paste(project_dir, '/2_annotation/results_Stromal/datasets/stromal_cells.h5Seurat', sep='/'))

# Myeloid cells:
myeloid_cells = subset(CRCatlas_integrated,
                       cells=rownames(CRCatlas_integrated@meta.data)[CRCatlas_integrated$Level_0=='Myeloid cells'])
SeuratDisk::SaveH5Seurat(myeloid_cells, paste(project_dir, '/2_annotation/results_Myeloid/datasets/myeloid_cells.h5Seurat', sep='/'))

# Bcells:
Bcells = subset(CRCatlas_integrated,
                cells=rownames(CRCatlas_integrated@meta.data)[CRCatlas_integrated$Level_0=='Bcells'])
SeuratDisk::SaveH5Seurat(Bcells, paste(project_dir, '/2_annotation/results_Bcells/datasets/Bcells.h5Seurat', sep='/'))

# Tcells:
Tcells = subset(CRCatlas_integrated, cells=rownames(CRCatlas_integrated@meta.data)[CRCatlas_integrated$Level_0=='Tcells'])
SeuratDisk::SaveH5Seurat(Tcells, '/home/scardoso/Documents/PhD/CRC_ATLAS/2_annotation/data/Tcells.h5Seurat')

