project_dir = '/home/scardoso/Documents/PhD/CRC_ATLAS'
source(paste(project_dir, 'utils/modified_plots.R', sep='/'))

# Load subset of epithelial cells:
Stromal = SeuratDisk::LoadH5Seurat(paste(project_dir, '2_annotation/results_Stromal/datasets/stromal_cells.h5Seurat', sep='/'))





# -----
# - Find clusters
# -----


# 1. Run PCA:
Stromal = Seurat::RunPCA(Stromal, assay='integrated')
invisible(gc())
Seurat::DefaultAssay(Stromal)='integrated'

# 2. Choose number of PCs to use (where the elbow occurs):
pct = Stromal[["pca"]]@stdev /sum(Stromal[["pca"]]@stdev) * 100
cumu = cumsum(pct)
co1 = which(cumu > 90 & pct < 5)[1]
co2 = sort(which((pct[1:length(pct) - 1] - pct[2:length(pct)]) > 0.1), decreasing = T)[1] + 1
# 2.1 Number of PCs to use
elbow = min(co1, co2) # 14
# 2.2. Visual representation of the PCs to use
plot_df = data.frame(pct=pct, cumu=cumu, rank=1:length(pct))
ggplot2::ggplot(plot_df, ggplot2::aes(cumu, pct, label = rank, color = rank > elbow)) + 
  ggplot2::geom_text() + 
  ggplot2::geom_vline(xintercept = 90, color = "grey") + 
  ggplot2::geom_hline(yintercept = min(pct[pct > 5]), color = "grey") +
  ggplot2::theme_bw()
# 2.3. Heatmap of the first 15 PCs
Seurat::DimHeatmap(Stromal, dims=1:15, cells=500, balanced=TRUE)


# 3. Find clusters
# 3.1. Determine the K-nearest neighbor graph
Stromal = Seurat::FindNeighbors(Stromal, dims=1:elbow)
# 3.2. Find clusters:
Stromal = Seurat::FindClusters(Stromal, resolution=seq(0.1, 1, by=.1))
# 3.3. UMAP:
Stromal = Seurat::RunUMAP(Stromal, dims=1:elbow, reduction='pca')
invisible(gc())


# 4. Choose best resolution:
library(ggraph)
clust_tree = clustree::clustree(Stromal, prefix='integrated_snn_res.')
clust_tree
# 4.1. Visualize different resolutions (first resolution was chosen to be the best for now):
clusters_01 = Seurat::DimPlot(Stromal, reduction="umap", group.by='integrated_snn_res.0.1', label=TRUE, label.size=6, pt.size=.2)
clusters_03 = Seurat::DimPlot(Stromal, reduction="umap", group.by='integrated_snn_res.0.3', label=TRUE, label.size=6, pt.size=.2)
clusters_05 = Seurat::DimPlot(Stromal, reduction="umap", group.by='integrated_snn_res.0.5', label=TRUE, label.size=6, pt.size=.2)
clusters_07 = Seurat::DimPlot(Stromal, reduction="umap", group.by='integrated_snn_res.0.7', label=TRUE, label.size=6, pt.size=.2)
((clusters_01 | clusters_03) / (clusters_05 | clusters_07)) | clust_tree
# 4.2. Plot Metrics
metrics =  c("nUMI", "nGene", "S.Score", "G2M.Score", "percent.mitochondrial_RNA")
feature_plots(Stromal, metrics, ncol=3, with_dimplot=TRUE, dimplot_group='integrated_snn_res.0.5')
# 4.3. Plot metadata
state = Seurat::DimPlot(Stromal, reduction="umap", group.by='state', label=FALSE, label.size=6) + Seurat::NoLegend()
patients = Seurat::DimPlot(Stromal, reduction="umap", group.by='patient', label=FALSE, label.size=6, pt.size=.1) + Seurat::NoLegend()
datasets = Seurat::DimPlot(Stromal, reduction="umap", group.by='dataset', label=FALSE, label.size=6, pt.size=.1) + Seurat::NoLegend()
(clusters_05 | patients) / (state | datasets)
# 4.4. Normal/ tumour/ healthy distribution across clusters:
clusters = c()
percs = c()
fill = c()
for(cluster in unique(as.character(Stromal@meta.data[,'integrated_snn_res.0.5']))){
  x = table(Stromal$state[Stromal@meta.data[,'integrated_snn_res.0.5']==cluster])[c('Tumor', 'Normal', 'Healthy')]
  x[is.na(x)] = 0
  y = sum(x)
  percs = c(percs, x/y)
  fill = c(fill, c('Tumor', 'Normal', 'Healthy'))
  clusters = c(clusters, rep(cluster, 3))
}
dist_nt = data.frame(percs=percs, fill=fill, clusters=clusters)
dist_nt$clusters = factor(dist_nt$clusters,
                          levels=as.character(0:(length(unique(Stromal@meta.data[,'integrated_snn_res.0.5']))-1)))
dist_state = ggplot2::ggplot(dist_nt, mapping = ggplot2::aes(clusters, percs, fill=fill)) +
  ggplot2::geom_bar(position='stack', stat='identity')
dist_state
# 4.5. Feature Plots:
feature_plots(Stromal, c('rna_VWF', 'rna_PLVAP', 'rna_CDH5', 'rna_ENG', 'rna_LYVE1', 'rna_PROX1'),
              ncol=3, with_dimplot=TRUE, dimplot_group='integrated_snn_res.0.5') # Vascular and lymphatic endothelial cells
feature_plots(Stromal, c('rna_PLP1', 'rna_S100B'),
              ncol=3, with_dimplot=TRUE, dimplot_group='integrated_snn_res.0.5') # Enteric glia cells
feature_plots(Stromal, c('rna_SYNPO2', 'rna_CNN1', 'rna_PDGFRB', 'rna_RGS5', 'rna_ABCC9', 'rna_KCNJ8'),
              ncol=3, with_dimplot=TRUE, dimplot_group='integrated_snn_res.0.5') # VSMCs + pericytes
feature_plots(Stromal, c('rna_COL1A1', 'rna_COL1A2', 'rna_COL6A1', 'rna_COL6A2', 'rna_COL3A1', 'rna_DCN'),
              ncol=3, with_dimplot=TRUE, dimplot_group='integrated_snn_res.0.5') # fibroblasts
# 4.6. Resolution 0.5 is chosen

# 5. Save object:
SeuratDisk::SaveH5Seurat(Stromal, paste(project_dir, '2_annotation/results_Stromal/datasets/Stromal_clustered.h5Seurat', sep='/'))





# -----
# - Annotate Clusters
# -----

# 1. Assess markers
# 1.1. Vascular (stalk- and tip- like) and lymphatic endothelial cells
feature_plots(Stromal, c('rna_VWF', 'rna_PLVAP', 'rna_CDH5', 'rna_ENG',
                         'rna_LYVE1', 'rna_PROX1',
                         'rna_ACKR1', 'rna_SELP',
                         'rna_RGCC', 'rna_RAMP3'),
              ncol=3, with_dimplot=TRUE, dimplot_group='integrated_snn_res.0.5')
violin_plots(Stromal, c('rna_VWF', 'rna_PLVAP', 'rna_CDH5', 'rna_ENG',
                         'rna_LYVE1', 'rna_PROX1',
                         'rna_ACKR1', 'rna_SELP',
                         'rna_RGCC', 'rna_RAMP3'),
              ncol=4, with_dimplot=TRUE, group.by='integrated_snn_res.0.5')

# 1.2. Enteric glia cells
feature_plots(Stromal, c('rna_PLP1', 'rna_S100B'), ncol=3, with_dimplot=TRUE, dimplot_group='integrated_snn_res.0.5')
violin_plots(Stromal, c('rna_PLP1', 'rna_S100B'), ncol=3, with_dimplot=TRUE, group.by='integrated_snn_res.0.5')

# 1.3. VSMCs + Pericytes
feature_plots(Stromal, c('rna_SYNPO2', 'rna_CNN1', 'rna_PDGFRB',
                         'rna_RGS5', 'rna_ABCC9', 'rna_KCNJ8'),
              ncol=3, with_dimplot=TRUE, dimplot_group='integrated_snn_res.0.5')
violin_plots(Stromal, c('rna_SYNPO2', 'rna_CNN1', 'rna_PDGFRB',
                         'rna_RGS5', 'rna_ABCC9', 'rna_KCNJ8'),
              ncol=4, with_dimplot=TRUE, group.by='integrated_snn_res.0.5')

# 1.4. Myofibroblasts
feature_plots(Stromal, c('rna_TAGLN', 'rna_ACTA2', 'rna_ACTG2', 'rna_MYH11', 'rna_MYLK', 'rna_DES'),
              ncol=3, with_dimplot=TRUE, dimplot_group='integrated_snn_res.0.5')
violin_plots(Stromal, c('rna_TAGLN', 'rna_ACTA2', 'rna_ACTG2', 'rna_MYH11', 'rna_MYLK', 'rna_DES'),
             ncol=4, with_dimplot=TRUE, group.by='integrated_snn_res.0.5')

# 1.5. Fibroblasts + CAFs
feature_plots(Stromal, c('rna_COL1A1', 'rna_COL1A2', 'rna_COL6A1', 'rna_COL6A2', 'rna_COL3A1', 'rna_DCN',
                         'rna_THY1', 'rna_FAP'),
              ncol=3, with_dimplot=TRUE, dimplot_group='integrated_snn_res.0.5')
violin_plots(Stromal, c('rna_COL1A1', 'rna_COL1A2', 'rna_COL6A1', 'rna_COL6A2', 'rna_COL3A1', 'rna_DCN',
                         'rna_THY1', 'rna_FAP'),
              ncol=3, with_dimplot=TRUE, group.by='integrated_snn_res.0.5')


# 2. Annotate cells:
# Sub-cluster  0           --> Fibroblasts_1
# Sub-cluster  1 + 10      --> CAFs_1
# Sub-cluster  2           --> Fibroblasts_2
# Sub-cluster  3           --> Fibroblasts_3
# Sub-cluster  4           --> Fibroblasts_4
# Sub-cluster  5           --> Pericytes
# Sub-clusters 6 + 8 + 15  --> Tip-like vascular ECs
# Sub-cluster  7           --> Stalk-like vascular ECs
# Sub-cluster  9           --> Enteric glia cells
# Sub-cluster  11          --> VSMCs
# Sub-cluster  12          --> CAFs_2
# Sub-cluster  13          --> Myofibroblasts
# Sub-cluster  14          --> Lymphatic ECs

# 2.1. Annotation level 3:
Stromal[['Annotation_Level_3']] = rep('', length(Stromal$integrated_snn_res.0.5))
fib1_cells = rownames(Stromal@meta.data)[Stromal$integrated_snn_res.0.5=='0']
Stromal@meta.data[fib1_cells, 'Annotation_Level_3'] = 'Fibroblasts_1'
cafs1_cells = rownames(Stromal@meta.data)[Stromal$integrated_snn_res.0.5%in%c('1', '10')]
Stromal@meta.data[cafs1_cells, 'Annotation_Level_3'] = 'CAFs_1'
fib2_cells = rownames(Stromal@meta.data)[Stromal$integrated_snn_res.0.5=='2']
Stromal@meta.data[fib2_cells, 'Annotation_Level_3'] = 'Fibroblasts_2'
fib3_cells = rownames(Stromal@meta.data)[Stromal$integrated_snn_res.0.5=='3']
Stromal@meta.data[fib3_cells, 'Annotation_Level_3'] = 'Fibroblasts_3'
fib4_cells = rownames(Stromal@meta.data)[Stromal$integrated_snn_res.0.5=='4']
Stromal@meta.data[fib4_cells, 'Annotation_Level_3'] = 'Fibroblasts_4'
pericytes_cells = rownames(Stromal@meta.data)[Stromal$integrated_snn_res.0.5=='5']
Stromal@meta.data[pericytes_cells, 'Annotation_Level_3'] = 'Pericytes'
tip_cells = rownames(Stromal@meta.data)[Stromal$integrated_snn_res.0.5%in%c('6', '8', '15')]
Stromal@meta.data[tip_cells, 'Annotation_Level_3'] = 'Tip-like vascular ECs'
stalk_cells = rownames(Stromal@meta.data)[Stromal$integrated_snn_res.0.5=='7']
Stromal@meta.data[stalk_cells, 'Annotation_Level_3'] = 'Stalk-like vascular ECs'
enteric_cells = rownames(Stromal@meta.data)[Stromal$integrated_snn_res.0.5=='9']
Stromal@meta.data[enteric_cells, 'Annotation_Level_3'] = 'Enteric glia cells'
vsmcs_cells = rownames(Stromal@meta.data)[Stromal$integrated_snn_res.0.5=='11']
Stromal@meta.data[vsmcs_cells, 'Annotation_Level_3'] = 'VSMCs'
cafs2_cells = rownames(Stromal@meta.data)[Stromal$integrated_snn_res.0.5=='12']
Stromal@meta.data[cafs2_cells, 'Annotation_Level_3'] = 'CAFs_2'
myofibroblasts_cells = rownames(Stromal@meta.data)[Stromal$integrated_snn_res.0.5=='13']
Stromal@meta.data[myofibroblasts_cells, 'Annotation_Level_3'] = 'Myofibroblasts'
lymphatic_cells = rownames(Stromal@meta.data)[Stromal$integrated_snn_res.0.5=='14']
Stromal@meta.data[lymphatic_cells, 'Annotation_Level_3'] = 'Lymphatic ECs'

# 2.2. Annotation level 1:
Stromal[['Annotation_Level_2']] = as.character(Stromal$Annotation_Level_3)
Stromal@meta.data[Stromal$Annotation_Level_3%in%c('Fibroblasts_1', 'Fibroblasts_2', 'Fibroblasts_3', 'Fibroblasts_4'),
                  'Annotation_Level_2'] = 'Fibroblasts'
Stromal@meta.data[Stromal$Annotation_Level_3%in%c('CAFs_1', 'CAFs_2'), 'Annotation_Level_2'] = 'CAFs'

# 2.3. Annotation level 1:
Stromal[['Annotation_Level_1']] = as.character(Stromal$Annotation_Level_2)
Stromal@meta.data[Stromal$Annotation_Level_2%in%c('Tip-like vascular ECs', 'Stalk-like vascular ECs'), 'Annotation_Level_1'] = 'Vascular ECs'

# 2.4. Annotation level 0:
Stromal[['Annotation_Level_0']] = factor(rep('Stromal cells', dim(Stromal@meta.data)[1]))

# 2.5. Remove [integrated...] metadata variables
Stromal@meta.data = Stromal@meta.data[, !colnames(Stromal@meta.data) %in%
                                      c(grep('^in', colnames(Stromal@meta.data), value=T), 'seurat_clusters')]

# 2.6. Visualize annotations:
Seurat::DimPlot(Stromal, reduction="umap", group.by='Annotation_Level_3', label=TRUE, label.size=3, pt.size=.5, repel=T) +
  ggplot2::theme_minimal()
Seurat::DimPlot(Stromal, reduction="umap", group.by='Annotation_Level_2', label=TRUE, label.size=3, pt.size=.5, repel=T) +
  ggplot2::theme_minimal()
Seurat::DimPlot(Stromal, reduction="umap", group.by='Annotation_Level_1', label=TRUE, label.size=3, pt.size=.5, repel=T) +
  ggplot2::theme_minimal()

# 2.7. Save Seurat object
SeuratDisk::SaveH5Seurat(Stromal, paste(project_dir, '2_annotation/results_Stromal/datasets/Stromal_finalAnnots.h5Seurat', sep='/'))


# 3. Get Markers between the final annotations:
# 3.1. Annotation Level 3
Seurat::Idents(Stromal) = 'Annotation_Level_3'
Stromal_Level_3_markers = Seurat::FindAllMarkers(Stromal, assay='RNA', slot='data', logfc.threshold=.8, min.pct = 0.3, only.pos=TRUE)
View(Stromal_Level_3_markers)
write.csv(Stromal_Level_3_markers, paste(project_dir, '2_annotation/results_Stromal/markers/markers_Level_3.csv', sep='/'))
invisible(gc())

# 3.2. Annotation Level 2
Seurat::Idents(Stromal) = 'Annotation_Level_2'
Stromal_Level_2_markers = Seurat::FindAllMarkers(Stromal, assay='RNA', slot='data', logfc.threshold=.8, min.pct = 0.3, only.pos=TRUE)
View(Stromal_Level_2_markers)
write.csv(Stromal_Level_2_markers, paste(project_dir, '2_annotation/results_Stromal/markers/markers_Level_2.csv', sep='/'))
invisible(gc())

# 3.3. Annotation Level 1
Seurat::Idents(Stromal) = 'Annotation_Level_1'
Stromal_Level_1_markers = Seurat::FindAllMarkers(Stromal, assay='RNA', slot='data', logfc.threshold=.8, min.pct = 0.3, only.pos=TRUE)
View(Stromal_Level_1_markers)
write.csv(Stromal_Level_1_markers, paste(project_dir, '2_annotation/results_Stromal/markers/markers_Level_1.csv', sep='/'))
invisible(gc())


