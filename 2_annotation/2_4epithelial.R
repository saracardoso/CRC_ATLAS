project_dir = '/home/scardoso/Documents/PhD/CRC_ATLAS'
source(paste(project_dir, 'utils/modified_plots.R', sep='/'))

# ######################################### #
# #### Annotate Normal epithelial cells ### #
# ######################################### #

# Load dataset of epithelial cells:
Epithelial_normal = SeuratDisk::LoadH5Seurat(paste(project_dir, '2_annotation/results_Epithelial/datasets/Epithelial_normal.h5Seurat', sep='/'))
invisible(gc())





# -----
# - Find clusters
# -----


# 1. Run PCA:
Epithelial_normal = Seurat::RunPCA(Epithelial_normal, assay='integrated')
invisible(gc())
Seurat::DefaultAssay(Epithelial_normal)='integrated'

# 2. Choose number of PCs to use (where the elbow occurs):
pct = Epithelial_normal[["pca"]]@stdev /sum(Epithelial_normal[["pca"]]@stdev) * 100
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
# 2.3. Heatmap of the first 15 PCs
Seurat::DimHeatmap(Epithelial_normal, dims=1:15, cells=500, balanced=TRUE)


# 3. Find clusters
# 3.1. Determine the K-nearest neighbor graph
Epithelial_normal = Seurat::FindNeighbors(Epithelial_normal, dims=1:elbow)
# 3.2. Find clusters:
Epithelial_normal = Seurat::FindClusters(Epithelial_normal, resolution=seq(0.1, 1, by=.1))
# 3.3. UMAP:
Epithelial_normal = Seurat::RunUMAP(Epithelial_normal, dims=1:elbow, reduction='pca')
invisible(gc())


# 4. Choose best resolution:
library(ggraph)
clust_tree = clustree::clustree(Epithelial_normal, prefix='integrated_snn_res.')
clust_tree

# 4.1. Visualize different resolutions (first resolution was chosen to be the best for now):
clusters_01 = Seurat::DimPlot(Epithelial_normal, reduction="umap", group.by='integrated_snn_res.0.1', label=TRUE, label.size=4, pt.size=.2) +
  ggplot2::theme_minimal()
clusters_02 = Seurat::DimPlot(Epithelial_normal, reduction="umap", group.by='integrated_snn_res.0.2', label=TRUE, label.size=4, pt.size=.2) +
  ggplot2::theme_minimal()
clusters_03 = Seurat::DimPlot(Epithelial_normal, reduction="umap", group.by='integrated_snn_res.0.3', label=TRUE, label.size=4, pt.size=.2) +
  ggplot2::theme_minimal()
clusters_04 = Seurat::DimPlot(Epithelial_normal, reduction="umap", group.by='integrated_snn_res.0.4', label=TRUE, label.size=4, pt.size=.2) +
  ggplot2::theme_minimal()
((clusters_01 | clusters_02) / (clusters_03 | clusters_04)) | clust_tree

# 4.2. Plot Metrics
metrics =  c("nUMI", "nGene", "S.Score", "G2M.Score", "percent.mitochondrial_RNA")
feature_plots(Epithelial_normal, metrics, ncol=3, with_dimplot=TRUE, dimplot_group='integrated_snn_res.0.4')

# 4.3. Plot metadata
state = Seurat::DimPlot(Epithelial_normal, reduction="umap", group.by='state', label=FALSE, label.size=4) +
  ggplot2::theme_minimal() + Seurat::NoLegend()
patients = Seurat::DimPlot(Epithelial_normal, reduction="umap", group.by='patient', label=FALSE, pt.size=.1) +
  ggplot2::theme_minimal() + Seurat::NoLegend()
datasets = Seurat::DimPlot(Epithelial_normal, reduction="umap", group.by='dataset', label=FALSE, pt.size=.1) +
  ggplot2::theme_minimal() + Seurat::NoLegend()
((clusters_04 + Seurat::NoLegend()) | patients) / (state | datasets)

# 4.4. Normal/ tumour/ healthy distribution across clusters:
clusters = c()
percs = c()
fill = c()
for(cluster in unique(as.character(Epithelial_normal@meta.data[,'integrated_snn_res.0.4']))){
  x = table(Epithelial_normal$state[Epithelial_normal@meta.data[,'integrated_snn_res.0.4']==cluster])[c('Tumor', 'Normal', 'Healthy')]
  x[is.na(x)] = 0
  y = sum(x)
  percs = c(percs, x/y)
  fill = c(fill, c('Tumor', 'Normal', 'Healthy'))
  clusters = c(clusters, rep(cluster, 3))
}
dist_nt = data.frame(percs=percs, fill=fill, clusters=clusters)
dist_nt$clusters = factor(dist_nt$clusters,
                          levels=as.character(0:(length(unique(Epithelial_normal@meta.data[,'integrated_snn_res.0.4']))-1)))
dist_state = ggplot2::ggplot(dist_nt, mapping = ggplot2::aes(clusters, percs, fill=fill)) +
  ggplot2::geom_bar(position='stack', stat='identity')
dist_state

# 4.5. Patient distribution across clusters:
clusters = c()
percs = c()
fill = c()
patients = unique(Epithelial_normal$patient)
for(cluster in unique(as.character(Epithelial_normal@meta.data[,'integrated_snn_res.0.4']))){
  x = table(Epithelial_normal$patient[Epithelial_normal@meta.data[,'integrated_snn_res.0.4']==cluster])[patients]
  x[is.na(x)] = 0
  y = sum(x)
  percs = c(percs, x/y)
  fill = c(fill, patients)
  clusters = c(clusters, rep(cluster, length(patients)))
}
dist_nt = data.frame(percs=percs, fill=fill, clusters=clusters)
dist_nt$clusters = factor(dist_nt$clusters,
                          levels=as.character(0:(length(unique(Epithelial_normal@meta.data[,'integrated_snn_res.0.4']))-1)))
dist_patient = ggplot2::ggplot(dist_nt, mapping = ggplot2::aes(clusters, percs, fill=fill)) +
  ggplot2::geom_bar(position='stack', stat='identity')
dist_patient # Cluster 10 is mostly N51 cells (~95%, others are from N46, N13, 35, 38). I have decided to remove this cluster.


# 5. Remove cells from cluster 10
Epithelial_normal = subset(Epithelial_normal, cells=rownames(Epithelial_normal[[]])[Epithelial_normal$integrated_snn_res.0.4!='10'])
Seurat::DimPlot(Epithelial_normal, reduction="umap", group.by='integrated_snn_res.0.4', label=TRUE, label.size=4, pt.size=.2) +
  ggplot2::theme_minimal()


# 6. Save object:
SeuratDisk::SaveH5Seurat(Epithelial_normal, paste(project_dir, '2_annotation/results_Epithelial/datasets/Epithelial_normal_clustered.h5Seurat', sep='/'), overwrite=TRUE)




# -----
# - First annotations
# -----


# 1. Stem Cells:
feature_plots(Epithelial_normal, c('rna_LGR5', 'rna_SMOC2', 'rna_ASCL2'),
              ncol=2, with_dimplot=TRUE, dimplot_group='integrated_snn_res.0.4')
violin_plots(Epithelial_normal, c('rna_LGR5', 'rna_SMOC2', 'rna_ASCL2'),
              ncol=2, with_dimplot=TRUE, group.by='integrated_snn_res.0.4')


# 2. Tuft cells
feature_plots(Epithelial_normal, c('rna_POU2F3', 'rna_TRPM5', 'rna_SPIB', 'rna_IL17RB', 'rna_HTR3E'),
              ncol=3, with_dimplot=TRUE, dimplot_group='integrated_snn_res.0.4')
violin_plots(Epithelial_normal, c('rna_POU2F3', 'rna_TRPM5', 'rna_SPIB', 'rna_IL17RB', 'rna_HTR3E'),
             ncol=3, with_dimplot=TRUE, group.by='integrated_snn_res.0.4')


# 3. Enteroendocrine cells (EECs)
feature_plots(Epithelial_normal, c('rna_CHGA', 'rna_CHGB', 'rna_CPE', 'rna_NEUROD1', 'rna_PYY'),
              ncol=3, with_dimplot=TRUE, dimplot_group='integrated_snn_res.0.4')
violin_plots(Epithelial_normal, c('rna_CHGA', 'rna_CHGB', 'rna_CPE', 'rna_NEUROD1', 'rna_PYY'),
             ncol=3, with_dimplot=TRUE, group.by='integrated_snn_res.0.4')


# 4. Goblet cells
feature_plots(Epithelial_normal, c('rna_MUC2', 'rna_CLCA1'),
              ncol=2, with_dimplot=TRUE, dimplot_group='integrated_snn_res.0.4')
violin_plots(Epithelial_normal, c('rna_MUC2', 'rna_CLCA1'),
             ncol=2, with_dimplot=TRUE, group.by='integrated_snn_res.0.4')


# 5. Secretory Progenitors:
feature_plots(Epithelial_normal, c('rna_MUC2', 'rna_ITLN1', 'rna_CLCA1'),
              ncol=2, with_dimplot=TRUE, dimplot_group='integrated_snn_res.0.4')
violin_plots(Epithelial_normal, c('rna_MUC2', 'rna_ITLN1', 'rna_CLCA1'),
             ncol=2, with_dimplot=TRUE, group.by='integrated_snn_res.0.4')


# 6. Transit-Amplifying (TA) cells
feature_plots(Epithelial_normal, c('rna_TOP2A', 'rna_CCNA2', 'rna_MCM5',
                                   'rna_OLFM4', 'rna_SLC12A2'),
              ncol=3, with_dimplot=TRUE, dimplot_group='integrated_snn_res.0.4')
violin_plots(Epithelial_normal, c('rna_TOP2A', 'rna_CCNA2', 'rna_MCM5',
                                  'rna_OLFM4', 'rna_SLC12A2'),
             ncol=3, with_dimplot=TRUE, group.by='integrated_snn_res.0.4')


# 7. Cycling cells
feature_plots(Epithelial_normal, c('S.Score', 'G2M.Score', 'rna_PCNA', 'rna_MKI67'),
              ncol=3, with_dimplot=TRUE, dimplot_group='integrated_snn_res.0.4')
violin_plots(Epithelial_normal, c('S.Score', 'G2M.Score', 'rna_PCNA', 'rna_MKI67'),
             ncol=3, with_dimplot=TRUE, group.by='integrated_snn_res.0.4')


# 8. Enterocyte/colonocyte cells
feature_plots(Epithelial_normal, c('rna_FABP1', 'rna_SLC26A3', 'rna_TMEM37', 'rna_BEST4'),
              ncol=3, with_dimplot=TRUE, dimplot_group='integrated_snn_res.0.4')
violin_plots(Epithelial_normal, c('rna_FABP1', 'rna_SLC26A3', 'rna_TMEM37', 'rna_BEST4'),
             ncol=3, with_dimplot=TRUE, group.by='integrated_snn_res.0.4')


# 9. Paneth-like cells
feature_plots(Epithelial_normal, c('rna_LYZ', 'rna_CA7', 'rna_CA4', 'rna_SPIB', 'rna_FKBP1A'),
              ncol=3, with_dimplot=TRUE, dimplot_group='integrated_snn_res.0.4')
violin_plots(Epithelial_normal, c('rna_LYZ', 'rna_CA7', 'rna_CA4', 'rna_SPIB', 'rna_FKBP1A'),
             ncol=3, with_dimplot=TRUE, group.by='integrated_snn_res.0.4')


# 10. Progenitor cells
feature_plots(Epithelial_normal, c('rna_SOX9', 'rna_CDK6', 'rna_MUC4', 'rna_FABP5', 'rna_PLA2G2A', 'rna_LCN2'),
             ncol=3, with_dimplot=TRUE, dimplot_group='integrated_snn_res.0.4')
violin_plots(Epithelial_normal, c('rna_SOX9', 'rna_CDK6', 'rna_MUC4', 'rna_FABP5', 'rna_PLA2G2A', 'rna_LCN2'),
              ncol=3, with_dimplot=TRUE, group.by='integrated_snn_res.0.4')


# 11. Microfold (M) cells
feature_plots(Epithelial_normal, c('rna_CCL20', 'rna_SPIB', 'rna_NTRK2', 'rna_TNFRSF11A', 'rna_CCL23', 'rna_SOX8', 'rna_TNFRSF11B',
                                   'rna_FOLR3'),
              ncol=3, with_dimplot=TRUE, dimplot_group='integrated_snn_res.0.4')
violin_plots(Epithelial_normal, c('rna_CCL20', 'rna_SPIB', 'rna_NTRK2', 'rna_TNFRSF11A', 'rna_CCL23', 'rna_SOX8', 'rna_TNFRSF11B',
                                   'rna_FOLR3'),
              ncol=3, with_dimplot=TRUE, group.by='integrated_snn_res.0.4')


# 12. First Annotations
# Cluster 0      --> Immature Enterocytes
# Cluster 1      --> Progenitor cells
# Cluster 2 + 9  --> TA
# Cluster 3      --> Secretory Progenitor
# Cluster 4 + 6  --> Enterocytes
# Cluster 5      --> Immature Goblet cells
# Cluster 8      --> Stem cells
# Cluster 11     --> Goblet cells
# 12.1. Give first annotations
Epithelial_normal[['First_Annotations']] = as.character(Epithelial_normal$integrated_snn_res.0.4)
immature_enterocytes = rownames(Epithelial_normal@meta.data)[Epithelial_normal$integrated_snn_res.0.4=='0']
Epithelial_normal@meta.data[immature_enterocytes, 'First_Annotations'] = 'Immature Enterocytes'
progenitor = rownames(Epithelial_normal@meta.data)[Epithelial_normal$integrated_snn_res.0.4=='1']
Epithelial_normal@meta.data[progenitor, 'First_Annotations'] = 'Progenitor cells'
ta = rownames(Epithelial_normal@meta.data)[Epithelial_normal$integrated_snn_res.0.4%in%c('2','9')]
Epithelial_normal@meta.data[ta, 'First_Annotations'] = 'TA'
secretory = rownames(Epithelial_normal@meta.data)[Epithelial_normal$integrated_snn_res.0.4=='3']
Epithelial_normal@meta.data[secretory, 'First_Annotations'] = 'Secretory Progenitor'
enterocytes = rownames(Epithelial_normal@meta.data)[Epithelial_normal$integrated_snn_res.0.4%in%c('4', '6')]
Epithelial_normal@meta.data[enterocytes, 'First_Annotations'] = 'Enterocytes'
immature_goblet = rownames(Epithelial_normal@meta.data)[Epithelial_normal$integrated_snn_res.0.4=='5']
Epithelial_normal@meta.data[immature_goblet, 'First_Annotations'] = 'Immature Goblet cells'
stem = rownames(Epithelial_normal@meta.data)[Epithelial_normal$integrated_snn_res.0.4=='8']
Epithelial_normal@meta.data[stem, 'First_Annotations'] = 'Stem cells'
goblet = rownames(Epithelial_normal@meta.data)[Epithelial_normal$integrated_snn_res.0.4=='11']
Epithelial_normal@meta.data[goblet, 'First_Annotations'] = 'Goblet cells'
# 12.2. Visualize annotations:
Seurat::DimPlot(Epithelial_normal, reduction="umap", group.by='First_Annotations', label=TRUE, label.size=3, pt.size=.5, repel=T) +
  ggplot2::theme_minimal() + Seurat::NoLegend()





# -----
# - Sub-cluster cluster 7 from resolution 0.4 separately
# -----


# 1. Sub-cluster:
Epithelial_normal_7 = subset(Epithelial_normal, cells=rownames(Epithelial_normal[[]])[Epithelial_normal$integrated_snn_res.0.4=='7'])
(Seurat::DimPlot(Epithelial_normal, group.by = 'integrated_snn_res.0.4', label=T, label.size=4) + ggplot2::theme_minimal() + Seurat::NoLegend()) |
  (Seurat::DimPlot(Epithelial_normal_7, group.by = 'integrated_snn_res.0.4', label=T, label.size=4) + ggplot2::theme_minimal() + Seurat::NoLegend())
# 1.1. PCA
Epithelial_normal_7 = Seurat::RunPCA(Epithelial_normal_7, assay='integrated')
# 1.2. Choose number of PCs to use
pct = Epithelial_normal_7[["pca"]]@stdev /sum(Epithelial_normal_7[["pca"]]@stdev) * 100
cumu = cumsum(pct)
co1 = which(cumu > 90 & pct < 5)[1]
co2 = sort(which((pct[1:length(pct) - 1] - pct[2:length(pct)]) > 0.1), decreasing = T)[1] + 1
elbow = min(co1, co2) # 8
# 1.3. Visual representation of the PCs to use
plot_df = data.frame(pct=pct, cumu=cumu, rank=1:length(pct))
ggplot2::ggplot(plot_df, ggplot2::aes(cumu, pct, label = rank, color = rank > elbow)) + 
  ggplot2::geom_text() + 
  ggplot2::geom_vline(xintercept = 90, color = "grey") + 
  ggplot2::geom_hline(yintercept = min(pct[pct > 5]), color = "grey") +
  ggplot2::theme_bw()
# 1.4. Heatmap of the first 9 PCs
Seurat::DimHeatmap(Epithelial_normal_7, dims=1:9, cells=500, balanced=TRUE)
# 1.5. Determine the K-nearest neighbor graph
Epithelial_normal_7 = Seurat::FindNeighbors(Epithelial_normal_7, dims=1:elbow)
# 1.6. Find clusters:
Epithelial_normal_7 = Seurat::FindClusters(Epithelial_normal_7, resolution=seq(0.1, 1, by=.1))
# 1.7. UMAP:
Epithelial_normal_7 = Seurat::RunUMAP(Epithelial_normal_7, dims=1:elbow, reduction='pca')
invisible(gc())
# 1.8. Save object
SeuratDisk::SaveH5Seurat(Epithelial_normal_7, paste(project_dir, '2_annotation/results_Epithelial/datasets/Epithelial_normal_7.h5Seurat', sep='/'), overwrite=TRUE)


# 2. Choose best resolution using clustree:
library(ggraph)
clust_tree = clustree::clustree(Epithelial_normal_7, prefix='integrated_snn_res.')
clust_tree
invisible(gc())

# 2.1. Visualize different resolutions (first resolution was chosen to be the best for now):
clusters_01_0 = Seurat::DimPlot(Epithelial_normal_7, reduction="umap", group.by='integrated_snn_res.0.1', label=TRUE, label.size=6, pt.size=.2)
clusters_02_0 = Seurat::DimPlot(Epithelial_normal_7, reduction="umap", group.by='integrated_snn_res.0.2', label=TRUE, label.size=6, pt.size=.2)
clusters_05_0 = Seurat::DimPlot(Epithelial_normal_7, reduction="umap", group.by='integrated_snn_res.0.5', label=TRUE, label.size=6, pt.size=.2)
clusters_06_0 = Seurat::DimPlot(Epithelial_normal_7, reduction="umap", group.by='integrated_snn_res.0.6', label=TRUE, label.size=6, pt.size=.2)
((clusters_01_0 | clusters_02_0) / (clusters_05_0 | clusters_06_0)) | clust_tree

# 2.2. Resolution 0.1 will be chosen


# 3. Check gene markers
# 3.2. Violin Plots:
violin_plots(Epithelial_normal_7, c('rna_FABP1', 'rna_SLC26A3', 'rna_TMEM37', 'rna_BEST4',
                                    'rna_LYZ', 'rna_CA7', 'rna_CA4', 'rna_SPIB', 'rna_FKBP1A'),
             ncol=4, with_dimplot=TRUE, group.by='integrated_snn_res.0.1')


# 4. Annotation:
# Sub-cluster  0  --> BEST4+ Enterocytes
# Sub-cluster  1  --> Paneth-like cells
# 4.1. Insert the annotations in the metadata:
Epithelial_normal_7[['Final_Annotation']] = rep('', length(Epithelial_normal_7$integrated_snn_res.1))
best4 = rownames(Epithelial_normal_7@meta.data)[Epithelial_normal_7$integrated_snn_res.0.1=='1']
Epithelial_normal_7@meta.data[best4, 'Final_Annotation'] = 'BEST4+ Enterocytes'
paneth = rownames(Epithelial_normal_7@meta.data)[Epithelial_normal_7$integrated_snn_res.0.1=='0']
Epithelial_normal_7@meta.data[paneth, 'Final_Annotation'] = 'Paneth-like cells'
# 4.2. Visualize annotations:
Seurat::DimPlot(Epithelial_normal_7, group.by='Final_Annotation', pt.size=0.5, label=T, label.size=4) +
  ggplot2::ggtitle('') + ggplot2::theme_minimal() + Seurat::NoLegend()
# 4.3. Save changes:
SeuratDisk::SaveH5Seurat(Epithelial_normal_7, paste(project_dir, '2_annotation/results_Epithelial/datasets/Epithelial_normal_7.h5Seurat', sep='/'), overwrite=TRUE)





# -----
# - Sub-cluster cluster 12 from resolution 0.4 separately
# -----


# 1. Sub-cluster:
Epithelial_normal_12 = subset(Epithelial_normal, cells=rownames(Epithelial_normal[[]])[Epithelial_normal$integrated_snn_res.0.4=='12'])
(Seurat::DimPlot(Epithelial_normal, group.by = 'integrated_snn_res.0.4', label=T, label.size=4) + ggplot2::theme_minimal() + Seurat::NoLegend()) |
  (Seurat::DimPlot(Epithelial_normal_12, group.by = 'integrated_snn_res.0.4', label=T, label.size=4) + ggplot2::theme_minimal() + Seurat::NoLegend())
# 1.1. PCA
Epithelial_normal_12 = Seurat::RunPCA(Epithelial_normal_12, assay='integrated')
# 1.2. Choose number of PCs to use
pct = Epithelial_normal_12[["pca"]]@stdev /sum(Epithelial_normal_12[["pca"]]@stdev) * 100
cumu = cumsum(pct)
co1 = which(cumu > 90 & pct < 5)[1]
co2 = sort(which((pct[1:length(pct) - 1] - pct[2:length(pct)]) > 0.1), decreasing = T)[1] + 1
elbow = min(co1, co2) # 11
# 1.3. Visual representation of the PCs to use
plot_df = data.frame(pct=pct, cumu=cumu, rank=1:length(pct))
ggplot2::ggplot(plot_df, ggplot2::aes(cumu, pct, label = rank, color = rank > elbow)) + 
  ggplot2::geom_text() + 
  ggplot2::geom_vline(xintercept = 90, color = "grey") + 
  ggplot2::geom_hline(yintercept = min(pct[pct > 5]), color = "grey") +
  ggplot2::theme_bw()
# 1.4. Heatmap of the first 12 PCs
Seurat::DimHeatmap(Epithelial_normal_12, dims=1:12, cells=500, balanced=TRUE)
# 1.5. Determine the K-nearest neighbor graph
Epithelial_normal_12 = Seurat::FindNeighbors(Epithelial_normal_12, dims=1:elbow)
# 1.6. Find clusters:
Epithelial_normal_12 = Seurat::FindClusters(Epithelial_normal_12, resolution=seq(0.1, 1, by=.1))
# 1.7. UMAP:
Epithelial_normal_12 = Seurat::RunUMAP(Epithelial_normal_12, dims=1:elbow, reduction='pca')
invisible(gc())
# 1.8. Save object
SeuratDisk::SaveH5Seurat(Epithelial_normal_12, paste(project_dir, '2_annotation/results_Epithelial/datasets/Epithelial_normal_12.h5Seurat', sep='/'), overwrite=TRUE)


# 2. Choose best resolution using clustree:
library(ggraph)
clust_tree = clustree::clustree(Epithelial_normal_12, prefix='integrated_snn_res.')
clust_tree
invisible(gc())

# 2.1. Visualize different resolutions (first resolution was chosen to be the best for now):
clusters_01_0 = Seurat::DimPlot(Epithelial_normal_12, reduction="umap", group.by='integrated_snn_res.0.1', label=TRUE, label.size=6, pt.size=.2)
clusters_02_0 = Seurat::DimPlot(Epithelial_normal_12, reduction="umap", group.by='integrated_snn_res.0.2', label=TRUE, label.size=6, pt.size=.2)
clusters_06_0 = Seurat::DimPlot(Epithelial_normal_12, reduction="umap", group.by='integrated_snn_res.0.6', label=TRUE, label.size=6, pt.size=.2)
clusters_09_0 = Seurat::DimPlot(Epithelial_normal_12, reduction="umap", group.by='integrated_snn_res.0.9', label=TRUE, label.size=6, pt.size=.2)
((clusters_01_0 | clusters_02_0) / (clusters_06_0 | clusters_09_0)) | clust_tree

# 2.2. Resolution 0.1 will be chosen


# 3. Check gene markers
# 3.1. Feature Plots:
# 3.2. Violin Plots:
feature_plots(Epithelial_normal_12, c('rna_POU2F3', 'rna_TRPM5', 'rna_SPIB', 'rna_IL17RB', 'rna_HTR3E',
                                     'rna_CHGA', 'rna_CHGB', 'rna_CPE', 'rna_NEUROD1', 'rna_PYY'),
             ncol=4, with_dimplot=TRUE, dimplot_group='integrated_snn_res.0.1')
# 3.2. Violin Plots:
violin_plots(Epithelial_normal_12, c('rna_POU2F3', 'rna_TRPM5', 'rna_SPIB', 'rna_IL17RB', 'rna_HTR3E',
                                     'rna_CHGA', 'rna_CHGB', 'rna_CPE', 'rna_NEUROD1', 'rna_PYY'),
             ncol=4, with_dimplot=TRUE, group.by='integrated_snn_res.0.1')


# 4. Annotation:
# Sub-clusters 0 + 2  --> Tuft cells
# Sub-cluster  1      --> EECs
# 4.1. Insert the annotations in the metadata:
Epithelial_normal_12[['Final_Annotation']] = rep('', length(Epithelial_normal_12$integrated_snn_res.1))
tuft = rownames(Epithelial_normal_12@meta.data)[Epithelial_normal_12$integrated_snn_res.0.1%in%c('0', '2')]
Epithelial_normal_12@meta.data[tuft, 'Final_Annotation'] = 'Tuft cells'
eecs = rownames(Epithelial_normal_12@meta.data)[Epithelial_normal_12$integrated_snn_res.0.1=='1']
Epithelial_normal_12@meta.data[eecs, 'Final_Annotation'] = 'EECs'
# 4.2. Visualize annotations:
Seurat::DimPlot(Epithelial_normal_12, group.by='Final_Annotation', pt.size=0.5, label=T, label.size=4) +
  ggplot2::ggtitle('') + ggplot2::theme_minimal() + Seurat::NoLegend()
# 4.3. Save changes:
SeuratDisk::SaveH5Seurat(Epithelial_normal_12, paste(project_dir, '2_annotation/results_Epithelial/datasets/Epithelial_normal_12.h5Seurat', sep='/'), overwrite=TRUE)





# ---
# - Finish Annotations
# ---

# Epithelial_normal_110 = SeuratDisk::LoadH5Seurat(paste(project_dir, '2_annotation/results_Epithelial/datasets/Epithelial_normal_110.h5Seurat', sep='/'))
# Epithelial_normal_7 = SeuratDisk::LoadH5Seurat(paste(project_dir, '2_annotation/results_Epithelial/datasets/Epithelial_normal_7.h5Seurat', sep='/'))
# Epithelial_normal_12 = SeuratDisk::LoadH5Seurat(paste(project_dir, '2_annotation/results_Epithelial/datasets/Epithelial_normal_12.h5Seurat', sep='/'))

# 1. Give remaining annotations
# 1.1. From sub-clustering of cluster 7
best4_enterocytes = rownames(Epithelial_normal_7@meta.data)[Epithelial_normal_7$Final_Annotation=='BEST4+ Enterocytes']
Epithelial_normal@meta.data[best4_enterocytes, 'First_Annotations'] = 'BEST4+ Enterocytes'
panethlike = rownames(Epithelial_normal_7@meta.data)[Epithelial_normal_7$Final_Annotation=='Paneth-like cells']
Epithelial_normal@meta.data[panethlike, 'First_Annotations'] = 'Paneth-like cells'

# 1.2. From sub-clustering of cluster 12
tuft = rownames(Epithelial_normal_12@meta.data)[Epithelial_normal_12$Final_Annotation=='Tuft cells']
Epithelial_normal@meta.data[tuft, 'First_Annotations'] = 'Tuft cells'
eecs = rownames(Epithelial_normal_12@meta.data)[Epithelial_normal_12$Final_Annotation=='EECs']
Epithelial_normal@meta.data[eecs, 'First_Annotations'] = 'EECs'


# 2. Visualize annotations:
Seurat::DimPlot(Epithelial_normal, group.by='First_Annotations', pt.size=0.5, label=T, label.size=4) +
  ggplot2::ggtitle('') + ggplot2::theme_minimal() + Seurat::NoLegend()


# 3. Annotation Level 1
Epithelial_normal[['Annotation_Level_1']] = Epithelial_normal$First_Annotations


# 4. Annotation Level 0
Epithelial_normal[['Annotation_Level_0']] = factor(rep('Normal Epithelial cells', dim(Epithelial_normal@meta.data)[1]))


# 5. Remove [integrated...] metadata variables
Epithelial_normal@meta.data = Epithelial_normal@meta.data[, !colnames(Epithelial_normal@meta.data) %in%
                                                            c(grep('^in', colnames(Epithelial_normal@meta.data), value=T), 'seurat_clusters')]


# 6. Save Seurat object
SeuratDisk::SaveH5Seurat(Epithelial_normal, paste(project_dir, '2_annotation/results_Epithelial/datasets/Epithelial_normal_finalAnnots.h5Seurat', sep='/'))


# 7. Get Markers between the final annotations:
# 7.1. Annotation Level 1
Seurat::Idents(Epithelial_normal) = 'Annotation_Level_1'
Epithelial_normal_Level_1_markers = Seurat::FindAllMarkers(Epithelial_normal, assay='RNA', slot='data', logfc.threshold=.8, min.pct = 0.3, only.pos=TRUE)
View(Epithelial_normal_Level_1_markers)
write.csv(Epithelial_normal_Level_1_markers, paste(project_dir, '2_annotation/results_Epithelial/markers/markers_Level_1.csv', sep='/'))
invisible(gc())







