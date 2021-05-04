project_dir = '/home/scardoso/Documents/PhD/CRC_ATLAS'

# Load subset of epithelial cells:
Epithelial = SeuratDisk::LoadH5Seurat('/home/scardoso/Documents/PhD/CRC_ATLAS/2_annotation/data/epithelial_cells.h5Seurat')



# -----
# - Integrate epithelial cells
# -----

# Separate epithelial cells by dataset:
Epithelial = Seurat::SplitObject(Epithelial, split.by='dataset')
# SCTransform data:
for(dts_name in names(Epithelial)) {
  Epithelial[[dts_name]] <- Seurat::SCTransform(Epithelial[[dts_name]], vars.to.regress = c('CC.Difference', 'percent.mitochondrial_RNA'))
  gc()
}
gc()
# Select integration features and prepare for integration:
features = Seurat::SelectIntegrationFeatures(Epithelial, nfeatures = 3000)
Epithelial = Seurat::PrepSCTIntegration(Epithelial, anchor.features=features)
# Find integration achors and integrate the data:
anchors = Seurat::FindIntegrationAnchors(Epithelial, normalization.method='SCT', anchor.features=features)
gc()
Epithelial = Seurat::IntegrateData(anchorset=anchors, normalization.method="SCT")
gc()



# -----
# - Prepare data to find clusters
# -----

# Run PCA:
Epithelial = Seurat::RunPCA(Epithelial, assay='integrated')
# Find the number of PCs to use:
pct = Epithelial[["pca"]]@stdev /sum(Epithelial[["pca"]]@stdev) * 100
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

Seurat::DimHeatmap(Epithelial, dims=1:24, cells=500, balanced=TRUE)



# -----
# - Find Clusters
# -----

# Determine the K-nearest neighbor graph
Epithelial = Seurat::FindNeighbors(Epithelial, dims=1:elbow)

# Find clusters:
Epithelial = Seurat::FindClusters(Epithelial, resolution=c(0.4, 0.6, 1))

# UMAP:
Epithelial = Seurat::RunUMAP(Epithelial, dims=1:elbow, reduction='pca')

# Save object with clusters:
SeuratDisk::SaveH5Seurat(Epithelial, paste(project_dir, '2_annotation/data/Epithelial_integratedClusters.h5Seurat', sep='/'))



# -----
# - Visualize Clusters
# -----

# Cluster from the three resolutions:
clusters_04 = Seurat::DimPlot(Epithelial, reduction="umap", group.by='integrated_snn_res.0.4', label=TRUE, label.size=6)
clusters_06 = Seurat::DimPlot(Epithelial, reduction="umap", group.by='integrated_snn_res.0.6', label=TRUE, label.size=6)
clusters_1 = Seurat::DimPlot(Epithelial, reduction="umap", group.by='integrated_snn_res.1', label=TRUE, label.size=6)
(clusters_04 | clusters_06) / (clusters_1 | ggplot2::ggplot())

# Choose resolution 0.6
res = 'integrated_snn_res.0.6'

# Metadata:
state = Seurat::DimPlot(Epithelial, reduction="umap", group.by='state', label=FALSE, label.size=6)
patients = Seurat::DimPlot(Epithelial, reduction="umap", group.by='patient', label=FALSE, label.size=6, pt.size=.1)
datasets = Seurat::DimPlot(Epithelial, reduction="umap", group.by='dataset', label=FALSE, label.size=6, pt.size=.1)
(clusters_06 | patients) / (state | datasets)

# Normal/ tumour/ healthy distribution across clusters:
clusters = c()
percs = c()
fill = c()
for(cluster in unique(as.character(Epithelial@meta.data[,res]))){
  x = table(Epithelial$state[Epithelial@meta.data[,res]==cluster])[c('Tumor', 'Normal', 'Healthy')]
  x[is.na(x)] = 0
  y = sum(x)
  percs = c(percs, x/y)
  fill = c(fill, c('Tumor', 'Normal', 'Healthy'))
  clusters = c(clusters, rep(cluster, 3))
}
dist_nt = data.frame(percs=percs, fill=fill, clusters=clusters)
dist_nt$clusters = factor(dist_nt$clusters,
                          levels=as.character(0:(length(unique(Epithelial@meta.data[,res]))-1)))
dist_state = ggplot2::ggplot(dist_nt, mapping = ggplot2::aes(clusters, percs, fill=fill)) +
  ggplot2::geom_bar(position='stack', stat='identity')
dist_state

# Some gene expression:
vlngenePlot = Seurat::VlnPlot(Epithelial, features=c('rna_MYC', 'rna_AXIN2', 'rna_RNF43', 'rna_AREG', 'rna_FABP2',
                                                  'rna_ZG16', 'rna_CPE', 'rna_POU2F3', 'rna_PCNA', 'rna_SOX9', 'rna_LGR5'),
                              pt.size = 0, ncol=2, group.by = res)
clusters_06 | vlngenePlot

vlngenePlot2 = Seurat::VlnPlot(Epithelial, features=c('rna_MYC', 'rna_AXIN2', 'rna_RNF43', 'rna_AREG'),
                              pt.size = 0, ncol=2, group.by = res)
dotplot = Seurat::DotPlot(Epithelial, features=c('rna_MYC', 'rna_AXIN2', 'rna_RNF43', 'rna_AREG'), group.by = res)
rgplot = Seurat::RidgePlot(Epithelial, features=c('rna_MYC', 'rna_AXIN2', 'rna_RNF43', 'rna_AREG'), ncol=2)
(clusters_06 | vlngenePlot2) / (dotplot | rgplot)

# Cells with expression in all oncogenes:
Seurat::FeatureScatter(Epithelial, feature1='rna_MYC', feature2 = 'rna_AREG', group.by='state')
myc_areg_cells = subset(Epithelial, subset= rna_MYC>0 & rna_AREG>0)
sort(table(myc_areg_cells$state) / sum(table(myc_areg_cells$integrated_snn_res.0.6)) * 100, decreasing=TRUE)
sort(table(myc_areg_cells$integrated_snn_res.0.6) / sum(table(myc_areg_cells$integrated_snn_res.0.6)) * 100, decreasing=TRUE)
Seurat::FeatureScatter(myc_areg_cells, feature1='rna_AXIN2', feature2 = 'rna_RNF43', group.by='state')
Seurat::FeatureScatter(myc_areg_cells, feature1='rna_AXIN2', feature2 = 'rna_RNF43', group.by=res)
myc_areg_axin_rnf_cells = subset(Epithelial, subset= rna_AXIN2>0 & rna_RNF43>0)
sort(table(myc_areg_axin_rnf_cells$state) / sum(table(myc_areg_axin_rnf_cells$integrated_snn_res.0.6)) * 100, decreasing=TRUE)
sort(table(myc_areg_axin_rnf_cells$integrated_snn_res.0.6) / sum(table(myc_areg_axin_rnf_cells$integrated_snn_res.0.6)) * 100, decreasing=TRUE)
cancer_genes_cells = colnames(myc_areg_axin_rnf_cells@assays$RNA@counts)

# Highlight of cells with all oncogenes expressed:
clusters_06 | Seurat::DimPlot(Epithelial, cells.highlight = cancer_genes_cells)
# Distribution of cells with all oncogenes expressed accross clusters:
Epithelial[['cells_all_oncogenes']] = rownames(Epithelial@meta.data) %in% cancer_genes_cells
clusters = c()
percs = c()
fill = c()
for(cluster in unique(as.character(Epithelial@meta.data[,res]))){
  x = table(Epithelial$cells_all_oncogenes[Epithelial@meta.data[,res]==cluster])[c('FALSE', 'TRUE')]
  x[is.na(x)] = 0
  y = sum(x)
  percs = c(percs, x/y)
  fill = c(fill, c('FALSE', 'TRUE'))
  clusters = c(clusters, rep(cluster, 2))
}
dist_nt = data.frame(percs=percs, fill=fill, clusters=clusters)
dist_nt$clusters = factor(dist_nt$clusters,
                          levels=as.character(0:(length(unique(Epithelial@meta.data[,res]))-1)))
dist_oncogenes = ggplot2::ggplot(dist_nt, mapping = ggplot2::aes(clusters, percs, fill=fill)) +
  ggplot2::geom_bar(position='stack', stat='identity')
dist_oncogenes





# -----
# - Find markers
# -----

Seurat::Idents(Epithelial) = res
Epithelial_markers = Seurat::FindAllMarkers(Epithelial, assay='RNA', slot='data', logfc.threshold=1, only.pos=TRUE)
View(Epithelial_markers)
write.csv(Epithelial_markers, paste(project_dir, '2_annotation/markers/Epithelial/Epithelial_res06.csv', sep='/'))

# Get top 10 markers for each cluster:
Epithelial_markers_top20 = dplyr::top_n(dplyr::group_by(Epithelial_markers[Epithelial_markers$p_val_adj<0.05,], cluster), n=20, wt=avg_log2FC)
View(Epithelial_markers_top20)
write.csv(Epithelial_markers_top20, paste(project_dir, '2_annotation/markers/Epithelial/Epithelial_res06_top20.csv', sep='/'))

# Clusters highlighted in the full atlas:
CRCatlas_integrated = SeuratDisk::LoadH5Seurat('/home/scardoso/Documents/PhD/CRC_ATLAS/2_annotation/data/CRC_annotations.h5Seurat')
for(cluster in as.character(0:19)){
  cluster_cells = rownames(Epithelial@meta.data)[Epithelial@meta.data$integrated_snn_res.0.6==cluster]
  grDevices::png(filename = paste(project_dir, '/2_annotation/figures/Epithelial/Cluster',
                                  cluster, '.png', sep=''), width=805, height=461)
  print(Seurat::DimPlot(CRCatlas_integrated, cells.highlight = cluster_cells) +
          ggplot2::ggtitle(paste('Cluster', cluster)) + Seurat::NoLegend())
  dev.off()
}

# Heatmap with the gene markers:
Epithelial = Seurat::NormalizeData(Epithelial, assay='RNA')
Epithelial = Seurat::ScaleData(Epithelial, assay='RNA')
gene_markers = unique(Epithelial_markers_top20$gene)
sub_cells = subset(Epithelial, downsample=600)
sub_cells = Seurat::NormalizeData(sub_cells, assay='RNA')
sub_cells = Seurat::ScaleData(sub_cells, assay='RNA')
Seurat::DoHeatmap(sub_cells, features=gene_markers[1:60], assay='RNA')
Seurat::DoHeatmap(sub_cells, features=gene_markers[61:120], assay='RNA')
Seurat::DoHeatmap(sub_cells, features=gene_markers[121:180], assay='RNA')
Seurat::DoHeatmap(sub_cells, features=gene_markers[181:200], assay='RNA')
Seurat::DoHeatmap(sub_cells, features=gene_markers[201:274], assay='RNA')





# -----
# - CNV Analysis
# -----
project_dir = '/home/scardoso/Documents/PhD/CRC_ATLAS'
#load(paste(project_dir, '2_annotation/CNV/cnv_environment.RData', sep='/'))

Epithelial = SeuratDisk::LoadH5Seurat(paste(project_dir, '2_annotation/data/Epithelial_integratedClusters.h5Seurat', sep='/'))
epithelial_m = Seurat::GetAssayData(Epithelial, slot='counts', assay='RNA')
epithelial_meta = Epithelial@meta.data
epithelial_meta[['state2']] = epithelial_meta$state
epithelial_meta$state2[epithelial_meta$state2=='Healthy'] = 'Normal'
epithelial_meta[epithelial_meta$state2=='Tumor', 'state2'] = epithelial_meta[epithelial_meta$state2=='Tumor', 'patient']
ids_use = c(sample(which(epithelial_meta$state=='Healthy'), 2500), sample(which(epithelial_meta$state=='Normal'), 2500),
            which(epithelial_meta$state=='Tumor'))
epithelial_meta = epithelial_meta[ids_use,]
epithelial_m = epithelial_m[, ids_use]
write.table(cbind(rownames(epithelial_meta), epithelial_meta[,'state2']),
            paste(project_dir, '2_annotation/CNV/annotation.txt', sep='/'), sep='\t', col.names=FALSE, row.names=FALSE)

# Create inferCNV object:
infercnv_obj = infercnv::CreateInfercnvObject(epithelial_m,
                                              annotations_file=paste(project_dir, '2_annotation/CNV/annotation.txt', sep='/'),
                                              gene_order_file = paste(project_dir, 'utils/gencode_v19_gene_pos.txt', sep='/'),
                                              delim='\t', ref_group_names = c('Normal'))
invisible(gc())

# Run inferCNV:
infercnv_obj = infercnv::run(infercnv_obj, cutoff=0.1, out_dir=paste(project_dir, '2_annotation/CNV', sep='/'),
                             cluster_by_groups=T, cluster_references=F,
                             analysis_mode='subclusters', tumor_subcluster_partition_method = 'qnorm',
                             denoise=T, HMM=T, scale_data=T,
                             BayesMaxPNormal=0)

# Tranfer CNV predictions to seurat dataset
Epithelial = SeuratDisk::LoadH5Seurat(paste(project_dir, '2_annotation/data/Epithelial_integratedClusters.h5Seurat', sep='/'))
infercnv_obj = readRDS(paste(project_dir, '2_annotation/CNV/run.final.infercnv_obj', sep='/'))
cells_names = colnames(infercnv_obj@expr.data)
remove(infercnv_obj)
invisible(gc())
Epithelial = subset(Epithelial, cells = cells_names)
invisible(gc())
Epithelial = infercnv::add_to_seurat(seurat_obj = Epithelial, infercnv_output_path = paste(project_dir, '2_annotation/CNV', sep='/'))
invisible(gc())
#SeuratDisk::SaveH5Seurat(Epithelial, paste(project_dir, '2_annotation/CNV/Epithelial_cnvs.h5Seurat', sep='/'))


# Presence of relevant CNVs in CRC literature:
Seurat::FeaturePlot(Epithelial, features=c('has_dupli_chr7', 'has_dupli_chr8', 'has_dupli_chr13', 'has_dupli_chr19', 'has_dupli_chr20'),
                    cols=c("lightgrey", "red")) |
  Seurat::FeaturePlot(Epithelial, features=c('has_loss_chr1', 'has_loss_chr4', 'has_loss_chr5', 'has_loss_chr8', 'has_loss_chr10',
                                             'has_loss_chr14', 'has_loss_chr15', 'has_loss_chr17', 'has_loss_chr18'),
                      cols=c("lightgrey", "blue"))

total_cnvs = Epithelial@meta.data[,c(grep('has_loss_chr', colnames(Epithelial@meta.data), value=TRUE),
                                     grep('has_dupli_chr', colnames(Epithelial@meta.data), value=TRUE))]
total_cnvs[total_cnvs==TRUE] = 1
total_cnvs[total_cnvs==FALSE] = 0
total_cnvs_tmr = total_cnvs[rownames(Epithelial@meta.data)[Epithelial@meta.data$state=='Tumor'],]

relevant_cnvs = Epithelial@meta.data[,c('has_dupli_chr7', 'has_dupli_chr8', 'has_dupli_chr13', 'has_dupli_chr19', 'has_dupli_chr20',
                                        'has_loss_chr1', 'has_loss_chr4', 'has_loss_chr5', 'has_loss_chr8', 'has_loss_chr10',
                                        'has_loss_chr14', 'has_loss_chr15', 'has_loss_chr17', 'has_loss_chr18')]
relevant_cnvs[relevant_cnvs==TRUE] = 1
relevant_cnvs[relevant_cnvs==FALSE] = 0
relevant_cnvs_tmr = relevant_cnvs[rownames(Epithelial@meta.data)[Epithelial@meta.data$state=='Tumor'],]

plot(table(rowSums(total_cnvs_tmr)), lty=1, type='l', xlab='Nº of Chromosomes with CNVs', ylab='Nº Cells', main='All Cells')
plot(table(rowSums(relevant_cnvs_tmr)), lty=1, type='l', xlab='Nº of Chromosomes with Relevant CNVs', ylab='Nº Cells', main='All Cells')

# Normal Cells:
normal_cells = rownames(Epithelial@meta.data)[Epithelial@meta.data$state!='Tumor']
# Putative normal cells 1
putative_normal_1 = rownames(total_cnvs_tmr)[rowSums(total_cnvs_tmr)==0]
# Putative tumour:
putative_tumour = rownames(relevant_cnvs_tmr)[rowSums(relevant_cnvs_tmr)>0]
# Putative normal cells 2:
putative_normal_2 = rownames(relevant_cnvs_tmr)[!rownames(relevant_cnvs_tmr)%in%c(normal_cells, putative_normal_1, putative_tumour)]

# Nº of CNVs present in the cells with no relevant CNV:
plot(table(rowSums(total_cnvs_tmr[putative_normal_2,])), lty=1, type='l',
     xlab='Nº of Chromosomes with CNVs present', ylab='Nº Cells', main='Putative Normal cells')

# Initial CNV characterization
Epithelial[['cnv_groups_initial']] = rep('Normal', dim(Epithelial@meta.data)[1])
Epithelial$cnv_groups_initial[putative_normal_1] = 'putative_normal_1'
Epithelial$cnv_groups_initial[putative_tumour] = 'putative_tumour'
Epithelial$cnv_groups_initial[putative_normal_2] = 'putative_normal_2'
Seurat::DimPlot(Epithelial, group.by = 'cnv_groups_initial', pt.size = 0.3)
Seurat::DimPlot(Epithelial, split.by = 'cnv_groups_initial', pt.size = 0.3)

Seurat::DimPlot(Epithelial, group.by = 'original_annotation_1', split.by = 'cnv_groups_initial', pt.size = 0.3)


# Cluster Analysis:
observations_cnv = read.table(paste(project_dir, '2_annotation/CNV/infercnv.17_HMM_predHMMi6.qnorm.hmm_mode-subclusters.observations.txt', sep='/'), header=TRUE)
references_cnv = read.table(paste(project_dir, '2_annotation/CNV/infercnv.17_HMM_predHMMi6.qnorm.hmm_mode-subclusters.references.txt', sep='/'), header=TRUE)
conv = cbind(unique(Epithelial$cnv_groups_initial), rainbow(length(unique(Epithelial$cnv_groups_initial))))
cols_column = unlist(lapply(Epithelial$cnv_groups_initial, FUN=function(x){conv[which(conv[,1]==x),2]}))

data_matrix_p31 = cbind(references_cnv,
                        observations_cnv[,rownames(Epithelial@meta.data[Epithelial@meta.data$state=='Tumor' & Epithelial@meta.data$patient=='31',])])
cols_column_31 = cols_column[c(normal_cells,
                               rownames(Epithelial@meta.data[Epithelial@meta.data$state=='Tumor' & Epithelial@meta.data$patient=='31',]))]
heatmap(t(as.matrix(data_matrix_p31)), Rowv = NULL, Colv = NA, RowSideColors = cols_column_31,
        distfun = parallelDist::parallelDist , hclustfun = fastcluster::hclust, scale='none')
legend('right', legend = conv[,1], fill=conv[,2])

