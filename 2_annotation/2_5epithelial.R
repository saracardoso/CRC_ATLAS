project_dir = '/home/scardoso/Documents/PhD/CRC_ATLAS'
source(paste(project_dir, 'utils/modified_plots.R', sep='/'))

# ######################################### #
# #### Annotate Tumour epithelial cells ### #
# ######################################### #

# Load dataset of epithelial cells:
Epithelial_tumour = SeuratDisk::LoadH5Seurat(paste(project_dir, '2_annotation/results_Epithelial/datasets/Epithelial_tumour.h5Seurat', sep='/'))
invisible(gc())





# ---
# - Set UMAP
# ---


# 1. Run PCA:
Epithelial_tumour = Seurat::RunPCA(Epithelial_tumour, assay='integrated')
invisible(gc())
Seurat::DefaultAssay(Epithelial_tumour)='integrated'

# 2. Choose number of PCs to use (where the elbow occurs):
pct = Epithelial_tumour[["pca"]]@stdev /sum(Epithelial_tumour[["pca"]]@stdev) * 100
cumu = cumsum(pct)
co1 = which(cumu > 90 & pct < 5)[1]
co2 = sort(which((pct[1:length(pct) - 1] - pct[2:length(pct)]) > 0.1), decreasing = T)[1] + 1
# 2.1 Number of PCs to use
elbow = min(co1, co2) # 18
# 2.2. Visual representation of the PCs to use
plot_df = data.frame(pct=pct, cumu=cumu, rank=1:length(pct))
ggplot2::ggplot(plot_df, ggplot2::aes(cumu, pct, label = rank, color = rank > elbow)) + 
  ggplot2::geom_text() + 
  ggplot2::geom_vline(xintercept = 90, color = "grey") + 
  ggplot2::geom_hline(yintercept = min(pct[pct > 5]), color = "grey") +
  ggplot2::theme_bw()
# 2.3. Heatmap of the first 18 PCs
Seurat::DimHeatmap(Epithelial_tumour, dims=1:18, cells=500, balanced=TRUE)


# 3. UMAP:
Epithelial_tumour = Seurat::RunUMAP(Epithelial_tumour, dims=1:elbow, reduction='pca')
invisible(gc())





# ---
# - CMS classification
# ---

# CMScaller was installed using: remotes::install_github("Lothelab/CMScaller")


# 1. Classify tumour cells individually using CMScaller
# 1.1. Get matrix (genes x cells). Only the most variable genes were selected (those present in the integrated assay) to decrease noise
Seurat::DefaultAssay(Epithelial_tumour) = 'RNA'
epithelial_cell_matrix = t(SeuratObject::FetchData(Epithelial_tumour, rownames(Epithelial_tumour@assays$integrated@data), slot='counts'))
invisible(gc())

# 1.2. Run CMScaller
res_cms_cells = CMScaller::CMScaller(emat=epithelial_cell_matrix, rowNames='symbol',  RNAseq=TRUE, FDR=0.05)
write.csv(res_cms_cells, paste(paste(project_dir, '2_annotation/results_Epithelial/CMScaller_cells.csv', sep='/')))
invisible(gc())

# 1.3. Gene Set Analysis
par.old <- par()
par(mfrow=c(1,1), mar=par.old$mar+c(0,6,0,0))
CMScaller::subCamera(epithelial_cell_matrix, res_cms_cells$prediction,
                     geneList=CMScaller::geneSets.CRC, xKey=CMScaller::fromTo(key=rownames(epithelial_cell_matrix), id.in='symbol', id.out='entrez'),
                     doVoom=TRUE)
par(mar=par.old$mar)
invisible(gc())

# 1.4. Visualize UMAP 
Epithelial_tumour[['CMScaller_cells']] = res_cms_cells[rownames(Epithelial_tumour[[]]), 'prediction']
umap_cmscells = Seurat::DimPlot(Epithelial_tumour, group.by='CMScaller_cells', pt.size=.2) + ggplot2::theme_minimal() +
  ggplot2::scale_colour_hue(na.value = grDevices::rgb(0.75, 0.75, 0.75, alpha=0.2)) + ggplot2::theme(legend.position='bottom') + 
  ggplot2::ggtitle('')

# 1.5. CMS classification across patients (no classification counts for the ratio)
clusters = c()
percs = c()
fill = c()
cms = c('CMS1', 'CMS2', 'CMS3', 'CMS4')
for(cluster in unique(as.character(Epithelial_tumour@meta.data[,'patient']))){
  x = table(Epithelial_tumour$CMScaller_cells[Epithelial_tumour@meta.data[,'patient']==cluster], useNA='always')
  sumx = sum(x)
  y = x/sumx 
  percs = c(percs, y[cms])
  fill = c(fill, cms)
  clusters = c(clusters, rep(cluster, length(cms)))
}
dist_nt = data.frame(percs=percs, fill=fill, clusters=clusters)
dist_cmscells = ggplot2::ggplot(dist_nt, mapping = ggplot2::aes(clusters, percs, fill=fill)) +
  ggplot2::geom_bar(position='stack', stat='identity') + ggplot2::ylab('Ratio') + ggplot2::xlab('Patients') + ggplot2::ylim(c(0,1)) +
  ggplot2::theme(legend.position='bottom') + ggplot2::labs(fill='') + Seurat::RotatedAxis()

# 1.6. CMS classification across samples (no classification counts for the ratio + some patients have more than one sample)
clusters = c()
percs = c()
fill = c()
cms = c('CMS1', 'CMS2', 'CMS3', 'CMS4')
for(cluster in unique(as.character(Epithelial_tumour@meta.data[,'sample']))){
  x = table(Epithelial_tumour$CMScaller_cells[Epithelial_tumour@meta.data[,'sample']==cluster], useNA='always')
  sumx = sum(x)
  y = x/sumx 
  percs = c(percs, y[cms])
  fill = c(fill, cms)
  clusters = c(clusters, rep(cluster, length(cms)))
}
dist_nt = data.frame(percs=percs, fill=fill, clusters=clusters)
dist_cmscells_samples = ggplot2::ggplot(dist_nt, mapping = ggplot2::aes(clusters, percs, fill=fill)) +
  ggplot2::geom_bar(position='stack', stat='identity') + ggplot2::ylab('Ratio') + ggplot2::xlab('Samples') + ggplot2::ylim(c(0,1)) +
  ggplot2::theme(legend.position='bottom') + ggplot2::labs(fill='') + Seurat::RotatedAxis()

# 1.7. See both plots side-by-side
dist_cmscells | umap_cmscells
dist_cmscells_samples | umap_cmscells

# 1.8. Ratio of cells classified vs not-classified (total and by patient)
library(dplyr)
n_NA = sum(is.na(res_cms_cells$prediction))
n_class = sum(!is.na(res_cms_cells$prediction))
props = round(c(n_class, n_NA) / sum(c(n_class, n_NA)) *100, 2)
rt_df_class = data.frame(values=props, class=c('Classified', 'Not Classified'))
rt_df_class <- rt_df_class %>%
  arrange(desc(class)) %>%
  mutate(lab.ypos = cumsum(values) - 0.5*values)
ggplot2::ggplot(rt_df_class,  ggplot2::aes(x = "", y = values, fill = class)) +
  ggplot2::geom_bar(width = 1, stat = "identity", color = "white") +
  ggplot2::coord_polar("y", start = 0)+
  ggplot2::geom_text(ggplot2::aes(y = lab.ypos, label = values), color = "white")+
  ggplot2::scale_fill_manual(values =  c('#009E73', 'grey')) +
  ggplot2::theme_void() + ggplot2::labs(fill='')


# 2. Classify by sample using CMScaller
# 2.1. Get matrix (genes x samples). Only the most variable genes were selected (those present in the integrated assay) to decrease noise
Seurat::DefaultAssay(Epithelial_tumour) = 'RNA'
epithelial_sample_matrix = Seurat::AggregateExpression(Epithelial_tumour, assays='RNA', group.by='sample', slot='counts')$RNA
invisible(gc())

# 2.2. Run CMScaller
res_cms_samples = CMScaller::CMScaller(emat=epithelial_sample_matrix, rowNames='symbol',  RNAseq=TRUE, FDR=0.05)
write.csv(res_cms_samples, paste(paste(project_dir, '2_annotation/results_Epithelial/CMScaller_samples.csv', sep='/')))
invisible(gc())

# 2.3. Gene Set Analysis
par.old <- par()
par(mfrow=c(1,1), mar=par.old$mar+c(0,6,0,0))
CMScaller::subCamera(epithelial_sample_matrix, res_cms_samples$prediction,
                     geneList=CMScaller::geneSets.CRC, xKey=CMScaller::fromTo(key=rownames(epithelial_sample_matrix),
                                                                              id.in='symbol', id.out='entrez'),
                     doVoom=TRUE)
par(mar=par.old$mar)
invisible(gc())
sig_values = CMScaller::subCamera(epithelial_sample_matrix, res_cms_samples$prediction,
                                  geneList=CMScaller::geneSets.CRC, xKey=CMScaller::fromTo(key=rownames(epithelial_sample_matrix),
                                                                                           id.in='symbol', id.out='entrez'),
                                  doVoom=TRUE, doPlot=FALSE)
for(cms in names(sig_values)) write.csv(sig_values[[cms]], paste(paste(project_dir, '/2_annotation/results_Epithelial/GSEA_', cms,
                                                                       '.csv', sep='')))

# 2.4. Visualize UMAP 
Epithelial_tumour[['CMScaller_samples']] = rep(NA, dim(Epithelial_tumour[[]])[1])
for(samp in unique(Epithelial_tumour$sample)){
  Epithelial_tumour$CMScaller_samples[Epithelial_tumour$sample==samp] = as.character(res_cms_samples[samp,'prediction'])
}
umap_cmssamples = Seurat::DimPlot(Epithelial_tumour, group.by='CMScaller_samples', pt.size=.2) + ggplot2::theme_minimal() +
  ggplot2::scale_colour_hue(na.value = grDevices::rgb(0.75, 0.75, 0.75, alpha=0.2)) + ggplot2::theme(legend.position='bottom') + 
  ggplot2::ggtitle('')

# 2.5. Visualize PC plot
Seurat::DimPlot(Epithelial_tumour, reduction='pca', group.by='CMScaller_samples', pt.size=.5, dims=c(2,5)) + ggplot2::theme_minimal() +
  ggplot2::scale_colour_hue(na.value = grDevices::rgb(0.75, 0.75, 0.75, alpha=0.2)) + ggplot2::theme(legend.position='bottom') + 
  ggplot2::ggtitle('')

# 2.6. Compare annotations between samples' annotations and cells' annotations
umap_cmscells | umap_cmssamples

clusters = c()
percs = c()
fill = c()
labels = c()
cms = c('CMS1', 'CMS2', 'CMS3', 'CMS4')
for(cluster in unique(as.character(Epithelial_tumour@meta.data[,'sample']))){
  x = table(Epithelial_tumour$CMScaller_cells[Epithelial_tumour@meta.data[,'sample']==cluster], useNA='always')
  sumx = sum(x)
  y = x/sumx 
  percs = c(percs, y[cms])
  fill = c(fill, cms)
  clusters = c(clusters, rep(cluster, length(cms)))
  labels = c(labels, rep(as.character(res_cms_samples[cluster,'prediction']), length(cms)))
}
dist_nt = data.frame(percs=percs, fill=fill, clusters=clusters, label=labels)
dist_labels = unique(dist_nt[,c(3,4)])
for(i in 1:length(unique(dist_nt$clusters))){
  samp = unique(dist_nt$clusters)[i]
  dist_labels$y[i] =  sum(dist_nt$percs[dist_nt$clusters==samp])
}
dist_labels$text = rep('', dim(dist_labels)[1])
barplot_compare = ggplot2::ggplot(dist_nt, mapping = ggplot2::aes(clusters, percs, fill=fill)) +
  ggplot2::geom_bar(position='stack', stat='identity') + ggplot2::ylab('Ratio') + ggplot2::xlab('Samples') + ggplot2::ylim(c(0,1)) +
  ggplot2::theme(legend.position='bottom') + ggplot2::labs(fill='') + Seurat::RotatedAxis() +
  ggplot2::geom_label(data=dist_labels, ggplot2::aes(label=text, y=y, fill=label, fontface='bold'),
                      show.legend=FALSE, size=3, nudge_y = 0.02, label.padding=ggplot2::unit(0.2, 'lines'))


# 3. Check distribution of cell-groups proportions for each sample (immune - myeloids, B and T cells - and stromal subsets).
#    CMS1 classified samples should have the higher propotions for immune subset and CMS4 should have highest proportions for stromal subset
CRCatlas_global = SeuratDisk::LoadH5Seurat(paste(project_dir, '2_annotation/results_globalAnnotation/datasets/CRC_annotations.h5Seurat', sep='/'))
invisible(gc())

# 3.1. Boxplot of the distribution of all cell-types proportions per CMS classification
samples = c()
cell_types = c()
cell_types2 = c()
ct_proportions = c()
labels = c()
cts = c('Epithelial cells', 'Stromal cells', 'Myeloid cells', 'Tcells', 'Bcells')
for(samp in unique(as.character(Epithelial_tumour@meta.data[,'sample']))){
  x = table(CRCatlas_global$Level_0[CRCatlas_global@meta.data[,'sample']==samp])[cts]
  x[is.na(x)] = 0
  sumx = sum(x)
  y = x/sumx

  samples = c(samples, rep(samp, 5))
  labels = c(labels, rep(as.character(res_cms_samples[samp,'prediction']), 5))
  cell_types = c(cell_types, cts)
  ct_proportions = c(ct_proportions, y)
}
dist_ct = data.frame(sample=samples, classification=labels, cell_type=cell_types, proportion=ct_proportions)
dist_ct$classification[is.na(dist_ct$classification)] = 'Mixed'
ct_boxplot = ggplot2::ggplot(dist_ct, ggplot2::aes(classification, proportion, colour=classification)) +
  ggplot2::geom_boxplot() + ggplot2::facet_wrap(ggplot2::vars(cell_type)) +
  ggplot2::scale_color_manual(values=c('CMS1'='#f8766d', 'CMS2'='#7cae00', 'CMS3'='#00bfc4', 'CMS4'='#c77cff',
                                       'Mixed'='#cccccc')) + Seurat::RotatedAxis()

# 3.2. Boxplot of the distribution of epithelial, stromal and immune cells proportions per CMS classification
samples = c()
cell_types = c()
ct_proportions = c()
labels = c()
cts = c('Epithelial cells', 'Stromal cells', 'Myeloid cells', 'Tcells', 'Bcells')
for(samp in unique(as.character(Epithelial_tumour@meta.data[,'sample']))){
  x = table(CRCatlas_global$Level_0[CRCatlas_global@meta.data[,'sample']==samp])[cts]
  x[is.na(x)] = 0
  sumx = sum(x)
  y = x/sumx
  
  samples = c(samples, rep(samp, 3))
  labels = c(labels, rep(as.character(res_cms_samples[samp,'prediction']), 3))
  cell_types = c(cell_types, c(cts[1:2], 'Immune cells'))
  ct_proportions = c(ct_proportions, c(y[1:2], sum(y[3:5])))
}
dist_ct2 = data.frame(sample=samples, classification=labels, cell_type=cell_types, proportion=ct_proportions)
dist_ct2$classification[is.na(dist_ct2$classification)] = 'Mixed'
ct_boxplot2 = ggplot2::ggplot(dist_ct2, ggplot2::aes(classification, proportion, colour=classification)) +
  ggplot2::geom_boxplot() + ggplot2::facet_wrap(ggplot2::vars(cell_type)) +
  ggplot2::scale_color_manual(values=c('CMS1'='#f8766d', 'CMS2'='#7cae00', 'CMS3'='#00bfc4', 'CMS4'='#c77cff',
                                'Mixed'='#cccccc')) + Seurat::RotatedAxis()


# 4. Hierarchical clustering of samples (do samples from same CMS cluster together?)
# 4.1. Calculate distance matrices between clusters:
dist_samples = dist(t(epithelial_sample_matrix), method='euclidean')
# 4.2. Create dendogram:
hclust_samples = hclust(dist_samples, method='complete')
dend_samples = as.dendrogram(hclust_samples)
classes_vals = as.character(res_cms_samples[colnames(epithelial_sample_matrix),'prediction'])
classes_vals[is.na(classes_vals)] = 'Mixed'
colors_to_use = rep('#cccccc', length(classes_vals))
colors_to_use[classes_vals=='CMS1'] = '#f8766d'
colors_to_use[classes_vals=='CMS2'] = '#7cae00'
colors_to_use[classes_vals=='CMS3'] = '#00bfc4'
colors_to_use[classes_vals=='CMS4'] = '#c77cff'
colors_to_use = factor(colors_to_use)
colors_to_use = colors_to_use[order.dendrogram(dend_samples)]
dendextend::labels_colors(dend_samples) <- as.character(colors_to_use)
plot(dend_samples)


# 5. See results plots together:
par(mfrow=c(2,2))
plot.new()
plot.new()
plot(dend_samples)
vp = grid::viewport(height=ggplot2::unit(.5, 'npc'), width=ggplot2::unit(1, 'npc'), just=c('left', 'top'), y=1, x=0)
print(barplot_compare, vp=vp)
vp = grid::viewport(height=ggplot2::unit(.5, 'npc'), width=ggplot2::unit(.5, 'npc'), just=c('left', 'top'), y=.5, x=.5)
print(ct_boxplot, vp=vp)


# 6. Set classification from using all cells together as a bulk sample as the final classifications:
Epithelial_tumour[['Annotation_Level_1']] = Epithelial_tumour$CMScaller_samples
Epithelial_tumour$Annotation_Level_1[is.na(Epithelial_tumour$Annotation_Level_1)] = 'Mixed'
# 6.1. Visualize annotations in a UMAP
Seurat::DimPlot(Epithelial_tumour, group.by='Annotation_Level_1', pt.size=.2) + ggplot2::theme_minimal() +
  ggplot2::scale_colour_hue(na.value = grDevices::rgb(0.75, 0.75, 0.75, alpha=0.2)) + ggplot2::theme(legend.position='bottom') + 
  ggplot2::ggtitle('') + ggplot2::scale_color_manual(values=c('CMS1'='#f8766d', 'CMS2'='#7cae00', 'CMS3'='#00bfc4', 'CMS4'='#c77cff',
                                                              'Mixed'='#cccccc'))

# 7. Save Seurat object
SeuratDisk::SaveH5Seurat(Epithelial_tumour, paste(project_dir, '2_annotation/results_Epithelial/datasets/Epithelial_tumour_finalAnnots.h5Seurat', sep='/'))


# 8. Get Markers between the CMS annotations:
# 8.1. Calculate markers
Seurat::Idents(Epithelial_tumour) = 'Annotation_Level_1'
Epithelial_tumour_Level_1_markers = Seurat::FindAllMarkers(Epithelial_tumour, assay='RNA', slot='data', logfc.threshold=.8,
                                                           min.pct = 0.3, only.pos=TRUE)
View(Epithelial_tumour_Level_1_markers)
write.csv(Epithelial_tumour_Level_1_markers, paste(project_dir, '2_annotation/results_Epithelial/markers/markers_tumour_Level_1.csv', sep='/'))
invisible(gc())





# ---
# - Check CNAs
# ---


# 1. Check distribution of CNAs by CMS annotation
cms_vec = c()
chr_vec = c()
type_vec = c()
props_vec = c()
for(cms in unique(Epithelial_tumour$Annotation_Level_1)){
  cms_cells = Epithelial_tumour@meta.data[Epithelial_tumour$Annotation_Level_1==cms, ]
  for(chr in grep('chr', colnames(Epithelial_tumour[[]]), value=T)){
    chr_name = gsub('CNA_', '', chr)
    chr_losses = sum(cms_cells[,chr]=='Loss') / dim(cms_cells)[1]
    chr_gains = sum(cms_cells[,chr]=='Gain') / dim(cms_cells)[1]
    chr_neutrals = sum(cms_cells[,chr]=='Neutral') / dim(cms_cells)[1]
    
    cms_vec = c(cms_vec, rep(cms, 3))
    chr_vec = c(chr_vec, rep(chr_name, 3))
    type_vec = c(type_vec, c('Loss', 'Neutral', 'Gain'))
    props_vec = c(props_vec, c(chr_losses, chr_neutrals, chr_gains))
  }
}
cna_cms_df = data.frame(CMS=cms_vec, chr_arm=chr_vec, type=type_vec, proportions=props_vec)
cna_cms_df$chr_arm = factor(cna_cms_df$chr_arm, levels=c('chr1p', 'chr1q', 'chr2p', 'chr2q', 'chr3p', 'chr3q', 'chr4p', 'chr4q',
                                                         'chr5p', 'chr5q', 'chr6p', 'chr6q', 'chr7p', 'chr7q', 'chr8p', 'chr8q',
                                                         'chr9p', 'chr9q', 'chr10p', 'chr10q', 'chr11p', 'chr11q', 'chr12p', 'chr12q',
                                                         'chr13q', 'chr14q', 'chr15q', 'chr16p', 'chr16q', 'chr17p', 'chr17q', 'chr18p',
                                                         'chr18q', 'chr19p', 'chr19q', 'chr20p', 'chr20q', 'chr21q', 'chr22q'))
cna_cms_df$type = factor(cna_cms_df$type, levels=c('Loss', 'Neutral', 'Gain'))

ggplot2::ggplot(cna_cms_df, mapping = ggplot2::aes(x=CMS, y=proportions, fill=type)) +
  ggplot2::geom_bar(stat='identity') + ggplot2::ylab('Ratio') + ggplot2::xlab('CMS') + ggplot2::facet_wrap(ggplot2::vars(chr_arm)) + 
  ggplot2::scale_fill_manual(values=c('Loss'='blue', 'Neutral'='lightgrey', 'Gain'='red')) +
  ggplot2::theme_minimal() + ggplot2::theme(legend.position='bottom') + Seurat::RotatedAxis()
