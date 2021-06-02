project_dir = '/home/scardoso/Documents/PhD/CRC_ATLAS'

# ###################### #
# #### CNV ANALYSIS #### #
# ###################### #



# -----
# - Run Analysis
# -----

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



# ---
# - Organize data to use 
# ---

# Tranfer CNV predictions to seurat dataset
Epithelial = SeuratDisk::LoadH5Seurat(paste(project_dir, '2_annotation/data/Epithelial_integratedClusters.h5Seurat', sep='/'))
invisible(gc())
infercnv_obj = readRDS(paste(project_dir, '2_annotation/CNV/run.final.infercnv_obj', sep='/'))
cells_names = colnames(infercnv_obj@expr.data)
remove(infercnv_obj)
invisible(gc())
Epithelial = subset(Epithelial, cells = cells_names)
invisible(gc())
Epithelial = infercnv::add_to_seurat(seurat_obj = Epithelial, infercnv_output_path = paste(project_dir, '2_annotation/CNV', sep='/'))
invisible(gc())
epithelial_metadata = Epithelial@meta.data
remove(Epithelial)
invisible(gc())

# HMM predictions, from infercnv.17_HMM_predHMMi6.qnorm.hmm_mode-subclusters.observations.txt and
# infercnv.17_HMM_predHMMi6.qnorm.hmm_mode-subclusters.references.txt
predictions_obs = read.table(paste(project_dir, '2_annotation/CNV/infercnv.17_HMM_predHMMi6.qnorm.hmm_mode-subclusters.observations.txt',
                                   sep='/'),
                             header=TRUE)
predictions_obs = as(as.matrix(predictions_obs), "sparseMatrix") 
predictions_refs = read.table(paste(project_dir, '2_annotation/CNV/infercnv.17_HMM_predHMMi6.qnorm.hmm_mode-subclusters.references.txt',
                                    sep='/'),
                              header=TRUE)
predictions_refs = as(as.matrix(predictions_refs), "sparseMatrix")

# Information on which chromosomes the genes belong to:
chrs_genes = read.table(paste(project_dir, 'utils/gene_pos_wArms.txt', sep='/'), row.names=1)





# ---
# - Tumour vs Normal classification based on CNAs
# ---

# Classify patient cells (cells with more than 1/3 of genes in a relevant chr arm with the expected aberration will be considered as tumour)
cell_annotation = c()
for(cell in colnames(predictions_obs)){
  message(cell, ' ', round(which(colnames(predictions_obs)==cell) / dim(predictions_obs)[2] * 100, digits=2))
  decision = 'Putative Normal'
  for(rel_chr in relevant_chrs){
    message('-', rel_chr)
    genes_chr_arm = chrs_genes$chr_arm==rel_chr
    n_genes = sum(genes_chr_arm)
    if(names(relevant_chrs[relevant_chrs==rel_chr])=='Loss') as_expected = sum(predictions_obs[genes_chr_arm, cell] < 3)
    else as_expected = sum(predictions_obs[genes_chr_arm, cell] > 3)
    pct_expected = as_expected / n_genes
    if(pct_expected > (1/3)){
      decision = 'Putative Tumour'
      break
    }
  }
  cell_annotation = c(cell_annotation, decision)
  message('')
}
names(cell_annotation) = colnames(predictions_obs)






# ---
# - Hierarchical clusterings
# ---
heatmap_colors = colorRampPalette(c("blue", "white","red"))(n = 299)
heatmap_color_breaks = c(seq(1,3,length=150), seq(3,6,length=150))


# 1. Heatmap of predictions with all genes
conv = cbind(c('Reference', 'Patient'), c('grey', 'red'))
ronv = cbind(unique(chrs_genes$V2), rainbow(length(unique(chrs_genes$V2))))
pdf(paste(project_dir, '2_annotation/CNV/predictions_heatmaps_patient.pdf', sep='/'), width=15, height=10)
predictions_clusts = list()
for(patient in grep('N', unique(epithelial_metadata$patient), value=TRUE, invert=TRUE)){
  message(patient)
  message('- Getting patient data')
  patient_preds = predictions_obs[,gsub('-', '.', rownames(epithelial_metadata)[epithelial_metadata$state=='Tumor' &
                                                                                  epithelial_metadata$patient==patient])]
  probs = cbind(predictions_refs[,sample(1:dim(predictions_refs)[2], round(dim(patient_preds)[2]/2))], patient_preds)
  cols_column = unlist(lapply(c(rep('Reference', round(dim(patient_preds)[2]/2)), rep('Patient', dim(patient_preds)[2])),
                              FUN=function(x){conv[which(conv[,1]==x),2]}))
  rows_column = unlist(lapply(chrs_genes$V2,
                              FUN=function(x){ronv[which(ronv[,1]==x),2]}))
  message('- Performing heatmap')
  predictions_clusts[[patient]] = heatmap(t(as.matrix(probs)), Rowv=NULL, Colv=NA, RowSideColors=cols_column, ColSideColors=rows_column,
                                          distfun=parallelDist::parallelDist , hclustfun=fastcluster::hclust, scale='none',
                                          main=patient, col=heatmap_colors, breaks=heatmap_color_breaks, labRow=NA, labCol=NA,
                                          keep.dendro=TRUE)
  legend('topleft', legend = conv[,1], fill=conv[,2], title='Source (rows)', bty="n", border=F, cex=0.7)
  legend('bottomright', legend = ronv[,1], fill=ronv[,2], title='Chromosomes (cols)', bty="n", border=F, ncol=2, cex=0.7)
  legend("topright", legend=c("Complete Loss", "Neutral", "Adition of > 2 copies"),fill=c("blue", "white","red"),
         title='Values', bty="n", border=T, cex=0.7, text.width=max(strwidth(conv[,1])))
  invisible(gc())
  message('')
}
dev.off()


# 2. Heatmap of predictions with only genes from relevant chromosomes, and with annotation of patients' cells:
relevant_chrs = c('chr1p', 'chr1q', 'chr4p', 'chr4q', 'chr5q', 'chr7p', 'chr7q', 'chr8p', 'chr8q', 'chr10q', 'chr13q', 'chr14q',
                  'chr15q', 'chr17p', 'chr18q', 'chr19p', 'chr19q', 'chr20q')
names(relevant_chrs) = c('Loss', 'Gain', 'Loss', 'Loss', 'Loss', 'Gain', 'Gain', 'Loss', 'Gain', 'Loss', 'Gain', 'Loss', 'Loss', 'Loss',
                         'Loss', 'Gain', 'Gain', 'Gain')
chr_arm = paste(chrs_genes$V2, chrs_genes$arm, sep='')
chrs_genes = cbind(chrs_genes, chr_arm)
heatmap_colors = circlize::colorRamp2(c(1, 3, 6), c("blue", "white", "red"))
gene_to_use = chrs_genes$chr_arm%in%relevant_chrs
literature_colors = c(rep('Loss', sum(chrs_genes$chr_arm=='chr1p')), rep('Gain', sum(chrs_genes$chr_arm=='chr1q')),
                      rep('Loss', sum(chrs_genes$chr_arm=='chr4p')), rep('Loss', sum(chrs_genes$chr_arm=='chr4q')),
                      rep('Loss', sum(chrs_genes$chr_arm=='chr5q')), rep('Gain', sum(chrs_genes$chr_arm=='chr7p')),
                      rep('Gain', sum(chrs_genes$chr_arm=='chr7q')), rep('Loss', sum(chrs_genes$chr_arm=='chr8p')),
                      rep('Gain', sum(chrs_genes$chr_arm=='chr8q')), rep('Loss', sum(chrs_genes$chr_arm=='chr10q')),
                      rep('Gain', sum(chrs_genes$chr_arm=='chr13q')), rep('Loss', sum(chrs_genes$chr_arm=='chr14q')),
                      rep('Loss', sum(chrs_genes$chr_arm=='chr15q')), rep('Loss', sum(chrs_genes$chr_arm=='chr17p')),
                      rep('Loss', sum(chrs_genes$chr_arm=='chr18q')), rep('Gain', sum(chrs_genes$chr_arm=='chr19p')),
                      rep('Gain', sum(chrs_genes$chr_arm=='chr19q')), rep('Gain', sum(chrs_genes$chr_arm=='chr20q')))
chrs_text_colors = c('white', 'black', 'white', 'black', 'white', 'black', 'white', 'white', 'white',
                     'black', 'white', 'black', 'white', 'black', 'white', 'white', 'white', 'white')
predictions_clusts_relchrs = list()
for(patient in grep('N', unique(epithelial_metadata$patient), value=TRUE, invert=TRUE)){
  message(patient)
  message('- Getting patient data')
  patient_preds = predictions_obs[gene_to_use, gsub('-', '.', rownames(epithelial_metadata)[epithelial_metadata$state=='Tumor' &
                                                                                              epithelial_metadata$patient==patient])]
  probs = cbind(predictions_refs[gene_to_use, sample(1:dim(predictions_refs)[2], round(dim(patient_preds)[2]/2))], patient_preds)
  cell_annotations = c(rep('Normal', round(dim(patient_preds)[2]/2)),
                       cell_annotation[gsub('-', '.', rownames(epithelial_metadata)[epithelial_metadata$state=='Tumor' &
                                                                                      epithelial_metadata$patient==patient])])
  message('- Performing heatmap')
  col_annotation = ComplexHeatmap::HeatmapAnnotation(foo = ComplexHeatmap::anno_block(gp = grid::gpar(fill = 2:length(relevant_chrs)),
                                                                                      labels=relevant_chrs,
                                                                                      labels_gp=grid::gpar(col=chrs_text_colors,
                                                                                                           fontsize=7)),
                                                     Literature = literature_colors,
                                                     col= list(Literature=c('Gain'='red', 'Loss'='blue')), 
                                                     show_annotation_name=FALSE)
  row_annotation = ComplexHeatmap::rowAnnotation(Samples = c(rep('Reference', round(dim(patient_preds)[2]/2)), rep('Patient', dim(patient_preds)[2])),
                                                 Annotation_1 = cell_annotations,
                                                 col = list(Samples = c('Reference'='grey', 'Patient'='red'), 
                                                            Annotation_1 = c('Normal'='grey', 'Putative Normal'='forestgreen',
                                                                             'Putative Tumour'='chocolate2')),
                                                 show_annotation_name=FALSE)
  predictions_clusts_relchrs[[patient]] = ComplexHeatmap::Heatmap(t(as.matrix(probs)), cluster_rows=TRUE, cluster_columns=FALSE,
                                                                  col=heatmap_colors, name='Heatmap Values', show_row_names=FALSE,
                                                                  show_column_names=FALSE, column_title = patient,
                                                                  column_split=chrs_genes$chr_arm[gene_to_use], top_annotation=col_annotation,
                                                                  left_annotation=row_annotation,
                                                                  jitter=TRUE)
  invisible(gc())
  message('')
}
pdf(paste(project_dir, '2_annotation/CNV/predictions_heatmaps_patient_relevant_chrs.pdf', sep='/'), width=20, height=12)
for(patient in names(predictions_clusts_relchrs)){
  message(patient)
  ComplexHeatmap::draw(predictions_clusts_relchrs[[patient]])
}
dev.off()



# ---
# - Large-scale CNA matrix
# ---
# Cells with more than 1/3 of their genes in a chr arm with one of the aberrations - loss or gain -  will be considered as having that
# aberration


# 1. Construct matrix for cells from tumour samples:
cna_matrix = as.data.frame(matrix(rep('Neutral', dim(predictions_obs)[2]*length(unique(chrs_genes$chr_arm))), nrow=dim(predictions_obs)[2]))
colnames(cna_matrix) = unique(chrs_genes$chr_arm)
rownames(cna_matrix) = colnames(predictions_obs)
for(cell in colnames(predictions_obs)[12412:26609]){
  message(cell, ' ', round(which(colnames(predictions_obs)==cell) / dim(predictions_obs)[2] * 100, digits=2))
  for(chr in unique(chrs_genes$chr_arm)){
    message('-', chr)
    genes_chr_arm = chrs_genes$chr_arm==chr
    n_genes = sum(genes_chr_arm)
    n_gain_genes = sum(predictions_obs[genes_chr_arm, cell] > 3)
    n_loss_genes = sum(predictions_obs[genes_chr_arm, cell] < 3)
    pct_gain = n_gain_genes / n_genes
    pct_loss = n_loss_genes / n_genes
    if(sum(c(pct_gain, pct_loss) > (1/3)) == 2){
      if(max(c(pct_gain, pct_loss)) > 0.5) cna_matrix[cell, chr] = c('Gain', 'Loss')[which(c(pct_gain, pct_loss)>0.5)]
    }
    else{
      if(pct_gain > (1/3)) cna_matrix[cell, chr] = 'Gain'
      else if(pct_loss > (1/3))  cna_matrix[cell, chr] = 'Loss'
    }
    
  }
  message('')
}


# 2. Create heatmaps for each sample with all chromosomes:
heatmap_colors_discrete = circlize::colorRamp2(c(1, 2, 3), c("blue", "white", "red"))
all_chrs_text_colors = rep('black', ncol(cna_matrix))
all_chrs_text_colors[c(5, 8, 13, 16, 21, 24, 29, 32, 37)] = 'white'
predictions_clusts_largeScale_allchrs = list()
for(patient in grep('N', unique(epithelial_metadata$patient), value=TRUE, invert=TRUE)){
  message(patient)
  message('- Getting patient data')
  patient_preds = as.matrix(cna_matrix[gsub('-', '.', rownames(epithelial_metadata)[epithelial_metadata$state=='Tumor' &
                                                                                      epithelial_metadata$patient==patient]), ])
  patient_preds[patient_preds=='Gain'] = 3
  patient_preds[patient_preds=='Neutral'] = 2
  patient_preds[patient_preds=='Loss'] = 1
  patient_preds = matrix(as.numeric(patient_preds),
                         ncol=ncol(patient_preds))
  cell_annotations = cell_annotation[gsub('-', '.', rownames(epithelial_metadata)[epithelial_metadata$state=='Tumor' &
                                                                                    epithelial_metadata$patient==patient])]
  message('- Performing heatmap')
  col_annotation = ComplexHeatmap::HeatmapAnnotation(foo = ComplexHeatmap::anno_block(gp = grid::gpar(fill = 2:length(colnames(cna_matrix))),
                                                                                      labels=colnames(cna_matrix),
                                                                                      labels_gp=grid::gpar(col=all_chrs_text_colors,
                                                                                                           fontsize=7)),
                                                     show_annotation_name=FALSE)
  row_annotation = ComplexHeatmap::rowAnnotation(Annotation_1 = cell_annotations,
                                                 col = list(Annotation_1 = c('Normal'='grey', 'Putative Normal'='forestgreen',
                                                                             'Putative Tumour'='chocolate2')),
                                                 show_annotation_name=FALSE)
  predictions_clusts_largeScale_allchrs[[patient]] = ComplexHeatmap::Heatmap(patient_preds, cluster_rows=TRUE,
                                                                             cluster_columns=FALSE, col=heatmap_colors_discrete,
                                                                             name='Heatmap Values', show_row_names=FALSE,
                                                                             show_column_names=FALSE, column_title = patient,
                                                                             column_split=colnames(cna_matrix),
                                                                             top_annotation=col_annotation, left_annotation=row_annotation,
                                                                             jitter=TRUE)
  invisible(gc())
  message('')
}
pdf(paste(project_dir, '2_annotation/CNV/largeScaleCNA_matrix_patient_all_chrs.pdf', sep='/'), width=20, height=12)
for(patient in names(predictions_clusts_largeScale_allchrs)){
  message(patient)
  ComplexHeatmap::draw(predictions_clusts_largeScale_allchrs[[patient]])
}
dev.off()


# 3. Create heatmaps for each sample with only relevant chromosomes:
predictions_clusts_largeScale_relchrs = list()
for(patient in grep('N', unique(epithelial_metadata$patient), value=TRUE, invert=TRUE)){
  message(patient)
  message('- Getting patient data')
  patient_preds = as.matrix(cna_matrix[gsub('-', '.', rownames(epithelial_metadata)[epithelial_metadata$state=='Tumor' &
                                                                                      epithelial_metadata$patient==patient]), relevant_chrs])
  patient_preds[patient_preds=='Gain'] = 3
  patient_preds[patient_preds=='Neutral'] = 2
  patient_preds[patient_preds=='Loss'] = 1
  patient_preds = matrix(as.numeric(patient_preds),
                         ncol=ncol(patient_preds))
  cell_annotations = cell_annotation[gsub('-', '.', rownames(epithelial_metadata)[epithelial_metadata$state=='Tumor' &
                                                                                    epithelial_metadata$patient==patient])]
  message('- Performing heatmap')
  col_annotation = ComplexHeatmap::HeatmapAnnotation(foo = ComplexHeatmap::anno_block(gp = grid::gpar(fill = 2:length(relevant_chrs)),
                                                                                      labels=relevant_chrs,
                                                                                      labels_gp=grid::gpar(col=chrs_text_colors,
                                                                                                           fontsize=7)),
                                                     Literature = names(relevant_chrs),
                                                     col= list(Literature=c('Gain'='red', 'Loss'='blue')),
                                                     show_annotation_name=FALSE)
  row_annotation = ComplexHeatmap::rowAnnotation(Annotation_1 = cell_annotations,
                                                 col = list(Annotation_1 = c('Normal'='grey', 'Putative Normal'='forestgreen',
                                                                             'Putative Tumour'='chocolate2')),
                                                 show_annotation_name=FALSE)
  predictions_clusts_largeScale_relchrs[[patient]] = ComplexHeatmap::Heatmap(patient_preds, cluster_rows=TRUE,
                                                                             cluster_columns=FALSE, col=heatmap_colors_discrete,
                                                                             name='Heatmap Values', show_row_names=FALSE,
                                                                             show_column_names=FALSE, column_title = patient,
                                                                             column_split=relevant_chrs,
                                                                             top_annotation=col_annotation, left_annotation=row_annotation,
                                                                             jitter=TRUE)
  invisible(gc())
  message('')
}
pdf(paste(project_dir, '2_annotation/CNV/largeScaleCNA_matrix_patient_relevant_chrs.pdf', sep='/'), width=20, height=12)
for(patient in names(predictions_clusts_largeScale_relchrs)){
  message(patient)
  ComplexHeatmap::draw(predictions_clusts_largeScale_relchrs[[patient]])
}
dev.off()


# 4. Store CNA analysis information in metadata slot of Epithelial's seurat dataset:
Epithelial = SeuratDisk::LoadH5Seurat(paste(project_dir, '2_annotation/data/Epithelial_integratedClusters.h5Seurat', sep='/'))
# 4.1. CNA annotation
Epithelial[['CNA_annotation']] = rep('Normal', dim(Epithelial@meta.data)[1])
Epithelial$CNA_annotation[gsub('[.]', '-', names(cell_annotation))] = cell_annotation
# 4.2. CNA matrix data
for(chr in colnames(cna_matrix)){
  meta_name = paste('CNA', chr, sep='_')
  Epithelial[[meta_name]] = rep('Neutral', dim(Epithelial@meta.data)[1])
  Epithelial@meta.data[gsub('[.]', '-', rownames(cna_matrix)), meta_name] = cna_matrix[, chr]
}
# 4.3. Save seurat object with new information
SeuratDisk::SaveH5Seurat(Epithelial, paste(project_dir, '2_annotation/data/Epithelial_integratedClusters.h5Seurat', sep='/'))





# ---
# - Analyse CNA results' 'distribution':
# ---
Epithelial = SeuratDisk::LoadH5Seurat(paste(project_dir, '2_annotation/data/Epithelial_integratedClusters.h5Seurat', sep='/'))

# 1. Subset cells to only those from tumour samples and re-cluster them:
Epithelial_tumour = subset(Epithelial, cells=rownames(Epithelial@meta.data)[Epithelial@meta.data$state=='Tumor'])
remove(Epithelial)
invisible(gc())

# 2. Redo dimension reduction and clustering
Epithelial_tumour = Seurat::RunPCA(Epithelial_tumour, assay='integrated')
Epithelial_tumour = Seurat::FindNeighbors(Epithelial_tumour, dims=1:50)
Epithelial_tumour = Seurat::FindClusters(Epithelial_tumour, resolution=seq(0.1, 1, by=.1))
Epithelial_tumour = Seurat::RunUMAP(Epithelial_tumour, dims=1:50, reduction='pca')

# 3. Choose best resolution based solely on clustree
library(ggraph)
clust_tree = clustree::clustree(Epithelial_tumour, prefix='integrated_snn_res.')
clust_tree

# 4. Visualy assess CNA analysis
res_01 = Seurat::DimPlot(Epithelial_tumour, reduction="umap", group.by='integrated_snn_res.0.1', label=TRUE, label.size=6, pt.size=.2)
res_02 = Seurat::DimPlot(Epithelial_tumour, reduction="umap", group.by='integrated_snn_res.0.2', label=TRUE, label.size=6, pt.size=.2)
res_03 = Seurat::DimPlot(Epithelial_tumour, reduction="umap", group.by='integrated_snn_res.0.3', label=TRUE, label.size=6, pt.size=.2)
res_04 = Seurat::DimPlot(Epithelial_tumour, reduction="umap", group.by='integrated_snn_res.0.4', label=TRUE, label.size=6, pt.size=.2)
(res_01 | res_02) / (res_03 | res_04)

patient = Seurat::DimPlot(Epithelial_tumour, reduction="umap", group.by='patient', label=TRUE, label.size=3, pt.size=.2)
cna_annotation = Seurat::DimPlot(Epithelial_tumour, reduction="umap", group.by='CNA_annotation', label=F, label.size=3, pt.size=.2)
(res_01 / patient) | cna_annotation

not_help_cells = rownames(Epithelial_tumour@meta.data)[Epithelial_tumour$patient%in%c('SMC03', 'SMC10', 'SMC19', 'SMC22', 'SMC24')]
not_good_cna = Seurat::DimPlot(Epithelial_tumour, cells.highlight=not_help_cells)
(patient / cna_annotation) | not_good_cna



rgb.val <- grDevices::col2rgb('grey')
grey_transparent <- rgb(rgb.val[1], rgb.val[2], rgb.val[3],  max = 255,
             alpha = (100 - 50) * 255 / 100, names = 'greyTransparent')
invisible(grey_transparent)
color_vec = c('red', 'blue', grey_transparent)
relevant_chrs = c('chr1p', 'chr1q', 'chr4p', 'chr4q', 'chr5q', 'chr7p', 'chr7q', 'chr8p', 'chr8q', 'chr10q', 'chr13q', 'chr14q',
                  'chr15q', 'chr17p', 'chr18q', 'chr19p', 'chr19q', 'chr20q')
names(relevant_chrs) = c('Loss', 'Gain', 'Loss', 'Loss', 'Loss', 'Gain', 'Gain', 'Loss', 'Gain', 'Loss', 'Gain', 'Loss', 'Loss', 'Loss',
                         'Loss', 'Gain', 'Gain', 'Gain')
Seurat::DimPlot(Epithelial_tumour, reduction="umap", group.by=paste('CNA', relevant_chrs, sep='_'),
                label=FALSE, pt.size=.2, cols=color_vec, ncol=5)
















# 5. Try to see if cells/ patients can be divided into different signatures according to the CNAs
cna_matrix_signatures_allChrs = unique(cna_matrix)
cna_matrix_only_pt = cna_matrix[cell_annotation=='Putative Tumour',]
cna_matrix_signatures_relChrs = unique(cna_matrix_only_pt[, relevant_chrs])

signatures_relChrs = apply(cna_matrix_only_pt[, relevant_chrs], 1, function(x) paste(x, collapse='_'))
signatures_relChrs_dist = table(signatures_relChrs)
top50 = signatures_relChrs_dist[order(signatures_relChrs_dist, decreasing=TRUE)][1:50]
which(signatures_relChrs%in%names(top50))


cna_matrix_signatures_relChrs_top50 = unique(cna_matrix_only_pt[which(signatures_relChrs%in%names(top50)), relevant_chrs])
cna_matrix_signatures_relChrs_heatmap = as.matrix(cna_matrix_signatures_relChrs_top50)
cna_matrix_signatures_relChrs_heatmap[cna_matrix_signatures_relChrs_heatmap=='Gain'] = 3
cna_matrix_signatures_relChrs_heatmap[cna_matrix_signatures_relChrs_heatmap=='Neutral'] = 2
cna_matrix_signatures_relChrs_heatmap[cna_matrix_signatures_relChrs_heatmap=='Loss'] = 1
cna_matrix_signatures_relChrs_heatmap = matrix(as.numeric(cna_matrix_signatures_relChrs_heatmap),
                                               ncol=ncol(cna_matrix_signatures_relChrs_heatmap))
colnames(cna_matrix_signatures_relChrs_heatmap) = colnames(cna_matrix_signatures_relChrs)
rownames(cna_matrix_signatures_relChrs_heatmap) = rownames(cna_matrix_signatures_relChrs_top50)
col_annotation = ComplexHeatmap::HeatmapAnnotation(foo = ComplexHeatmap::anno_block(gp = grid::gpar(fill = 2:length(relevant_chrs)),
                                                                                    labels=relevant_chrs,
                                                                                    labels_gp=grid::gpar(col=chrs_text_colors,
                                                                                                         fontsize=7)),
                                                   Literature = names(relevant_chrs),
                                                   col= list(Literature=c('Gain'='red', 'Loss'='blue')),
                                                   show_annotation_name=FALSE)
#row_annotation = ComplexHeatmap::rowAnnotation(Annotation_1 = cell_annotation[rownames(cna_matrix_signatures_relChrs_heatmap)],
#                                               col = list(Annotation_1 = c('Normal'='grey', 'Putative Normal'='forestgreen',
#                                                                           'Putative Tumour'='chocolate2')),
#                                               show_annotation_name=FALSE)
right_row_annotation = ComplexHeatmap::rowAnnotation(axis_normal=ComplexHeatmap::anno_barplot(as.vector(top50/sum(signatures_relChrs_dist)),
                                                                                              bar_width=0.9),
                                                     show_annotation_name=FALSE)
ComplexHeatmap::Heatmap(cna_matrix_signatures_relChrs_heatmap, cluster_rows=FALSE, cluster_columns=FALSE,
                        col=heatmap_colors_discrete, name='Heatmap Values', show_row_names=FALSE,
                        show_column_names=FALSE, column_title = 'CNA signatures - relevant chromosomes',
                        column_split=relevant_chrs,
                        top_annotation=col_annotation, #left_annotation=row_annotation,
                        right_annotation = right_row_annotation,
                        jitter=TRUE)









# ---
# - Assess validity of classifications:
# ---
Epithelial = SeuratDisk::LoadH5Seurat(paste(project_dir, '2_annotation/data/Epithelial_integratedClusters.h5Seurat', sep='/'))

#1: Subset cells to only those from tumour samples and re-cluster them:
Epithelial_tumour = subset(Epithelial, cells=rownames(Epithelial@meta.data)[Epithelial@meta.data$state=='Tumor'])
remove(Epithelial)
invisible(gc())
Epithelial_tumour = Seurat::RunPCA(Epithelial_tumour, assay='integrated')
Epithelial_tumour = Seurat::FindNeighbors(Epithelial_tumour, dims=1:50)
Epithelial_tumour = Seurat::FindClusters(Epithelial_tumour, resolution=0.6)
Epithelial_tumour = Seurat::RunUMAP(Epithelial_tumour, dims=1:50, reduction='pca')
Seurat::DimPlot(Epithelial_tumour, reduction="umap", group.by='integrated_snn_res.0.6', label=TRUE, label.size=6)
Seurat::DimPlot(Epithelial_tumour, reduction="umap", group.by='CNV_classification', label=F, label.size=6)





















