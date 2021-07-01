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
#relevant_chrs = c('chr1p', 'chr1q', 'chr4p', 'chr4q', 'chr5q', 'chr6p', 'chr6q', 'chr7p', 'chr7q', 'chr8p', 'chr8q', 'chr11p', 'chr11q',
#                  'chr12p', 'chr12q', 'chr13q', 'chr14q', 'chr15q', 'chr17p', 'chr18p', 'chr18q', 'chr19p', 'chr19q', 'chr20p', 'chr20q',
#                  'chr21q', 'chr22q')
#names(relevant_chrs) = c('Loss', 'Gain', 'Loss', 'Loss', 'Loss', 'Gain', 'Gain', 'Gain', 'Gain', 'Loss', 'Both', 'Loss', 'Loss', 'Gain',
#                         'Gain', 'Gain', 'Loss', 'Loss', 'Loss', 'Loss', 'Loss', 'Both', 'Both', 'Both', 'Gain', 'Loss', 'Both')


# 1. Construct large-scale CNA matrix for cells from tumour samples:
#    Cells with more than 1/3 of their genes in a chr arm with one of the aberrations - loss or gain -  will be considered as having that
#    aberration
cna_matrix = as.data.frame(matrix(rep('Neutral', dim(predictions_obs)[2]*length(unique(chrs_genes$chr_arm))), nrow=dim(predictions_obs)[2]))
colnames(cna_matrix) = unique(chrs_genes$chr_arm)
rownames(cna_matrix) = colnames(predictions_obs)
for(chr in unique(chrs_genes$chr_arm)){
  message('-', chr)
  chr_gain_decision = rep('Neutral', dim(predictions_obs)[2])
  chr_loss_decision = rep('Neutral', dim(predictions_obs)[2])
  
  genes_chr_arm = chrs_genes$chr_arm==chr
  n_genes = sum(genes_chr_arm)
  n_gain_genes = Matrix::colSums(predictions_obs[genes_chr_arm, ] > 3)
  n_loss_genes = Matrix::colSums(predictions_obs[genes_chr_arm, ] < 3)
  pct_gain = n_gain_genes / n_genes
  pct_loss = n_loss_genes / n_genes
  
  chr_gain_decision[pct_gain>(1/3)] = 'Gain'
  chr_loss_decision[pct_loss>(1/3)] = 'Loss'
  
  chr_combined_decision = paste(chr_gain_decision, chr_loss_decision, sep='_')
  chr_combined_decision[chr_combined_decision=='Neutral_Neutral'] = 'Neutral'
  chr_combined_decision[chr_combined_decision=='Neutral_Loss'] = 'Loss'
  chr_combined_decision[chr_combined_decision=='Gain_Neutral'] = 'Gain'
  chr_combined_decision[chr_combined_decision=='Gain_Loss'] = 'Both'
  
  cna_matrix[,chr] = chr_combined_decision
}


# 2. Classify patient cells (cells with an expected aberration in a relevant chr arm will be considered as tumour)
cell_annotation = rep('Putative Normal', dim(cna_matrix)[1])
#sub_cna = matrix(rep(FALSE, dim(cna_matrix)[1]*dim(cna_matrix)[2]), ncol=dim(cna_matrix)[2])
#colnames(sub_cna) = colnames(cna_matrix)
#rownames(sub_cna) = rownames(cna_matrix)
#for(idx in 1:length(relevant_chrs)){
#  chr = relevant_chrs[idx]
#  expected = names(relevant_chrs)[idx]
#  
#  if(expected == 'Both') sub_cna[, chr] = (cna_matrix[, chr] != 'Neutral')
#  else sub_cna[, chr] = (cna_matrix[, chr] == expected | cna_matrix[, chr] == 'Both')
#  
#}
#cell_annotation[rowSums(sub_cna) > 0] = 'Putative Tumour'
cell_annotation[rowSums(cna_matrix!='Neutral') > 0] = 'Putative Tumour'
names(cell_annotation) = rownames(cna_matrix)







# ---
# - Hierarchical clusterings of inferCNV prediction values
# ---
all_chrs_text_colors = rep('black', ncol(cna_matrix))
all_chrs_text_colors[c(5, 8, 13, 16, 21, 24, 29, 32, 37)] = 'white'
heatmap_colors = circlize::colorRamp2(c(1, 3, 6), c("blue", "white", "red"))

# 1. Heatmap of predictions with all genes
predictions_clusts = list()
for(patient in grep('N', unique(epithelial_metadata$patient), value=TRUE, invert=TRUE)){
  message(patient)
  
  message('- Getting patient data')
  patient_preds = predictions_obs[, gsub('-', '.', rownames(epithelial_metadata)[epithelial_metadata$state=='Tumor' &
                                                                                   epithelial_metadata$patient==patient])]
  probs = as.matrix(cbind(predictions_refs[, sample(1:dim(predictions_refs)[2], round(dim(patient_preds)[2]/2))], patient_preds))
  cell_annotations = c(rep('Normal', round(dim(patient_preds)[2]/2)),
                       cell_annotation[gsub('-', '.', rownames(epithelial_metadata)[epithelial_metadata$state=='Tumor' &
                                                                                      epithelial_metadata$patient==patient])])
  message('- Performing heatmap')
  col_annotation = ComplexHeatmap::HeatmapAnnotation(foo = ComplexHeatmap::anno_block(gp = grid::gpar(fill = 2:length(unique(chrs_genes$chr_arm))),
                                                                                      labels=unique(chrs_genes$chr_arm),
                                                                                      labels_gp=grid::gpar(col=all_chrs_text_colors,
                                                                                                           fontsize=7)),
                                                     show_annotation_name=FALSE)
  row_annotation = ComplexHeatmap::rowAnnotation(Annotation_1 = cell_annotations,
                                                 col = list(Annotation_1 = c('Normal'='grey', 'Putative Normal'='forestgreen',
                                                                             'Putative Tumour'='chocolate2')),
                                                 show_annotation_name=FALSE)
  predictions_clusts[[patient]] = ComplexHeatmap::Heatmap(t(probs), cluster_rows=TRUE,
                                                          cluster_columns=FALSE, col=heatmap_colors,
                                                          name='Heatmap Values', show_row_names=FALSE,
                                                          show_column_names=FALSE, column_title = patient,
                                                          column_split=factor(chrs_genes$chr_arm, unique(chrs_genes$chr_arm)),
                                                          top_annotation=col_annotation, left_annotation=row_annotation,
                                                          jitter=TRUE)
  
  invisible(gc())
  message('')
}
pdf(paste(project_dir, '2_annotation/CNV/predictions_heatmaps_patient.pdf', sep='/'), width=20, height=12)
for(patient in names(predictions_clusts)){
  message(patient)
  ComplexHeatmap::draw(predictions_clusts[[patient]])
}
dev.off()


# 2. Heatmap of predictions with only genes from relevant chromosomes, and with annotation of patients' cells:
#chr_arm = paste(chrs_genes$V2, chrs_genes$arm, sep='')
#chrs_genes = cbind(chrs_genes, chr_arm)
#heatmap_colors = circlize::colorRamp2(c(1, 3, 6), c("blue", "white", "red"))
#gene_to_use = chrs_genes$chr_arm%in%relevant_chrs
#literature_colors = c(rep('Loss', sum(chrs_genes$chr_arm=='chr1p')), rep('Gain', sum(chrs_genes$chr_arm=='chr1q')),
#                      rep('Loss', sum(chrs_genes$chr_arm=='chr4p')), rep('Loss', sum(chrs_genes$chr_arm=='chr4q')),
#                      rep('Loss', sum(chrs_genes$chr_arm=='chr5q')), rep('Gain', sum(chrs_genes$chr_arm=='chr6p')),
#                      rep('Gain', sum(chrs_genes$chr_arm=='chr6q')), rep('Gain', sum(chrs_genes$chr_arm=='chr7p')),
#                      rep('Gain', sum(chrs_genes$chr_arm=='chr7q')), rep('Loss', sum(chrs_genes$chr_arm=='chr8p')),
#                      rep('Both', sum(chrs_genes$chr_arm=='chr8q')), rep('Loss', sum(chrs_genes$chr_arm=='chr11p')),
#                      rep('Loss', sum(chrs_genes$chr_arm=='chr11q')), rep('Gain', sum(chrs_genes$chr_arm=='chr12p')),
#                      rep('Gain', sum(chrs_genes$chr_arm=='chr12q')), rep('Gain', sum(chrs_genes$chr_arm=='chr13q')),
#                      rep('Loss', sum(chrs_genes$chr_arm=='chr14q')), rep('Loss', sum(chrs_genes$chr_arm=='chr15q')),
#                      rep('Loss', sum(chrs_genes$chr_arm=='chr17p')), rep('Loss', sum(chrs_genes$chr_arm=='chr18p')),
#                      rep('Loss', sum(chrs_genes$chr_arm=='chr18q')), rep('Both', sum(chrs_genes$chr_arm=='chr19p')),
#                      rep('Both', sum(chrs_genes$chr_arm=='chr19q')), rep('Both', sum(chrs_genes$chr_arm=='chr20p')),
#                      rep('Gain', sum(chrs_genes$chr_arm=='chr20q')), rep('Loss', sum(chrs_genes$chr_arm=='chr21q')),
#                      rep('Both', sum(chrs_genes$chr_arm=='chr22q')))
##chrs_text_colors = c('white', 'black', 'white', 'black', 'white', 'black', 'white', 'white', 'white',
##                     'black', 'white', 'black', 'white', 'black', 'white', 'white', 'white', 'white')
#predictions_clusts_relchrs = list()
#for(patient in grep('N', unique(epithelial_metadata$patient), value=TRUE, invert=TRUE)){
#  message(patient)
#  message('- Getting patient data')
#  patient_preds = predictions_obs[gene_to_use, gsub('-', '.', rownames(epithelial_metadata)[epithelial_metadata$state=='Tumor' &
#                                                                                              epithelial_metadata$patient==patient])]
#  probs = cbind(predictions_refs[gene_to_use, sample(1:dim(predictions_refs)[2], round(dim(patient_preds)[2]/2))], patient_preds)
#  cell_annotations = c(rep('Normal', round(dim(patient_preds)[2]/2)),
#                       cell_annotation[gsub('-', '.', rownames(epithelial_metadata)[epithelial_metadata$state=='Tumor' &
#                                                                                      epithelial_metadata$patient==patient])])
#  message('- Performing heatmap')
#  col_annotation = ComplexHeatmap::HeatmapAnnotation(foo = ComplexHeatmap::anno_block(gp = grid::gpar(fill = 2:length(relevant_chrs)),
#                                                                                      labels=relevant_chrs,
#                                                                                      labels_gp=grid::gpar(#col=chrs_text_colors,
#                                                                                                           fontsize=7)),
#                                                     Literature = literature_colors,
#                                                     col= list(Literature=c('Gain'='red', 'Loss'='blue', 'Both'='grey')), 
#                                                     show_annotation_name=FALSE)
#  row_annotation = ComplexHeatmap::rowAnnotation(Samples = c(rep('Reference', round(dim(patient_preds)[2]/2)), rep('Patient', dim(patient_preds)[2])),
#                                                 Annotation_1 = cell_annotations,
#                                                 col = list(Samples = c('Reference'='grey', 'Patient'='red'), 
#                                                            Annotation_1 = c('Normal'='grey', 'Putative Normal'='forestgreen',
#                                                                             'Putative Tumour'='chocolate2')),
#                                                 show_annotation_name=FALSE)
#  predictions_clusts_relchrs[[patient]] = ComplexHeatmap::Heatmap(t(as.matrix(probs)), cluster_rows=TRUE, cluster_columns=FALSE,
#                                                                  col=heatmap_colors, name='Heatmap Values', show_row_names=FALSE,
#                                                                  show_column_names=FALSE, column_title = patient,
#                                                                  column_split=factor(chrs_genes$chr_arm[gene_to_use],
#                                                                                      unique(chrs_genes$chr_arm[gene_to_use])),
#                                                                  top_annotation=col_annotation,
#                                                                  left_annotation=row_annotation,
#                                                                  jitter=TRUE)
#  invisible(gc())
#  message('')
#}
#pdf(paste(project_dir, '2_annotation/CNV/predictions_heatmaps_patient_relevant_chrs.pdf', sep='/'), width=20, height=12)
#for(patient in names(predictions_clusts_relchrs)){
#  message(patient)
#  ComplexHeatmap::draw(predictions_clusts_relchrs[[patient]])
#}
#dev.off()



# ---
# - Hierarchical clusterings of large-scale CNA matrix
# ---

# 1. Create heatmaps for each sample with all chromosomes:
heatmap_colors_discrete = circlize::colorRamp2(c(-1, 0, 1, 2), c("blue", "white", "red", "grey"))
predictions_clusts_largeScale_allchrs = list()
for(patient in grep('N', unique(epithelial_metadata$patient), value=TRUE, invert=TRUE)){
  message(patient)
  message('- Getting patient data')
  patient_preds = as.matrix(cna_matrix[gsub('-', '.', rownames(epithelial_metadata)[epithelial_metadata$state=='Tumor' &
                                                                                      epithelial_metadata$patient==patient]), ])
  patient_preds[patient_preds=='Both'] = 2
  patient_preds[patient_preds=='Gain'] = 1
  patient_preds[patient_preds=='Neutral'] = 0
  patient_preds[patient_preds=='Loss'] = -1
  patient_preds = matrix(as.numeric(patient_preds),
                         ncol=ncol(patient_preds))
  #rownames(patient_preds) = gsub('-', '.', rownames(epithelial_metadata)[epithelial_metadata$state=='Tumor' &
  #                                                                         epithelial_metadata$patient==patient])
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
  predictions_clusts_largeScale_allchrs[[patient]] = ComplexHeatmap::Heatmap(#patient_preds[rownames(predictions_clusts[[patient]]@matrix),],
                                                                             #cluster_rows=FALSE,
                                                                             patient_preds, cluster_rows=TRUE,
                                                                             cluster_columns=FALSE, col=heatmap_colors_discrete,
                                                                             name='Heatmap Values', show_row_names=FALSE,
                                                                             show_column_names=FALSE, column_title = patient,
                                                                             column_split=factor(colnames(cna_matrix), colnames(cna_matrix)),
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
#predictions_clusts_largeScale_relchrs = list()
#for(patient in grep('N', unique(epithelial_metadata$patient), value=TRUE, invert=TRUE)){
#  message(patient)
#  message('- Getting patient data')
#  patient_preds = as.matrix(cna_matrix[gsub('-', '.', rownames(epithelial_metadata)[epithelial_metadata$state=='Tumor' &
#                                                                                      epithelial_metadata$patient==patient]), relevant_chrs])
#  patient_preds[patient_preds=='Both'] = 3
#  patient_preds[patient_preds=='Gain'] = 2
#  patient_preds[patient_preds=='Neutral'] = 0
#  patient_preds[patient_preds=='Loss'] = 1
#  patient_preds = matrix(as.numeric(patient_preds),
#                         ncol=ncol(patient_preds))
#  colnames(patient_preds) = relevant_chrs
#  rownames(patient_preds) = gsub('-', '.', rownames(epithelial_metadata)[epithelial_metadata$state=='Tumor' &
#                                                                           epithelial_metadata$patient==patient])
#  cell_annotations = cell_annotation[gsub('-', '.', rownames(epithelial_metadata)[epithelial_metadata$state=='Tumor' &
#                                                                                    epithelial_metadata$patient==patient])]
#  message('- Performing heatmap')
#  col_annotation = ComplexHeatmap::HeatmapAnnotation(foo = ComplexHeatmap::anno_block(gp = grid::gpar(fill = 2:length(relevant_chrs)),
#                                                                                      labels=relevant_chrs,
#                                                                                      labels_gp=grid::gpar(#col=chrs_text_colors,
#                                                                                                           fontsize=7)),
#                                                     Literature = names(relevant_chrs),
#                                                     col= list(Literature=c('Gain'='red', 'Loss'='blue', 'Both'='grey')),
#                                                     show_annotation_name=FALSE)
#  row_annotation = ComplexHeatmap::rowAnnotation(Annotation_1 = cell_annotations,
#                                                 col = list(Annotation_1 = c('Normal'='grey', 'Putative Normal'='forestgreen',
#                                                                             'Putative Tumour'='chocolate2')),
#                                                 show_annotation_name=FALSE)
#  predictions_clusts_largeScale_relchrs[[patient]] = ComplexHeatmap::Heatmap(patient_preds, cluster_rows=TRUE,
#                                                                             cluster_columns=FALSE, col=heatmap_colors_discrete,
#                                                                             name='Heatmap Values', show_row_names=FALSE,
#                                                                             show_column_names=FALSE, column_title = patient,
#                                                                             column_split=factor(relevant_chrs, relevant_chrs),
#                                                                             top_annotation=col_annotation, left_annotation=row_annotation,
#                                                                             jitter=TRUE)
#  invisible(gc())
#  message('')
#}
#pdf(paste(project_dir, '2_annotation/CNV/largeScaleCNA_matrix_patient_relevant_chrs.pdf', sep='/'), width=20, height=12)
#for(patient in names(predictions_clusts_largeScale_relchrs)){
#  message(patient)
#  ComplexHeatmap::draw(predictions_clusts_largeScale_relchrs[[patient]])
#}
#dev.off()


# 4. Store CNA analysis information in metadata slot of Epithelial's seurat dataset:
Epithelial = SeuratDisk::LoadH5Seurat(paste(project_dir, '2_annotation/data/Epithelial_integratedClusters.h5Seurat', sep='/'))
invisible(gc())
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
invisible(gc())


# 1. Normal + Putative normal: Subset cells and run UMAP
# 1.1. Only those from normal samples and classified as putative normal
Epithelial_normal = subset(Epithelial, cells=rownames(Epithelial@meta.data)[Epithelial@meta.data$CNA_annotation!='Putative Tumour'])
Epithelial_normal = Seurat::RunPCA(Epithelial_normal, assay='integrated')
Epithelial_normal = Seurat::RunUMAP(Epithelial_normal, dims=1:50, reduction='pca')
invisible(gc())
# 1.2. Find clusters:
Epithelial_normal = Seurat::FindNeighbors(Epithelial_normal, dims=1:50)
Epithelial_normal = Seurat::FindClusters(Epithelial_normal, resolution=c(1))
invisible(gc())


# 2.  Normal + Putative normal: Check how cells cluster together
# 2.1. Normal and classified as normal:
normal_clusters = Seurat::DimPlot(Epithelial_normal, reduction="umap", group.by='integrated_snn_res.1', label=T,
                                  label.size=3, pt.size=.2) + ggplot2::theme_minimal() + Seurat::NoLegend()
normal_cna = Seurat::DimPlot(Epithelial_normal, reduction="umap", group.by='CNA_annotation', label=F,
                             label.size=3, pt.size=.2, order='Putative Normal', cols=c('Normal'='grey', 'Putative Normal'='red')) +
  ggplot2::theme_minimal() + ggplot2::ggtitle('Putative normal cells highlighted') + Seurat::NoLegend() 
normal_patients = Seurat::DimPlot(Epithelial_normal, reduction="umap", group.by='patient', label=T, label.size=3, pt.size=.2) +
  ggplot2::theme_minimal() + ggplot2::ggtitle('Colored by Patient') + Seurat::NoLegend()
(normal_cna | normal_patients)
normal_clusters | normal_cna | normal_patients
# 2.2. Check normal and putative normal distribution in clusters 21 and 29
table(Epithelial_normal@meta.data[Epithelial_normal$integrated_snn_res.1=='21','CNA_annotation'])
table(Epithelial_normal@meta.data[Epithelial_normal$integrated_snn_res.1=='29','CNA_annotation'])
# Out of the 903 cells classified as putative normal, 510 (~56%) were clustered into two clusters (21 and 29), while the other 393 (~44%)
# are distributed among other clusters.
# Cluster 21 still contains some normal cells (~21%, all of them being from healthy samples, no matched
# normal) clustered together with the putative normal ones. The healthy cells span 12 samples, while the putative normal span 5 patients.
# However, these healthy cells do not map together in the umap!
Seurat::DimPlot(subset(Epithelial_normal,
                       cells=rownames(Epithelial_normal@meta.data)[Epithelial_normal$integrated_snn_res.1%in%c('21', '29')]),
                reduction="umap", group.by='patient', split.by='integrated_snn_res.1', label=T, label.size=3, pt.size=.2) + Seurat::NoLegend()
# On the other hand, cluster 29 contains only putative normal cells, all from the same patient (128 cells from SMC10). However, not all cells
# classified as putative normal from this patient are in this cluster. Other 136 are seprated among other various clusters (including 21). 
# Could thus this cells, from cluster 29, be in fact tumour cells?
# 2.3. Where are putative normal mapped?
patients_mapped = list()
for(patient in c('SMC07', 'KUL01', 'SMC03', '31', 'SMC10', 'SMC19', 'SMC24')){
  patient_tumour_cells = rownames(Epithelial_normal@meta.data)[Epithelial_normal$patient==patient & Epithelial_normal$state=='Tumor']
  patients_mapped[[patient]] = Seurat::DimPlot(Epithelial_normal, cells.highlight = patient_tumour_cells) + 
    ggplot2::theme_minimal() + Seurat::NoLegend() + ggplot2::ggtitle(patient)
}
patchwork::wrap_plots(patients_mapped, ncol=3)

# 3.  Normal + Putative normal: Find Markers for clusters 21 and 29
markers_21 = Seurat::FindMarkers(Epithelial_normal, assay='integrated', ident.1='21', logfc.threshold=0.8)
markers_29 = Seurat::FindMarkers(Epithelial_normal, assay='integrated', ident.1='29', logfc.threshold=0.8)











# 3. Cells were classified into normal/tumour based on relevant chromosomes, lets check the other ones:
color_vec = c('Gain'='red', 'Loss'='blue', 'Neutral'='lightgrey')
relevant_chrs = c('chr1p', 'chr1q', 'chr4p', 'chr4q', 'chr5q', 'chr7p', 'chr7q', 'chr8p', 'chr8q', 'chr10q', 'chr13q', 'chr14q',
                  'chr15q', 'chr17p', 'chr18q', 'chr19p', 'chr19q', 'chr20q')
all_chrs = grep('CNA_chr', colnames(Epithelial_normal@meta.data), value=TRUE)
not_relevant_chrs = all_chrs[!all_chrs%in%paste('CNA', relevant_chrs, sep='_')]
Seurat::DimPlot(Epithelial_normal, reduction="umap", group.by=not_relevant_chrs[1:9],
                label=FALSE, pt.size=.2, cols=color_vec, ncol=3, order=c('Gain', 'Loss'))
Seurat::DimPlot(Epithelial_normal, reduction="umap", group.by=not_relevant_chrs[10:18],
                label=FALSE, pt.size=.2, cols=color_vec, ncol=3, order=c('Gain', 'Loss'))
Seurat::DimPlot(Epithelial_normal, reduction="umap", group.by=not_relevant_chrs[19:21],
                label=FALSE, pt.size=.2, cols=color_vec, ncol=3, order=c('Gain', 'Loss'))

















# 4. Visualy assess CNA analysis
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