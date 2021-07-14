project_dir = '/home/scardoso/Documents/PhD/CRC_ATLAS'

# ###################### #
# #### CNV ANALYSIS #### #
# ###################### #



# -----
# - Run Analysis
# -----

Epithelial = SeuratDisk::LoadH5Seurat(paste(project_dir, '2_annotation/results_Epithelial/datasets/Epithelial2.h5Seurat', sep='/'))
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
            paste(project_dir, '2_annotation/results_Epithelial/inferCNV/annotation.txt', sep='/'), sep='\t', col.names=FALSE, row.names=FALSE)

# Create inferCNV object:
infercnv_obj = infercnv::CreateInfercnvObject(epithelial_m,
                                              annotations_file=paste(project_dir, '2_annotation/results_Epithelial/inferCNV/annotation.txt', sep='/'),
                                              gene_order_file = paste(project_dir, 'utils/gencode_v19_gene_pos.txt', sep='/'),
                                              delim='\t', ref_group_names = c('Normal'))
invisible(gc())

# Run inferCNV:
infercnv_obj = infercnv::run(infercnv_obj, cutoff=0.1, out_dir=paste(project_dir, '2_annotation/results_Epithelial/inferCNV', sep='/'),
                             cluster_by_groups=T, cluster_references=F,
                             analysis_mode='subclusters', tumor_subcluster_partition_method = 'qnorm',
                             denoise=T, HMM=T, scale_data=T,
                             BayesMaxPNormal=0)



# ---
# - Organize data to use 
# ---

# Tranfer CNV predictions to seurat dataset
Epithelial = SeuratDisk::LoadH5Seurat(paste(project_dir, '2_annotation/results_Epithelial/inferCNV/Epithelial2.h5Seurat', sep='/'))
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
predictions_obs = read.table(paste(project_dir, '2_annotation/results_Epithelial/inferCNV/infercnv.17_HMM_predHMMi6.qnorm.hmm_mode-subclusters.observations.txt',
                                   sep='/'),
                             header=TRUE)
predictions_obs = as(as.matrix(predictions_obs), "sparseMatrix") 
predictions_refs = read.table(paste(project_dir, '2_annotation/results_Epithelial/inferCNV/infercnv.17_HMM_predHMMi6.qnorm.hmm_mode-subclusters.references.txt',
                                    sep='/'),
                              header=TRUE)
predictions_refs = as(as.matrix(predictions_refs), "sparseMatrix")

# Information on which chromosomes the genes belong to:
chrs_genes = read.table(paste(project_dir, 'utils/gene_pos_wArms.txt', sep='/'), row.names=1)





# ---
# - Tumour vs Normal classification based on CNAs
# ---


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


# 2. Classify patient cells (cells with at least one aberration in a chr arm will be considered as tumour)
cell_annotation = rep('Putative Normal', dim(cna_matrix)[1])
cell_annotation[rowSums(cna_matrix!='Neutral') > 0] = 'Putative Tumour'
names(cell_annotation) = rownames(cna_matrix)





# ---
# - Hierarchical clusterings of inferCNV results
# ---

all_chrs_text_colors = rep('black', ncol(cna_matrix))
all_chrs_text_colors[c(5, 8, 13, 16, 21, 24, 29, 32, 37)] = 'white'
heatmap_colors = circlize::colorRamp2(c(1, 3, 6), c("blue", "white", "red"))

# 1. Heatmap of predictions
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
# 1.1. Store plots
pdf(paste(project_dir, '2_annotation/results_Epithelial/inferCNV/predictions_heatmaps_patient.pdf', sep='/'), width=20, height=12)
for(patient in names(predictions_clusts)){
  message(patient)
  ComplexHeatmap::draw(predictions_clusts[[patient]])
}
dev.off()


# 2. Heatmaps of large-scale CNA
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
# 2.1. Store plots
pdf(paste(project_dir, '2_annotation/results_Epithelial/inferCNV/largeScaleCNA_matrix_patient_all_chrs.pdf', sep='/'), width=20, height=12)
for(patient in names(predictions_clusts_largeScale_allchrs)){
  message(patient)
  ComplexHeatmap::draw(predictions_clusts_largeScale_allchrs[[patient]])
}
dev.off()


# 3. Store CNA analysis information in metadata slot of Epithelial's seurat dataset:
Epithelial = SeuratDisk::LoadH5Seurat(paste(project_dir, '2_annotation/results_Epithelial/datasets/Epithelial2.h5Seurat', sep='/'))
invisible(gc())
# 3.1. CNA annotation
Epithelial[['CNA_annotation']] = rep('Normal', dim(Epithelial@meta.data)[1])
Epithelial$CNA_annotation[gsub('[.]', '-', names(cell_annotation))] = cell_annotation
# 3.2. CNA matrix data
for(chr in colnames(cna_matrix)){
  meta_name = paste('CNA', chr, sep='_')
  Epithelial[[meta_name]] = rep('Neutral', dim(Epithelial@meta.data)[1])
  Epithelial@meta.data[gsub('[.]', '-', rownames(cna_matrix)), meta_name] = cna_matrix[, chr]
}
# 3.3. Save seurat object with new information
SeuratDisk::SaveH5Seurat(Epithelial, paste(project_dir, '2_annotation/results_Epithelial/datasets/Epithelial3.h5Seurat', sep='/'))




