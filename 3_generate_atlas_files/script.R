project_dir = '/home/scardoso/Documents/PhD/CRC_ATLAS'





# -----
# - Organize metadata
# -----

# 1. Normal Epithelial
# 1.1. Load dataset
Epithelial_normal = SeuratDisk::LoadH5Seurat(paste(project_dir,
                                                   '2_annotation/results_Epithelial/datasets/Epithelial_normal_finalAnnots.h5Seurat', sep='/'))
# 1.2. Extract metadata
Epithelial_normal_meta = Epithelial_normal@meta.data
remove(Epithelial_normal)
invisible(gc())
# 1.3. Organize metadata
Epithelial_normal_meta = Epithelial_normal_meta[, c(4:8, 10, 14:21, 23:61, 64, 65)]
colnames(Epithelial_normal_meta) = c('Sample', 'nUMI', 'nGene', 'Individual', 'Sample.Source', 'Dataset', 'Percent.Mito.RNA',
                                     'log10.Complexity', 'S.Score', 'G2M.Score', 'Phase', 'CC.Difference', 'Sample.State', 'Annotation_0',
                                     gsub('CNA_', 'CNA.', grep('CNA_chr', colnames(Epithelial_normal_meta), value=T)),
                                     'Annotation_2', 'Annotation_1')
Epithelial_normal_meta$Sample.Source[Epithelial_normal_meta$Sample.Source=='Tumor'] = 'Tumour'
Epithelial_normal_meta$Sample.Source[Epithelial_normal_meta$Sample.Source%in%c('Normal Tissue', 'Normal')] = 'Normal Matched'
Epithelial_normal_meta$Sample.Source[Epithelial_normal_meta$Sample.Source=='Healthy'] = 'Healthy Donors'
Epithelial_normal_meta = Epithelial_normal_meta[, c(2,3,7,8,9,10,12,11,6,4,1,5,13,14,55,54,15:53)]
# 1.4. Save metadata
write.csv(Epithelial_normal_meta, paste(project_dir, 'FINAL_ATLAS/metadata/Epithelial_normal.csv', sep='/'))


# 2. Tumour Epithelial
# 2.1. Load dataset
Epithelial_tumour = SeuratDisk::LoadH5Seurat(paste(project_dir,
                                                   '2_annotation/results_Epithelial/datasets/Epithelial_tumour_finalAnnots.h5Seurat', sep='/'))
# 2.2. Extract metadata
Epithelial_tumour_meta = Epithelial_tumour@meta.data
remove(Epithelial_tumour)
invisible(gc())
# 2.3. Organize metadata
Epithelial_tumour_meta = Epithelial_tumour_meta[, c(4:8, 10, 14:21, 23:61, 65)]
colnames(Epithelial_tumour_meta) = c('Sample', 'nUMI', 'nGene', 'Individual', 'Sample.Source', 'Dataset', 'Percent.Mito.RNA',
                                     'log10.Complexity', 'S.Score', 'G2M.Score', 'Phase', 'CC.Difference', 'Sample.State', 'Annotation_0',
                                     gsub('CNA_', 'CNA.', grep('CNA_chr', colnames(Epithelial_tumour_meta), value=T)),
                                     'Annotation_2')
Epithelial_tumour_meta$Sample.Source[Epithelial_tumour_meta$Sample.Source=='Tumor'] = 'Tumour'
Epithelial_tumour_meta$Sample.Source[Epithelial_tumour_meta$Sample.Source%in%c('Normal Tissue', 'Normal')] = 'Normal Matched'
Epithelial_tumour_meta$Sample.Source[Epithelial_tumour_meta$Sample.Source=='Healthy'] = 'Healthy Donors'
Epithelial_tumour_meta$Annotation_1 = rep('Tumour Epithelial cells', dim(Epithelial_tumour_meta)[1])
Epithelial_tumour_meta = Epithelial_tumour_meta[, c(2,3,7,8,9,10,12,11,6,4,1,5,13,14,55,54,15:53)]
# 2.4. Save metadata
write.csv(Epithelial_tumour_meta, paste(project_dir, 'FINAL_ATLAS/metadata/Epithelial_tumour.csv', sep='/'))


# 3. Stromal cells
# 3.1. Load dataset
Stromal = SeuratDisk::LoadH5Seurat(paste(project_dir, '2_annotation/results_Stromal/datasets/Stromal_finalAnnots.h5Seurat', sep='/'))
# 3.2. Extract metadata
Stromal_meta = Stromal@meta.data
remove(Stromal)
invisible(gc())
# 3.3. Organize metadata
Stromal_meta = Stromal_meta[, c(4:8, 10, 14:24)]
colnames(Stromal_meta) = c('Sample', 'nUMI', 'nGene', 'Individual', 'Sample.Source', 'Dataset', 'Percent.Mito.RNA', 'log10.Complexity',
                           'S.Score', 'G2M.Score', 'Phase', 'CC.Difference', 'Sample.State', 'Annotation_0', 'Annotation_3', 'Annotation_2',
                           'Annotation_1')
Stromal_meta$Sample.Source[Stromal_meta$Sample.Source=='Tumor'] = 'Tumour'
Stromal_meta$Sample.Source[Stromal_meta$Sample.Source%in%c('Normal Tissue', 'Normal')] = 'Normal Matched'
Stromal_meta$Sample.Source[Stromal_meta$Sample.Source=='Healthy'] = 'Healthy Donors'
Stromal_meta = Stromal_meta[, c(2,3,7,8,9,10,12,11,6,4,1,5,13,14,17,16,15)]
# 3.4. Save metadata
write.csv(Stromal_meta, paste(project_dir, 'FINAL_ATLAS/metadata/Stromal_cells.csv', sep='/'))


# 4. Myeloid cells
# 4.1. Load dataset
Myeloid = SeuratDisk::LoadH5Seurat(paste(project_dir, '2_annotation/results_Myeloid/datasets/Myeloid_finalAnnots.h5Seurat', sep='/'))
# 4.2. Extract metadata
Myeloid_meta = Myeloid@meta.data
remove(Myeloid)
invisible(gc())
# 4.3. Organize metadata
Myeloid_meta = Myeloid_meta[, c(4:8, 10, 14:23)]
colnames(Myeloid_meta) = c('Sample', 'nUMI', 'nGene', 'Individual', 'Sample.Source', 'Dataset', 'Percent.Mito.RNA', 'log10.Complexity',
                           'S.Score', 'G2M.Score', 'Phase', 'CC.Difference', 'Sample.State', 'Annotation_0', 'Annotation_2', 'Annotation_1')
Myeloid_meta$Sample.Source[Myeloid_meta$Sample.Source=='Tumor'] = 'Tumour'
Myeloid_meta$Sample.Source[Myeloid_meta$Sample.Source%in%c('Normal Tissue', 'Normal')] = 'Normal Matched'
Myeloid_meta$Sample.Source[Myeloid_meta$Sample.Source=='Healthy'] = 'Healthy Donors'
Myeloid_meta = Myeloid_meta[, c(2,3,7,8,9,10,12,11,6,4,1,5,13,14,16,15)]
# 4.4. Save metadata
write.csv(Myeloid_meta, paste(project_dir, 'FINAL_ATLAS/metadata/Myeloid_cells.csv', sep='/'))


# 5. B cells
# 5.1. Load dataset
Bcells = SeuratDisk::LoadH5Seurat(paste(project_dir, '2_annotation/results_Bcells/datasets/Bcells_finalAnnots.h5Seurat', sep='/'))
# 5.2. Extract metadata
Bcells_meta = Bcells@meta.data
remove(Bcells)
invisible(gc())
# 5.3. Organize metadata
Bcells_meta = Bcells_meta[, c(4:8, 10, 14:21, 25, 26)]
colnames(Bcells_meta) = c('Sample', 'nUMI', 'nGene', 'Individual', 'Sample.Source', 'Dataset', 'Percent.Mito.RNA', 'log10.Complexity',
                          'S.Score', 'G2M.Score', 'Phase', 'CC.Difference', 'Sample.State', 'Annotation_0', 'Annotation_2', 'Annotation_1')
Bcells_meta$Sample.Source[Bcells_meta$Sample.Source=='Tumor'] = 'Tumour'
Bcells_meta$Sample.Source[Bcells_meta$Sample.Source%in%c('Normal Tissue', 'Normal')] = 'Normal Matched'
Bcells_meta$Sample.Source[Bcells_meta$Sample.Source=='Healthy'] = 'Healthy Donors'
Bcells_meta = Bcells_meta[, c(2,3,7,8,9,10,12,11,6,4,1,5,13,14,16,15)]
# 5.4. Save metadata
write.csv(Bcells_meta, paste(project_dir, 'FINAL_ATLAS/metadata/B_cells.csv', sep='/'))


# 6. T cells
# 6.1. Load dataset
Tcells = SeuratDisk::LoadH5Seurat(paste(project_dir, '2_annotation/results_Tcells/datasets/Tcells_finalAnnots.h5Seurat', sep='/'))
# 6.2. Extract metadata
Tcells_meta = Tcells@meta.data
remove(Tcells)
invisible(gc())
# 6.3. Organize metadata
Tcells_meta = Tcells_meta[, c(4:8, 10, 14:21, 25:30)]
colnames(Tcells_meta) = c('Sample', 'nUMI', 'nGene', 'Individual', 'Sample.Source', 'Dataset', 'Percent.Mito.RNA', 'log10.Complexity',
                          'S.Score', 'G2M.Score', 'Phase', 'CC.Difference', 'Sample.State', 'Annotation_0', 'Annotation_6', 'Annotation_5',
                          'Annotation_4', 'Annotation_3', 'Annotation_2', 'Annotation_1')
Tcells_meta$Sample.Source[Tcells_meta$Sample.Source=='Tumor'] = 'Tumour'
Tcells_meta$Sample.Source[Tcells_meta$Sample.Source%in%c('Normal Tissue', 'Normal')] = 'Normal Matched'
Tcells_meta$Sample.Source[Tcells_meta$Sample.Source=='Healthy'] = 'Healthy Donors'
Tcells_meta = Tcells_meta[, c(2,3,7,8,9,10,12,11,6,4,1,5,13,14,20,19,18,17,16,15)]
# 6.4. Save metadata
write.csv(Tcells_meta, paste(project_dir, 'FINAL_ATLAS/metadata/T_cells.csv', sep='/'))


# 7. Join metadata from the different groups together to make the global metadata:
# 7.1. Bind first the columns that are common to all metadatas
all_meta = rbind(Epithelial_normal_meta[,1:15], Epithelial_tumour_meta[,1:15], Stromal_meta[,1:15],
                 Myeloid_meta[,1:15], Bcells_meta[,1:15], Tcells_meta[,1:15])
# 7.2. Add annotation level 2
all_meta$Annotation_2 = rep(NA, dim(all_meta)[1])
all_meta[rownames(Epithelial_normal_meta), 'Annotation_2'] = as.character(Epithelial_normal_meta$Annotation_1)
all_meta[rownames(Epithelial_tumour_meta), 'Annotation_2'] = as.character(Epithelial_tumour_meta$Annotation_1)
all_meta[rownames(Stromal_meta), 'Annotation_2'] = as.character(Stromal_meta$Annotation_1)
all_meta[rownames(Myeloid_meta), 'Annotation_2'] = as.character(Myeloid_meta$Annotation_1)
all_meta[rownames(Bcells_meta), 'Annotation_2'] = as.character(Bcells_meta$Annotation_1)
all_meta[rownames(Tcells_meta), 'Annotation_2'] = as.character(Tcells_meta$Annotation_2)
# 7.3. Add annotation level 3
all_meta$Annotation_3 = rep(NA, dim(all_meta)[1])
all_meta[rownames(Epithelial_normal_meta), 'Annotation_3'] = as.character(Epithelial_normal_meta$Annotation_2)
all_meta[rownames(Epithelial_tumour_meta), 'Annotation_3'] = as.character(Epithelial_tumour_meta$Annotation_1)
all_meta[rownames(Stromal_meta), 'Annotation_3'] = as.character(Stromal_meta$Annotation_1)
all_meta[rownames(Myeloid_meta), 'Annotation_3'] = as.character(Myeloid_meta$Annotation_1)
all_meta[rownames(Bcells_meta), 'Annotation_3'] = as.character(Bcells_meta$Annotation_1)
all_meta[rownames(Tcells_meta), 'Annotation_3'] = as.character(Tcells_meta$Annotation_3)
# 7.4. Add annotation level 4
all_meta$Annotation_4 = all_meta$Annotation_3
all_meta[rownames(Tcells_meta), 'Annotation_4'] = as.character(Tcells_meta$Annotation_4)
# 7.5. Add annotation level 5
all_meta$Annotation_5 = rep(NA, dim(all_meta)[1])
all_meta[rownames(Epithelial_normal_meta), 'Annotation_5'] = as.character(Epithelial_normal_meta$Annotation_2)
all_meta[rownames(Epithelial_tumour_meta), 'Annotation_5'] = as.character(Epithelial_tumour_meta$Annotation_2)
all_meta[rownames(Stromal_meta), 'Annotation_5'] = as.character(Stromal_meta$Annotation_2)
all_meta[rownames(Myeloid_meta), 'Annotation_5'] = as.character(Myeloid_meta$Annotation_2)
all_meta[rownames(Bcells_meta), 'Annotation_5'] = as.character(Bcells_meta$Annotation_1)
all_meta[rownames(Tcells_meta), 'Annotation_5'] = as.character(Tcells_meta$Annotation_5)
# 7.5. Add annotation level 6
all_meta$Annotation_6 = rep(NA, dim(all_meta)[1])
all_meta[rownames(Epithelial_normal_meta), 'Annotation_6'] = as.character(Epithelial_normal_meta$Annotation_2)
all_meta[rownames(Epithelial_tumour_meta), 'Annotation_6'] = as.character(Epithelial_tumour_meta$Annotation_2)
all_meta[rownames(Stromal_meta), 'Annotation_6'] = as.character(Stromal_meta$Annotation_3)
all_meta[rownames(Myeloid_meta), 'Annotation_6'] = as.character(Myeloid_meta$Annotation_2)
all_meta[rownames(Bcells_meta), 'Annotation_6'] = as.character(Bcells_meta$Annotation_2)
all_meta[rownames(Tcells_meta), 'Annotation_6'] = as.character(Tcells_meta$Annotation_6)
# 7.6. Add CNAs
n_notEpithelial = dim(all_meta)[1] - (dim(Epithelial_normal_meta)[1] + dim(Epithelial_tumour_meta)[1])
cnas_notEpithelial = matrix(rep('Neutral', n_notEpithelial*39), ncol=39)
dimnames(cnas_notEpithelial) = list(rownames(all_meta)[all_meta$Annotation_0!='Epithelial cells'], colnames(Epithelial_normal_meta)[17:55])
cnas_df = rbind(Epithelial_normal_meta[,17:55], Epithelial_tumour_meta[,17:55], cnas_notEpithelial)
all_meta = cbind(all_meta[rownames(cnas_df),], cnas_df)
# 7.7. Save metadata
write.csv(all_meta, paste(project_dir, 'FINAL_ATLAS/metadata/all_cells.csv', sep='/'))
# 7.8. Remove objects that are not necessary anymore
remove(Bcells_meta, cnas_df, cnas_notEpithelial, Epithelial_normal_meta, Epithelial_tumour_meta, Myeloid_meta, Stromal_meta, Tcells_meta,
       n_notEpithelial)
invisible(gc())




# -----
# - Create Full Atlas (Patients + Healthy Donors)
# -----
library(SeuratObject)


# 1. Load full dataset
CRCatlas_PatientsHealthyDonors = SeuratDisk::LoadH5Seurat(
  paste(project_dir, '2_annotation/results_globalAnnotation/datasets/CRC_annotations.h5Seurat', sep='/'))
invisible(gc())


# 2. Remove cells that were removed while performing annotations
dim(CRCatlas_PatientsHealthyDonors@meta.data)
CRCatlas_PatientsHealthyDonors = subset(CRCatlas_PatientsHealthyDonors, cells=rownames(all_meta))
invisible(gc())


# 3. Update metadata
CRCatlas_PatientsHealthyDonors@meta.data = all_meta[rownames(CRCatlas_PatientsHealthyDonors[[]]), ]


# 4. Create files for all cells:
# 4.1. Seurat, only for further use when performing tumour deconvolution (integrated assay necessary): version 3.1.5.9900
SeuratDisk::SaveH5Seurat(CRCatlas_PatientsHealthyDonors,
                         paste(project_dir, 'FINAL_ATLAS/internal/CRCatlas.h5Seurat', sep='/'))
invisible(gc())
# 4.2. Seurat only with RNA assay: version 3.1.5.9900
Seurat::DefaultAssay(CRCatlas_PatientsHealthyDonors) = 'RNA'
CRCatlas_PatientsHealthyDonors_RNAassay = Seurat::DietSeurat(CRCatlas_PatientsHealthyDonors, assays='RNA')
remove(CRCatlas_PatientsHealthyDonors)
invisible(gc())
SeuratDisk::SaveH5Seurat(CRCatlas_PatientsHealthyDonors_RNAassay,
                         paste(project_dir, 'FINAL_ATLAS/0_Patients_HealthyDonors/All_cells/CRCatlas_PH.h5Seurat', sep='/'))
invisible(gc())

# 4.3. H5AD
SeuratDisk::Convert(paste(project_dir, 'FINAL_ATLAS/0_Patients_HealthyDonors/All_cells/CRCatlas_PH.h5Seurat', sep='/'),
                    dest='h5ad', assay='RNA')
invisible(gc())

# 4.4. MTX
PH_matrix = SeuratObject::GetAssayData(CRCatlas_PatientsHealthyDonors_RNAassay, slot='counts', assay='RNA')
dir.create(paste(project_dir, 'FINAL_ATLAS/0_Patients_HealthyDonors/All_cells/mtx', sep='/'))
Matrix::writeMM(obj=PH_matrix,
                file=paste(project_dir, 'FINAL_ATLAS/0_Patients_HealthyDonors/All_cells/mtx/matrix.mtx', sep='/'))
write(x=rownames(PH_matrix), file=paste(project_dir, 'FINAL_ATLAS/0_Patients_HealthyDonors/All_cells/mtx/genes.tsv', sep='/'))
write(x=colnames(PH_matrix), file=paste(project_dir, 'FINAL_ATLAS/0_Patients_HealthyDonors/All_cells/mtx/barcodes.tsv', sep='/'))
invisible(gc())

# 4.5. Metadata CSV
write.csv(CRCatlas_PatientsHealthyDonors_RNAassay@meta.data,
          paste(project_dir, 'FINAL_ATLAS/0_Patients_HealthyDonors/All_cells/metadata.csv', sep='/'),
          row.names=TRUE)


# 5. Create files for epithelial cells:
CRCatlas_PH_epithelial = subset(CRCatlas_PatientsHealthyDonors_RNAassay,
                                cells=rownames(CRCatlas_PatientsHealthyDonors_RNAassay[[]])[CRCatlas_PatientsHealthyDonors_RNAassay$Annotation_0=='Epithelial cells'])
invisible(gc())
# 5.1. Seurat only with RNA assay: version 3.1.5.9900
SeuratDisk::SaveH5Seurat(CRCatlas_PH_epithelial,
                         paste(project_dir, 'FINAL_ATLAS/0_Patients_HealthyDonors/Epithelial_cells/CRCatlas_PH_epithelial.h5Seurat', sep='/'))
invisible(gc())

# 5.2. H5AD
SeuratDisk::Convert(paste(project_dir, 'FINAL_ATLAS/0_Patients_HealthyDonors/Epithelial_cells/CRCatlas_PH_epithelial.h5Seurat', sep='/'),
                    dest='h5ad', assay='RNA')
invisible(gc())

# 5.3. MTX
PH_matrix = SeuratObject::GetAssayData(CRCatlas_PH_epithelial, slot='counts', assay='RNA')
dir.create(paste(project_dir, 'FINAL_ATLAS/0_Patients_HealthyDonors/Epithelial_cells/mtx', sep='/'))
Matrix::writeMM(obj=PH_matrix,
                file=paste(project_dir, 'FINAL_ATLAS/0_Patients_HealthyDonors/Epithelial_cells/mtx/matrix.mtx', sep='/'))
write(x=rownames(PH_matrix), file=paste(project_dir, 'FINAL_ATLAS/0_Patients_HealthyDonors/Epithelial_cells/mtx/genes.tsv', sep='/'))
write(x=colnames(PH_matrix), file=paste(project_dir, 'FINAL_ATLAS/0_Patients_HealthyDonors/Epithelial_cells/mtx/barcodes.tsv', sep='/'))
invisible(gc())

# 5.4. Metadata CSV
write.csv(CRCatlas_PH_epithelial@meta.data,
          paste(project_dir, 'FINAL_ATLAS/0_Patients_HealthyDonors/Epithelial_cells/metadata.csv', sep='/'),
          row.names=TRUE)

# 5.5. Remove epithelial object
remove(CRCatlas_PH_epithelial, PH_matrix)
invisible(gc())


# 6. Create files for Stromal cells:
CRCatlas_PH_stromal = subset(CRCatlas_PatientsHealthyDonors_RNAassay,
                             cells=rownames(CRCatlas_PatientsHealthyDonors_RNAassay[[]])[CRCatlas_PatientsHealthyDonors_RNAassay$Annotation_0=='Stromal cells'])
invisible(gc())
# 6.1. Seurat only with RNA assay: version 3.1.5.9900
SeuratDisk::SaveH5Seurat(CRCatlas_PH_stromal,
                         paste(project_dir, 'FINAL_ATLAS/0_Patients_HealthyDonors/Stromal_cells/CRCatlas_PH_stromal.h5Seurat', sep='/'))
invisible(gc())

# 6.2. H5AD
SeuratDisk::Convert(paste(project_dir, 'FINAL_ATLAS/0_Patients_HealthyDonors/Stromal_cells/CRCatlas_PH_stromal.h5Seurat', sep='/'),
                    dest='h5ad', assay='RNA')
invisible(gc())

# 6.3. MTX
PH_matrix = SeuratObject::GetAssayData(CRCatlas_PH_stromal, slot='counts', assay='RNA')
dir.create(paste(project_dir, 'FINAL_ATLAS/0_Patients_HealthyDonors/Stromal_cells/mtx', sep='/'))
Matrix::writeMM(obj=PH_matrix,
                file=paste(project_dir, 'FINAL_ATLAS/0_Patients_HealthyDonors/Stromal_cells/mtx/matrix.mtx', sep='/'))
write(x=rownames(PH_matrix), file=paste(project_dir, 'FINAL_ATLAS/0_Patients_HealthyDonors/Stromal_cells/mtx/genes.tsv', sep='/'))
write(x=colnames(PH_matrix), file=paste(project_dir, 'FINAL_ATLAS/0_Patients_HealthyDonors/Stromal_cells/mtx/barcodes.tsv', sep='/'))
invisible(gc())

# 6.4. Metadata CSV
write.csv(CRCatlas_PH_stromal@meta.data,
          paste(project_dir, 'FINAL_ATLAS/0_Patients_HealthyDonors/Stromal_cells/metadata.csv', sep='/'),
          row.names=TRUE)

# 6.5. Remove epithelial object
remove(CRCatlas_PH_stromal, PH_matrix)
invisible(gc())


# 7. Create files for Myeloid cells:
CRCatlas_PH_myeloid = subset(CRCatlas_PatientsHealthyDonors_RNAassay,
                                cells=rownames(CRCatlas_PatientsHealthyDonors_RNAassay[[]])[CRCatlas_PatientsHealthyDonors_RNAassay$Annotation_0=='Myeloid cells'])
invisible(gc())
# 7.1. Seurat only with RNA assay: version 3.1.5.9900
SeuratDisk::SaveH5Seurat(CRCatlas_PH_myeloid,
                         paste(project_dir, 'FINAL_ATLAS/0_Patients_HealthyDonors/Myeloid_cells/CRCatlas_PH_myeloid.h5Seurat', sep='/'))
invisible(gc())

# 7.2. H5AD
SeuratDisk::Convert(paste(project_dir, 'FINAL_ATLAS/0_Patients_HealthyDonors/Myeloid_cells/CRCatlas_PH_myeloid.h5Seurat', sep='/'),
                    dest='h5ad', assay='RNA')
invisible(gc())

# 7.3. MTX
PH_matrix = SeuratObject::GetAssayData(CRCatlas_PH_myeloid, slot='counts', assay='RNA')
dir.create(paste(project_dir, 'FINAL_ATLAS/0_Patients_HealthyDonors/Myeloid_cells/mtx', sep='/'))
Matrix::writeMM(obj=PH_matrix,
                file=paste(project_dir, 'FINAL_ATLAS/0_Patients_HealthyDonors/Myeloid_cells/mtx/matrix.mtx', sep='/'))
write(x=rownames(PH_matrix), file=paste(project_dir, 'FINAL_ATLAS/0_Patients_HealthyDonors/Myeloid_cells/mtx/genes.tsv', sep='/'))
write(x=colnames(PH_matrix), file=paste(project_dir, 'FINAL_ATLAS/0_Patients_HealthyDonors/Myeloid_cells/mtx/barcodes.tsv', sep='/'))
invisible(gc())

# 7.4. Metadata CSV
write.csv(CRCatlas_PH_myeloid@meta.data,
          paste(project_dir, 'FINAL_ATLAS/0_Patients_HealthyDonors/Myeloid_cells/metadata.csv', sep='/'),
          row.names=TRUE)

# 7.5. Remove epithelial object
remove(CRCatlas_PH_myeloid, PH_matrix)
invisible(gc())


# 8. Create files for B cells:
CRCatlas_PH_Bcells = subset(CRCatlas_PatientsHealthyDonors_RNAassay,
                             cells=rownames(CRCatlas_PatientsHealthyDonors_RNAassay[[]])[CRCatlas_PatientsHealthyDonors_RNAassay$Annotation_0=='Bcells'])
invisible(gc())
# 8.1. Seurat only with RNA assay: version 3.1.5.9900
SeuratDisk::SaveH5Seurat(CRCatlas_PH_Bcells,
                         paste(project_dir, 'FINAL_ATLAS/0_Patients_HealthyDonors/B_cells/CRCatlas_PH_Bcells.h5Seurat', sep='/'))
invisible(gc())

# 8.2. H5AD
SeuratDisk::Convert(paste(project_dir, 'FINAL_ATLAS/0_Patients_HealthyDonors/B_cells/CRCatlas_PH_Bcells.h5Seurat', sep='/'),
                    dest='h5ad', assay='RNA')
invisible(gc())

# 8.3. MTX
PH_matrix = SeuratObject::GetAssayData(CRCatlas_PH_Bcells, slot='counts', assay='RNA')
dir.create(paste(project_dir, 'FINAL_ATLAS/0_Patients_HealthyDonors/B_cells/mtx', sep='/'))
Matrix::writeMM(obj=PH_matrix,
                file=paste(project_dir, 'FINAL_ATLAS/0_Patients_HealthyDonors/B_cells/mtx/matrix.mtx', sep='/'))
write(x=rownames(PH_matrix), file=paste(project_dir, 'FINAL_ATLAS/0_Patients_HealthyDonors/B_cells/mtx/genes.tsv', sep='/'))
write(x=colnames(PH_matrix), file=paste(project_dir, 'FINAL_ATLAS/0_Patients_HealthyDonors/B_cells/mtx/barcodes.tsv', sep='/'))
invisible(gc())

# 8.4. Metadata CSV
write.csv(CRCatlas_PH_Bcells@meta.data,
          paste(project_dir, 'FINAL_ATLAS/0_Patients_HealthyDonors/B_cells/metadata.csv', sep='/'),
          row.names=TRUE)

# 8.5. Remove epithelial object
remove(CRCatlas_PH_Bcells, PH_matrix)
invisible(gc())


# 9. Create files for T cells:
CRCatlas_PH_Tcells = subset(CRCatlas_PatientsHealthyDonors_RNAassay,
                            cells=rownames(CRCatlas_PatientsHealthyDonors_RNAassay[[]])[CRCatlas_PatientsHealthyDonors_RNAassay$Annotation_0=='Tcells'])
invisible(gc())
# 9.1. Seurat only with RNA assay: version 3.1.5.9900
SeuratDisk::SaveH5Seurat(CRCatlas_PH_Tcells,
                         paste(project_dir, 'FINAL_ATLAS/0_Patients_HealthyDonors/T_cells/CRCatlas_PH_Tcells.h5Seurat', sep='/'))
invisible(gc())

# 9.2. H5AD
SeuratDisk::Convert(paste(project_dir, 'FINAL_ATLAS/0_Patients_HealthyDonors/T_cells/CRCatlas_PH_Tcells.h5Seurat', sep='/'),
                    dest='h5ad', assay='RNA')
invisible(gc())

# 9.3. MTX
PH_matrix = SeuratObject::GetAssayData(CRCatlas_PH_Tcells, slot='counts', assay='RNA')
dir.create(paste(project_dir, 'FINAL_ATLAS/0_Patients_HealthyDonors/T_cells/mtx', sep='/'))
Matrix::writeMM(obj=PH_matrix,
                file=paste(project_dir, 'FINAL_ATLAS/0_Patients_HealthyDonors/T_cells/mtx/matrix.mtx', sep='/'))
write(x=rownames(PH_matrix), file=paste(project_dir, 'FINAL_ATLAS/0_Patients_HealthyDonors/T_cells/mtx/genes.tsv', sep='/'))
write(x=colnames(PH_matrix), file=paste(project_dir, 'FINAL_ATLAS/0_Patients_HealthyDonors/T_cells/mtx/barcodes.tsv', sep='/'))
invisible(gc())

# 9.4. Metadata CSV
write.csv(CRCatlas_PH_Tcells@meta.data,
          paste(project_dir, 'FINAL_ATLAS/0_Patients_HealthyDonors/T_cells/metadata.csv', sep='/'),
          row.names=TRUE)

# 9.5. Remove epithelial object
remove(CRCatlas_PH_Tcells, PH_matrix)
invisible(gc())





# -----
# - Create Full Atlas (Patients: Tumour + Normal samples)
# -----
library(SeuratObject)


# 1. Load full dataset
CRCatlas_PH = SeuratDisk::LoadH5Seurat(paste(project_dir, 'FINAL_ATLAS/0_Patients_HealthyDonors/All_cells/CRCatlas_PH.h5Seurat', sep='/'))

# 2. Get only patients that have normal matched samples
pt_out = c()
for(pt in unique(CRCatlas_P_TN$Individual)) if(length(table(CRCatlas_P_TN$Sample.State[CRCatlas_P_TN$Individual==pt]))==1) pt_out=c(pt_out,pt)
CRCatlas_P_TN = subset(CRCatlas_PH, cells = rownames(CRCatlas_PH[[]])[CRCatlas_PH$Sample.State!='Healthy' & !CRCatlas_PH$Individual%in%pt_out])
remove(CRCatlas_PH)
invisible(gc())

# 3. Create files for all cells:
# 3.1. Seurat: version 3.1.5.9900
SeuratDisk::SaveH5Seurat(CRCatlas_P_TN,
                         paste(project_dir, 'FINAL_ATLAS/1_Patients_Tumour_NormalMatched/All_cells/CRCatlas_PTN.h5Seurat', sep='/'))
invisible(gc())

# 3.2. H5AD
SeuratDisk::Convert(paste(project_dir, 'FINAL_ATLAS/1_Patients_Tumour_NormalMatched/All_cells/CRCatlas_PTN.h5Seurat', sep='/'),
                    dest='h5ad', assay='RNA')
invisible(gc())

# 3.3. MTX
PTN_matrix = SeuratObject::GetAssayData(CRCatlas_P_TN, slot='counts', assay='RNA')
dir.create(paste(project_dir, 'FINAL_ATLAS/1_Patients_Tumour_NormalMatched/All_cells/mtx', sep='/'))
Matrix::writeMM(obj=PTN_matrix,
                file=paste(project_dir, 'FINAL_ATLAS/1_Patients_Tumour_NormalMatched/All_cells/mtx/matrix.mtx', sep='/'))
write(x=rownames(PTN_matrix), file=paste(project_dir, 'FINAL_ATLAS/1_Patients_Tumour_NormalMatched/All_cells/mtx/genes.tsv', sep='/'))
write(x=colnames(PTN_matrix), file=paste(project_dir, 'FINAL_ATLAS/1_Patients_Tumour_NormalMatched/All_cells/mtx/barcodes.tsv', sep='/'))
invisible(gc())

# 3.4. Metadata CSV
write.csv(CRCatlas_P_TN@meta.data,
          paste(project_dir, 'FINAL_ATLAS/1_Patients_Tumour_NormalMatched/All_cells/metadata.csv', sep='/'),
          row.names=TRUE)


# 4. Create files for epithelial cells:
CRCatlas_PTN_epithelial = subset(CRCatlas_P_TN,
                                cells=rownames(CRCatlas_P_TN[[]])[CRCatlas_P_TN$Annotation_0=='Epithelial cells'])
invisible(gc())
# 4.1. Seurat only with RNA assay: version 3.1.5.9900
SeuratDisk::SaveH5Seurat(CRCatlas_PTN_epithelial,
                         paste(project_dir, 'FINAL_ATLAS/1_Patients_Tumour_NormalMatched/Epithelial_cells/CRCatlas_PTN_epithelial.h5Seurat', sep='/'))
invisible(gc())

# 4.2. H5AD
SeuratDisk::Convert(paste(project_dir, 'FINAL_ATLAS/1_Patients_Tumour_NormalMatched/Epithelial_cells/CRCatlas_PTN_epithelial.h5Seurat', sep='/'),
                    dest='h5ad', assay='RNA')
invisible(gc())

# 4.3. MTX
PTN_matrix = SeuratObject::GetAssayData(CRCatlas_PTN_epithelial, slot='counts', assay='RNA')
dir.create(paste(project_dir, 'FINAL_ATLAS/1_Patients_Tumour_NormalMatched/Epithelial_cells/mtx', sep='/'))
Matrix::writeMM(obj=PTN_matrix,
                file=paste(project_dir, 'FINAL_ATLAS/1_Patients_Tumour_NormalMatched/Epithelial_cells/mtx/matrix.mtx', sep='/'))
write(x=rownames(PTN_matrix), file=paste(project_dir, 'FINAL_ATLAS/1_Patients_Tumour_NormalMatched/Epithelial_cells/mtx/genes.tsv', sep='/'))
write(x=colnames(PTN_matrix), file=paste(project_dir, 'FINAL_ATLAS/1_Patients_Tumour_NormalMatched/Epithelial_cells/mtx/barcodes.tsv', sep='/'))
invisible(gc())

# 4.4. Metadata CSV
write.csv(CRCatlas_PTN_epithelial@meta.data,
          paste(project_dir, 'FINAL_ATLAS/1_Patients_Tumour_NormalMatched/Epithelial_cells/metadata.csv', sep='/'),
          row.names=TRUE)

# 4.5. Remove epithelial object
remove(CRCatlas_PTN_epithelial, PTN_matrix)
invisible(gc())


# 5. Create files for Stromal cells:
CRCatlas_PTN_stromal = subset(CRCatlas_P_TN,
                             cells=rownames(CRCatlas_P_TN[[]])[CRCatlas_P_TN$Annotation_0=='Stromal cells'])
invisible(gc())
# 5.1. Seurat only with RNA assay: version 3.1.5.9900
SeuratDisk::SaveH5Seurat(CRCatlas_PTN_stromal,
                         paste(project_dir, 'FINAL_ATLAS/1_Patients_Tumour_NormalMatched/Stromal_cells/CRCatlas_PTN_stromal.h5Seurat', sep='/'))
invisible(gc())

# 5.2. H5AD
SeuratDisk::Convert(paste(project_dir, 'FINAL_ATLAS/1_Patients_Tumour_NormalMatched/Stromal_cells/CRCatlas_PTN_stromal.h5Seurat', sep='/'),
                    dest='h5ad', assay='RNA')
invisible(gc())

# 5.3. MTX
PTN_matrix = SeuratObject::GetAssayData(CRCatlas_PTN_stromal, slot='counts', assay='RNA')
dir.create(paste(project_dir, 'FINAL_ATLAS/1_Patients_Tumour_NormalMatched/Stromal_cells/mtx', sep='/'))
Matrix::writeMM(obj=PTN_matrix,
                file=paste(project_dir, 'FINAL_ATLAS/1_Patients_Tumour_NormalMatched/Stromal_cells/mtx/matrix.mtx', sep='/'))
write(x=rownames(PTN_matrix), file=paste(project_dir, 'FINAL_ATLAS/1_Patients_Tumour_NormalMatched/Stromal_cells/mtx/genes.tsv', sep='/'))
write(x=colnames(PTN_matrix), file=paste(project_dir, 'FINAL_ATLAS/1_Patients_Tumour_NormalMatched/Stromal_cells/mtx/barcodes.tsv', sep='/'))
invisible(gc())

# 5.4. Metadata CSV
write.csv(CRCatlas_PTN_stromal@meta.data,
          paste(project_dir, 'FINAL_ATLAS/1_Patients_Tumour_NormalMatched/Stromal_cells/metadata.csv', sep='/'),
          row.names=TRUE)

# 5.5. Remove epithelial object
remove(CRCatlas_PTN_stromal, PTN_matrix)
invisible(gc())


# 6. Create files for Myeloid cells:
CRCatlas_PTN_myeloid = subset(CRCatlas_P_TN,
                             cells=rownames(CRCatlas_P_TN[[]])[CRCatlas_P_TN$Annotation_0=='Myeloid cells'])
invisible(gc())
# 6.1. Seurat only with RNA assay: version 3.1.5.9900
SeuratDisk::SaveH5Seurat(CRCatlas_PTN_myeloid,
                         paste(project_dir, 'FINAL_ATLAS/1_Patients_Tumour_NormalMatched/Myeloid_cells/CRCatlas_PTN_myeloid.h5Seurat', sep='/'))
invisible(gc())

# 6.2. H5AD
SeuratDisk::Convert(paste(project_dir, 'FINAL_ATLAS/1_Patients_Tumour_NormalMatched/Myeloid_cells/CRCatlas_PTN_myeloid.h5Seurat', sep='/'),
                    dest='h5ad', assay='RNA')
invisible(gc())

# 6.3. MTX
PTN_matrix = SeuratObject::GetAssayData(CRCatlas_PTN_myeloid, slot='counts', assay='RNA')
dir.create(paste(project_dir, 'FINAL_ATLAS/1_Patients_Tumour_NormalMatched/Myeloid_cells/mtx', sep='/'))
Matrix::writeMM(obj=PTN_matrix,
                file=paste(project_dir, 'FINAL_ATLAS/1_Patients_Tumour_NormalMatched/Myeloid_cells/mtx/matrix.mtx', sep='/'))
write(x=rownames(PTN_matrix), file=paste(project_dir, 'FINAL_ATLAS/1_Patients_Tumour_NormalMatched/Myeloid_cells/mtx/genes.tsv', sep='/'))
write(x=colnames(PTN_matrix), file=paste(project_dir, 'FINAL_ATLAS/1_Patients_Tumour_NormalMatched/Myeloid_cells/mtx/barcodes.tsv', sep='/'))
invisible(gc())

# 6.4. Metadata CSV
write.csv(CRCatlas_PTN_myeloid@meta.data,
          paste(project_dir, 'FINAL_ATLAS/1_Patients_Tumour_NormalMatched/Myeloid_cells/metadata.csv', sep='/'),
          row.names=TRUE)

# 6.5. Remove epithelial object
remove(CRCatlas_PTN_myeloid, PTN_matrix)
invisible(gc())


# 7. Create files for B cells:
CRCatlas_PTN_Bcells = subset(CRCatlas_P_TN,
                            cells=rownames(CRCatlas_P_TN[[]])[CRCatlas_P_TN$Annotation_0=='Bcells'])
invisible(gc())
# 7.1. Seurat only with RNA assay: version 3.1.5.9900
SeuratDisk::SaveH5Seurat(CRCatlas_PTN_Bcells,
                         paste(project_dir, 'FINAL_ATLAS/1_Patients_Tumour_NormalMatched/B_cells/CRCatlas_PTN_Bcells.h5Seurat', sep='/'))
invisible(gc())

# 7.2. H5AD
SeuratDisk::Convert(paste(project_dir, 'FINAL_ATLAS/1_Patients_Tumour_NormalMatched/B_cells/CRCatlas_PTN_Bcells.h5Seurat', sep='/'),
                    dest='h5ad', assay='RNA')
invisible(gc())

# 7.3. MTX
PTN_matrix = SeuratObject::GetAssayData(CRCatlas_PTN_Bcells, slot='counts', assay='RNA')
dir.create(paste(project_dir, 'FINAL_ATLAS/1_Patients_Tumour_NormalMatched/B_cells/mtx', sep='/'))
Matrix::writeMM(obj=PTN_matrix,
                file=paste(project_dir, 'FINAL_ATLAS/1_Patients_Tumour_NormalMatched/B_cells/mtx/matrix.mtx', sep='/'))
write(x=rownames(PTN_matrix), file=paste(project_dir, 'FINAL_ATLAS/1_Patients_Tumour_NormalMatched/B_cells/mtx/genes.tsv', sep='/'))
write(x=colnames(PTN_matrix), file=paste(project_dir, 'FINAL_ATLAS/1_Patients_Tumour_NormalMatched/B_cells/mtx/barcodes.tsv', sep='/'))
invisible(gc())

# 7.4. Metadata CSV
write.csv(CRCatlas_PTN_Bcells@meta.data,
          paste(project_dir, 'FINAL_ATLAS/1_Patients_Tumour_NormalMatched/B_cells/metadata.csv', sep='/'),
          row.names=TRUE)

# 7.5. Remove epithelial object
remove(CRCatlas_PTN_Bcells, PTN_matrix)
invisible(gc())


# 8. Create files for T cells:
CRCatlas_PTN_Tcells = subset(CRCatlas_P_TN,
                            cells=rownames(CRCatlas_P_TN[[]])[CRCatlas_P_TN$Annotation_0=='Tcells'])
invisible(gc())
# 8.1. Seurat only with RNA assay: version 3.1.5.9900
SeuratDisk::SaveH5Seurat(CRCatlas_PTN_Tcells,
                         paste(project_dir, 'FINAL_ATLAS/1_Patients_Tumour_NormalMatched/T_cells/CRCatlas_PTN_Tcells.h5Seurat', sep='/'))
invisible(gc())

# 8.2. H5AD
SeuratDisk::Convert(paste(project_dir, 'FINAL_ATLAS/1_Patients_Tumour_NormalMatched/T_cells/CRCatlas_PTN_Tcells.h5Seurat', sep='/'),
                    dest='h5ad', assay='RNA')
invisible(gc())

# 8.3. MTX
PTN_matrix = SeuratObject::GetAssayData(CRCatlas_PTN_Tcells, slot='counts', assay='RNA')
dir.create(paste(project_dir, 'FINAL_ATLAS/1_Patients_Tumour_NormalMatched/T_cells/mtx', sep='/'))
Matrix::writeMM(obj=PTN_matrix,
                file=paste(project_dir, 'FINAL_ATLAS/1_Patients_Tumour_NormalMatched/T_cells/mtx/matrix.mtx', sep='/'))
write(x=rownames(PTN_matrix), file=paste(project_dir, 'FINAL_ATLAS/1_Patients_Tumour_NormalMatched/T_cells/mtx/genes.tsv', sep='/'))
write(x=colnames(PTN_matrix), file=paste(project_dir, 'FINAL_ATLAS/1_Patients_Tumour_NormalMatched/T_cells/mtx/barcodes.tsv', sep='/'))
invisible(gc())

# 8.4. Metadata CSV
write.csv(CRCatlas_PTN_Tcells@meta.data,
          paste(project_dir, 'FINAL_ATLAS/1_Patients_Tumour_NormalMatched/T_cells/metadata.csv', sep='/'),
          row.names=TRUE)

# 8.5. Remove epithelial object
remove(CRCatlas_PTN_Tcells, PTN_matrix)
invisible(gc())





# -----
# - Create Full Atlas (Patients: Tumour samples)
# -----
library(SeuratObject)


# 1. Load full dataset
CRCatlas_PH = SeuratDisk::LoadH5Seurat(paste(project_dir, 'FINAL_ATLAS/0_Patients_HealthyDonors/All_cells/CRCatlas_PH.h5Seurat', sep='/'))

# 2. Get only patients that have normal matched samples
CRCatlas_PT = subset(CRCatlas_PH, cells = rownames(CRCatlas_PH[[]])[CRCatlas_PH$Sample.State=='Tumor'])
remove(CRCatlas_PH)
invisible(gc())

# 3. Create files for all cells:
# 3.1. Seurat: version 3.1.5.9900
SeuratDisk::SaveH5Seurat(CRCatlas_PT,
                         paste(project_dir, 'FINAL_ATLAS/2_Patients_Tumour/All_cells/CRCatlas_PT.h5Seurat', sep='/'))
invisible(gc())

# 3.2. H5AD
SeuratDisk::Convert(paste(project_dir, 'FINAL_ATLAS/2_Patients_Tumour/All_cells/CRCatlas_PT.h5Seurat', sep='/'),
                    dest='h5ad', assay='RNA')
invisible(gc())

# 3.3. MTX
PT_matrix = SeuratObject::GetAssayData(CRCatlas_PT, slot='counts', assay='RNA')
dir.create(paste(project_dir, 'FINAL_ATLAS/2_Patients_Tumour/All_cells/mtx', sep='/'))
Matrix::writeMM(obj=PT_matrix,
                file=paste(project_dir, 'FINAL_ATLAS/2_Patients_Tumour/All_cells/mtx/matrix.mtx', sep='/'))
write(x=rownames(PT_matrix), file=paste(project_dir, 'FINAL_ATLAS/2_Patients_Tumour/All_cells/mtx/genes.tsv', sep='/'))
write(x=colnames(PT_matrix), file=paste(project_dir, 'FINAL_ATLAS/2_Patients_Tumour/All_cells/mtx/barcodes.tsv', sep='/'))
invisible(gc())

# 3.4. Metadata CSV
write.csv(CRCatlas_PT@meta.data,
          paste(project_dir, 'FINAL_ATLAS/2_Patients_Tumour/All_cells/metadata.csv', sep='/'),
          row.names=TRUE)


# 4. Create files for epithelial cells:
CRCatlas_PT_epithelial = subset(CRCatlas_PT,
                                 cells=rownames(CRCatlas_PT[[]])[CRCatlas_PT$Annotation_0=='Epithelial cells'])
invisible(gc())
# 4.1. Seurat only with RNA assay: version 3.1.5.9900
SeuratDisk::SaveH5Seurat(CRCatlas_PT_epithelial,
                         paste(project_dir, 'FINAL_ATLAS/2_Patients_Tumour/Epithelial_cells/CRCatlas_PT_epithelial.h5Seurat', sep='/'))
invisible(gc())

# 4.2. H5AD
SeuratDisk::Convert(paste(project_dir, 'FINAL_ATLAS/2_Patients_Tumour/Epithelial_cells/CRCatlas_PT_epithelial.h5Seurat', sep='/'),
                    dest='h5ad', assay='RNA')
invisible(gc())

# 4.3. MTX
PT_matrix = SeuratObject::GetAssayData(CRCatlas_PT_epithelial, slot='counts', assay='RNA')
dir.create(paste(project_dir, 'FINAL_ATLAS/2_Patients_Tumour/Epithelial_cells/mtx', sep='/'))
Matrix::writeMM(obj=PT_matrix,
                file=paste(project_dir, 'FINAL_ATLAS/2_Patients_Tumour/Epithelial_cells/mtx/matrix.mtx', sep='/'))
write(x=rownames(PT_matrix), file=paste(project_dir, 'FINAL_ATLAS/2_Patients_Tumour/Epithelial_cells/mtx/genes.tsv', sep='/'))
write(x=colnames(PT_matrix), file=paste(project_dir, 'FINAL_ATLAS/2_Patients_Tumour/Epithelial_cells/mtx/barcodes.tsv', sep='/'))
invisible(gc())

# 4.4. Metadata CSV
write.csv(CRCatlas_PT_epithelial@meta.data,
          paste(project_dir, 'FINAL_ATLAS/2_Patients_Tumour/Epithelial_cells/metadata.csv', sep='/'),
          row.names=TRUE)

# 4.5. Remove epithelial object
remove(CRCatlas_PT_epithelial, PT_matrix)
invisible(gc())


# 5. Create files for Stromal cells:
CRCatlas_PT_stromal = subset(CRCatlas_PT,
                              cells=rownames(CRCatlas_PT[[]])[CRCatlas_PT$Annotation_0=='Stromal cells'])
invisible(gc())
# 5.1. Seurat only with RNA assay: version 3.1.5.9900
SeuratDisk::SaveH5Seurat(CRCatlas_PT_stromal,
                         paste(project_dir, 'FINAL_ATLAS/2_Patients_Tumour/Stromal_cells/CRCatlas_PT_stromal.h5Seurat', sep='/'))
invisible(gc())

# 5.2. H5AD
SeuratDisk::Convert(paste(project_dir, 'FINAL_ATLAS/2_Patients_Tumour/Stromal_cells/CRCatlas_PT_stromal.h5Seurat', sep='/'),
                    dest='h5ad', assay='RNA')
invisible(gc())

# 5.3. MTX
PT_matrix = SeuratObject::GetAssayData(CRCatlas_PT_stromal, slot='counts', assay='RNA')
dir.create(paste(project_dir, 'FINAL_ATLAS/2_Patients_Tumour/Stromal_cells/mtx', sep='/'))
Matrix::writeMM(obj=PT_matrix,
                file=paste(project_dir, 'FINAL_ATLAS/2_Patients_Tumour/Stromal_cells/mtx/matrix.mtx', sep='/'))
write(x=rownames(PT_matrix), file=paste(project_dir, 'FINAL_ATLAS/2_Patients_Tumour/Stromal_cells/mtx/genes.tsv', sep='/'))
write(x=colnames(PT_matrix), file=paste(project_dir, 'FINAL_ATLAS/2_Patients_Tumour/Stromal_cells/mtx/barcodes.tsv', sep='/'))
invisible(gc())

# 5.4. Metadata CSV
write.csv(CRCatlas_PT_stromal@meta.data,
          paste(project_dir, 'FINAL_ATLAS/2_Patients_Tumour/Stromal_cells/metadata.csv', sep='/'),
          row.names=TRUE)

# 5.5. Remove epithelial object
remove(CRCatlas_PT_stromal, PT_matrix)
invisible(gc())


# 6. Create files for Myeloid cells:
CRCatlas_PT_myeloid = subset(CRCatlas_PT,
                              cells=rownames(CRCatlas_PT[[]])[CRCatlas_PT$Annotation_0=='Myeloid cells'])
invisible(gc())
# 6.1. Seurat only with RNA assay: version 3.1.5.9900
SeuratDisk::SaveH5Seurat(CRCatlas_PT_myeloid,
                         paste(project_dir, 'FINAL_ATLAS/2_Patients_Tumour/Myeloid_cells/CRCatlas_PT_myeloid.h5Seurat', sep='/'))
invisible(gc())

# 6.2. H5AD
SeuratDisk::Convert(paste(project_dir, 'FINAL_ATLAS/2_Patients_Tumour/Myeloid_cells/CRCatlas_PT_myeloid.h5Seurat', sep='/'),
                    dest='h5ad', assay='RNA')
invisible(gc())

# 6.3. MTX
PT_matrix = SeuratObject::GetAssayData(CRCatlas_PT_myeloid, slot='counts', assay='RNA')
dir.create(paste(project_dir, 'FINAL_ATLAS/2_Patients_Tumour/Myeloid_cells/mtx', sep='/'))
Matrix::writeMM(obj=PT_matrix,
                file=paste(project_dir, 'FINAL_ATLAS/2_Patients_Tumour/Myeloid_cells/mtx/matrix.mtx', sep='/'))
write(x=rownames(PT_matrix), file=paste(project_dir, 'FINAL_ATLAS/2_Patients_Tumour/Myeloid_cells/mtx/genes.tsv', sep='/'))
write(x=colnames(PT_matrix), file=paste(project_dir, 'FINAL_ATLAS/2_Patients_Tumour/Myeloid_cells/mtx/barcodes.tsv', sep='/'))
invisible(gc())

# 6.4. Metadata CSV
write.csv(CRCatlas_PT_myeloid@meta.data,
          paste(project_dir, 'FINAL_ATLAS/2_Patients_Tumour/Myeloid_cells/metadata.csv', sep='/'),
          row.names=TRUE)

# 6.5. Remove epithelial object
remove(CRCatlas_PT_myeloid, PT_matrix)
invisible(gc())


# 7. Create files for B cells:
CRCatlas_PT_Bcells = subset(CRCatlas_PT,
                             cells=rownames(CRCatlas_PT[[]])[CRCatlas_PT$Annotation_0=='Bcells'])
invisible(gc())
# 7.1. Seurat only with RNA assay: version 3.1.5.9900
SeuratDisk::SaveH5Seurat(CRCatlas_PT_Bcells,
                         paste(project_dir, 'FINAL_ATLAS/2_Patients_Tumour/B_cells/CRCatlas_PT_Bcells.h5Seurat', sep='/'))
invisible(gc())

# 7.2. H5AD
SeuratDisk::Convert(paste(project_dir, 'FINAL_ATLAS/2_Patients_Tumour/B_cells/CRCatlas_PT_Bcells.h5Seurat', sep='/'),
                    dest='h5ad', assay='RNA')
invisible(gc())

# 7.3. MTX
PT_matrix = SeuratObject::GetAssayData(CRCatlas_PT_Bcells, slot='counts', assay='RNA')
dir.create(paste(project_dir, 'FINAL_ATLAS/2_Patients_Tumour/B_cells/mtx', sep='/'))
Matrix::writeMM(obj=PT_matrix,
                file=paste(project_dir, 'FINAL_ATLAS/2_Patients_Tumour/B_cells/mtx/matrix.mtx', sep='/'))
write(x=rownames(PT_matrix), file=paste(project_dir, 'FINAL_ATLAS/2_Patients_Tumour/B_cells/mtx/genes.tsv', sep='/'))
write(x=colnames(PT_matrix), file=paste(project_dir, 'FINAL_ATLAS/2_Patients_Tumour/B_cells/mtx/barcodes.tsv', sep='/'))
invisible(gc())

# 7.4. Metadata CSV
write.csv(CRCatlas_PT_Bcells@meta.data,
          paste(project_dir, 'FINAL_ATLAS/2_Patients_Tumour/B_cells/metadata.csv', sep='/'),
          row.names=TRUE)

# 7.5. Remove epithelial object
remove(CRCatlas_PT_Bcells, PT_matrix)
invisible(gc())


# 8. Create files for T cells:
CRCatlas_PT_Tcells = subset(CRCatlas_PT,
                             cells=rownames(CRCatlas_PT[[]])[CRCatlas_PT$Annotation_0=='Tcells'])
invisible(gc())
# 8.1. Seurat only with RNA assay: version 3.1.5.9900
SeuratDisk::SaveH5Seurat(CRCatlas_PT_Tcells,
                         paste(project_dir, 'FINAL_ATLAS/2_Patients_Tumour/T_cells/CRCatlas_PT_Tcells.h5Seurat', sep='/'))
invisible(gc())

# 8.2. H5AD
SeuratDisk::Convert(paste(project_dir, 'FINAL_ATLAS/2_Patients_Tumour/T_cells/CRCatlas_PT_Tcells.h5Seurat', sep='/'),
                    dest='h5ad', assay='RNA')
invisible(gc())

# 8.3. MTX
PT_matrix = SeuratObject::GetAssayData(CRCatlas_PT_Tcells, slot='counts', assay='RNA')
dir.create(paste(project_dir, 'FINAL_ATLAS/2_Patients_Tumour/T_cells/mtx', sep='/'))
Matrix::writeMM(obj=PT_matrix,
                file=paste(project_dir, 'FINAL_ATLAS/2_Patients_Tumour/T_cells/mtx/matrix.mtx', sep='/'))
write(x=rownames(PT_matrix), file=paste(project_dir, 'FINAL_ATLAS/2_Patients_Tumour/T_cells/mtx/genes.tsv', sep='/'))
write(x=colnames(PT_matrix), file=paste(project_dir, 'FINAL_ATLAS/2_Patients_Tumour/T_cells/mtx/barcodes.tsv', sep='/'))
invisible(gc())

# 8.4. Metadata CSV
write.csv(CRCatlas_PT_Tcells@meta.data,
          paste(project_dir, 'FINAL_ATLAS/2_Patients_Tumour/T_cells/metadata.csv', sep='/'),
          row.names=TRUE)

# 8.5. Remove epithelial object
remove(CRCatlas_PT_Tcells, PT_matrix)
invisible(gc())





# -----
# - Create Full Atlas (Patients: Healthy + Normal matched)
# -----
library(SeuratObject)


# 1. Load full dataset
CRCatlas_PH = SeuratDisk::LoadH5Seurat(paste(project_dir, 'FINAL_ATLAS/0_Patients_HealthyDonors/All_cells/CRCatlas_PH.h5Seurat', sep='/'))

# 2. Get only patients that have normal matched samples
CRCatlas_H = subset(CRCatlas_PH, cells = rownames(CRCatlas_PH[[]])[CRCatlas_PH$Sample.State!='Tumor'])
remove(CRCatlas_PH)
invisible(gc())

# 3. Create files for all cells:
# 3.1. Seurat: version 3.1.5.9900
SeuratDisk::SaveH5Seurat(CRCatlas_H,
                         paste(project_dir, 'FINAL_ATLAS/3_HealthyDonors_NormalMatched/All_cells/CRCatlas_H.h5Seurat', sep='/'))
invisible(gc())

# 3.2. H5AD
SeuratDisk::Convert(paste(project_dir, 'FINAL_ATLAS/3_HealthyDonors_NormalMatched/All_cells/CRCatlas_H.h5Seurat', sep='/'),
                    dest='h5ad', assay='RNA')
invisible(gc())

# 3.3. MTX
H_matrix = SeuratObject::GetAssayData(CRCatlas_H, slot='counts', assay='RNA')
dir.create(paste(project_dir, 'FINAL_ATLAS/3_HealthyDonors_NormalMatched/All_cells/mtx', sep='/'))
Matrix::writeMM(obj=H_matrix,
                file=paste(project_dir, 'FINAL_ATLAS/3_HealthyDonors_NormalMatched/All_cells/mtx/matrix.mtx', sep='/'))
write(x=rownames(H_matrix), file=paste(project_dir, 'FINAL_ATLAS/3_HealthyDonors_NormalMatched/All_cells/mtx/genes.tsv', sep='/'))
write(x=colnames(H_matrix), file=paste(project_dir, 'FINAL_ATLAS/3_HealthyDonors_NormalMatched/All_cells/mtx/barcodes.tsv', sep='/'))
invisible(gc())

# 3.4. Metadata CSV
write.csv(CRCatlas_H@meta.data,
          paste(project_dir, 'FINAL_ATLAS/3_HealthyDonors_NormalMatched/All_cells/metadata.csv', sep='/'),
          row.names=TRUE)


# 4. Create files for epithelial cells:
CRCatlas_H_epithelial = subset(CRCatlas_H,
                                cells=rownames(CRCatlas_H[[]])[CRCatlas_H$Annotation_0=='Epithelial cells'])
invisible(gc())
# 4.1. Seurat only with RNA assay: version 3.1.5.9900
SeuratDisk::SaveH5Seurat(CRCatlas_H_epithelial,
                         paste(project_dir, 'FINAL_ATLAS/3_HealthyDonors_NormalMatched/Epithelial_cells/CRCatlas_H_epithelial.h5Seurat', sep='/'))
invisible(gc())

# 4.2. H5AD
SeuratDisk::Convert(paste(project_dir, 'FINAL_ATLAS/3_HealthyDonors_NormalMatched/Epithelial_cells/CRCatlas_H_epithelial.h5Seurat', sep='/'),
                    dest='h5ad', assay='RNA')
invisible(gc())

# 4.3. MTX
H_matrix = SeuratObject::GetAssayData(CRCatlas_H_epithelial, slot='counts', assay='RNA')
dir.create(paste(project_dir, 'FINAL_ATLAS/3_HealthyDonors_NormalMatched/Epithelial_cells/mtx', sep='/'))
Matrix::writeMM(obj=H_matrix,
                file=paste(project_dir, 'FINAL_ATLAS/3_HealthyDonors_NormalMatched/Epithelial_cells/mtx/matrix.mtx', sep='/'))
write(x=rownames(H_matrix), file=paste(project_dir, 'FINAL_ATLAS/3_HealthyDonors_NormalMatched/Epithelial_cells/mtx/genes.tsv', sep='/'))
write(x=colnames(H_matrix), file=paste(project_dir, 'FINAL_ATLAS/3_HealthyDonors_NormalMatched/Epithelial_cells/mtx/barcodes.tsv', sep='/'))
invisible(gc())

# 4.4. Metadata CSV
write.csv(CRCatlas_H_epithelial@meta.data,
          paste(project_dir, 'FINAL_ATLAS/3_HealthyDonors_NormalMatched/Epithelial_cells/metadata.csv', sep='/'),
          row.names=TRUE)

# 4.5. Remove epithelial object
remove(CRCatlas_H_epithelial, H_matrix)
invisible(gc())


# 5. Create files for Stromal cells:
CRCatlas_H_stromal = subset(CRCatlas_H,
                             cells=rownames(CRCatlas_H[[]])[CRCatlas_H$Annotation_0=='Stromal cells'])
invisible(gc())
# 5.1. Seurat only with RNA assay: version 3.1.5.9900
SeuratDisk::SaveH5Seurat(CRCatlas_H_stromal,
                         paste(project_dir, 'FINAL_ATLAS/3_HealthyDonors_NormalMatched/Stromal_cells/CRCatlas_H_stromal.h5Seurat', sep='/'))
invisible(gc())

# 5.2. H5AD
SeuratDisk::Convert(paste(project_dir, 'FINAL_ATLAS/3_HealthyDonors_NormalMatched/Stromal_cells/CRCatlas_H_stromal.h5Seurat', sep='/'),
                    dest='h5ad', assay='RNA')
invisible(gc())

# 5.3. MTX
H_matrix = SeuratObject::GetAssayData(CRCatlas_H_stromal, slot='counts', assay='RNA')
dir.create(paste(project_dir, 'FINAL_ATLAS/3_HealthyDonors_NormalMatched/Stromal_cells/mtx', sep='/'))
Matrix::writeMM(obj=H_matrix,
                file=paste(project_dir, 'FINAL_ATLAS/3_HealthyDonors_NormalMatched/Stromal_cells/mtx/matrix.mtx', sep='/'))
write(x=rownames(H_matrix), file=paste(project_dir, 'FINAL_ATLAS/3_HealthyDonors_NormalMatched/Stromal_cells/mtx/genes.tsv', sep='/'))
write(x=colnames(H_matrix), file=paste(project_dir, 'FINAL_ATLAS/3_HealthyDonors_NormalMatched/Stromal_cells/mtx/barcodes.tsv', sep='/'))
invisible(gc())

# 5.4. Metadata CSV
write.csv(CRCatlas_H_stromal@meta.data,
          paste(project_dir, 'FINAL_ATLAS/3_HealthyDonors_NormalMatched/Stromal_cells/metadata.csv', sep='/'),
          row.names=TRUE)

# 5.5. Remove epithelial object
remove(CRCatlas_H_stromal, PT_matrix)
invisible(gc())


# 6. Create files for Myeloid cells:
CRCatlas_H_myeloid = subset(CRCatlas_H,
                             cells=rownames(CRCatlas_H[[]])[CRCatlas_H$Annotation_0=='Myeloid cells'])
invisible(gc())
# 6.1. Seurat only with RNA assay: version 3.1.5.9900
SeuratDisk::SaveH5Seurat(CRCatlas_H_myeloid,
                         paste(project_dir, 'FINAL_ATLAS/3_HealthyDonors_NormalMatched/Myeloid_cells/CRCatlas_H_myeloid.h5Seurat', sep='/'))
invisible(gc())

# 6.2. H5AD
SeuratDisk::Convert(paste(project_dir, 'FINAL_ATLAS/3_HealthyDonors_NormalMatched/Myeloid_cells/CRCatlas_H_myeloid.h5Seurat', sep='/'),
                    dest='h5ad', assay='RNA')
invisible(gc())

# 6.3. MTX
H_matrix = SeuratObject::GetAssayData(CRCatlas_H_myeloid, slot='counts', assay='RNA')
dir.create(paste(project_dir, 'FINAL_ATLAS/3_HealthyDonors_NormalMatched/Myeloid_cells/mtx', sep='/'))
Matrix::writeMM(obj=H_matrix,
                file=paste(project_dir, 'FINAL_ATLAS/3_HealthyDonors_NormalMatched/Myeloid_cells/mtx/matrix.mtx', sep='/'))
write(x=rownames(H_matrix), file=paste(project_dir, 'FINAL_ATLAS/3_HealthyDonors_NormalMatched/Myeloid_cells/mtx/genes.tsv', sep='/'))
write(x=colnames(H_matrix), file=paste(project_dir, 'FINAL_ATLAS/3_HealthyDonors_NormalMatched/Myeloid_cells/mtx/barcodes.tsv', sep='/'))
invisible(gc())

# 6.4. Metadata CSV
write.csv(CRCatlas_H_myeloid@meta.data,
          paste(project_dir, 'FINAL_ATLAS/3_HealthyDonors_NormalMatched/Myeloid_cells/metadata.csv', sep='/'),
          row.names=TRUE)

# 6.5. Remove epithelial object
remove(CRCatlas_H_myeloid, H_matrix)
invisible(gc())


# 7. Create files for B cells:
CRCatlas_H_Bcells = subset(CRCatlas_H,
                            cells=rownames(CRCatlas_H[[]])[CRCatlas_H$Annotation_0=='Bcells'])
invisible(gc())
# 7.1. Seurat only with RNA assay: version 3.1.5.9900
SeuratDisk::SaveH5Seurat(CRCatlas_H_Bcells,
                         paste(project_dir, 'FINAL_ATLAS/3_HealthyDonors_NormalMatched/B_cells/CRCatlas_H_Bcells.h5Seurat', sep='/'))
invisible(gc())

# 7.2. H5AD
SeuratDisk::Convert(paste(project_dir, 'FINAL_ATLAS/3_HealthyDonors_NormalMatched/B_cells/CRCatlas_H_Bcells.h5Seurat', sep='/'),
                    dest='h5ad', assay='RNA')
invisible(gc())

# 7.3. MTX
H_matrix = SeuratObject::GetAssayData(CRCatlas_H_Bcells, slot='counts', assay='RNA')
dir.create(paste(project_dir, 'FINAL_ATLAS/3_HealthyDonors_NormalMatched/B_cells/mtx', sep='/'))
Matrix::writeMM(obj=H_matrix,
                file=paste(project_dir, 'FINAL_ATLAS/3_HealthyDonors_NormalMatched/B_cells/mtx/matrix.mtx', sep='/'))
write(x=rownames(H_matrix), file=paste(project_dir, 'FINAL_ATLAS/3_HealthyDonors_NormalMatched/B_cells/mtx/genes.tsv', sep='/'))
write(x=colnames(H_matrix), file=paste(project_dir, 'FINAL_ATLAS/3_HealthyDonors_NormalMatched/B_cells/mtx/barcodes.tsv', sep='/'))
invisible(gc())

# 7.4. Metadata CSV
write.csv(CRCatlas_H_Bcells@meta.data,
          paste(project_dir, 'FINAL_ATLAS/3_HealthyDonors_NormalMatched/B_cells/metadata.csv', sep='/'),
          row.names=TRUE)

# 7.5. Remove epithelial object
remove(CRCatlas_H_Bcells, H_matrix)
invisible(gc())


# 8. Create files for T cells:
CRCatlas_H_Tcells = subset(CRCatlas_H,
                            cells=rownames(CRCatlas_H[[]])[CRCatlas_H$Annotation_0=='Tcells'])
invisible(gc())
# 8.1. Seurat only with RNA assay: version 3.1.5.9900
SeuratDisk::SaveH5Seurat(CRCatlas_H_Tcells,
                         paste(project_dir, 'FINAL_ATLAS/3_HealthyDonors_NormalMatched/T_cells/CRCatlas_H_Tcells.h5Seurat', sep='/'))
invisible(gc())

# 8.2. H5AD
SeuratDisk::Convert(paste(project_dir, 'FINAL_ATLAS/3_HealthyDonors_NormalMatched/T_cells/CRCatlas_H_Tcells.h5Seurat', sep='/'),
                    dest='h5ad', assay='RNA')
invisible(gc())

# 8.3. MTX
H_matrix = SeuratObject::GetAssayData(CRCatlas_H_Tcells, slot='counts', assay='RNA')
dir.create(paste(project_dir, 'FINAL_ATLAS/3_HealthyDonors_NormalMatched/T_cells/mtx', sep='/'))
Matrix::writeMM(obj=H_matrix,
                file=paste(project_dir, 'FINAL_ATLAS/3_HealthyDonors_NormalMatched/T_cells/mtx/matrix.mtx', sep='/'))
write(x=rownames(H_matrix), file=paste(project_dir, 'FINAL_ATLAS/3_HealthyDonors_NormalMatched/T_cells/mtx/genes.tsv', sep='/'))
write(x=colnames(H_matrix), file=paste(project_dir, 'FINAL_ATLAS/3_HealthyDonors_NormalMatched/T_cells/mtx/barcodes.tsv', sep='/'))
invisible(gc())

# 8.4. Metadata CSV
write.csv(CRCatlas_H_Tcells@meta.data,
          paste(project_dir, 'FINAL_ATLAS/3_HealthyDonors_NormalMatched/T_cells/metadata.csv', sep='/'),
          row.names=TRUE)

# 8.5. Remove epithelial object
remove(CRCatlas_H_Tcells, H_matrix)
invisible(gc())








