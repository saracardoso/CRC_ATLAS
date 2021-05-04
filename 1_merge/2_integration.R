project_dir = '/home/scardoso/Documents/PhD/CRC_ATLAS'

# Cell-cycle markers:
cell_cycle_markers = paste(project_dir, 'utils/cycle_markers.rda', sep='/')
load(cell_cycle_markers)

# Load previously merged dataset and separate them by dataset:
datasets = Seurat::SplitObject(SeuratDisk::LoadH5Seurat(paste(project_dir, '1_merge/CRC.h5Seurat', sep='/')),
                               split.by='dataset')
gc()

# Normalize and find variable features for each dataset:
for(dts_name in names(datasets)) {
  datasets[[dts_name]] <- Seurat::NormalizeData(datasets[[dts_name]])
  datasets[[dts_name]] <- Seurat::FindVariableFeatures(datasets[[dts_name]], selection.method = "vst", nfeatures = 2000)
  gc()
}
gc()

# Select integration features:
features <- Seurat::SelectIntegrationFeatures(object.list = datasets)
gc()

# Scale data and run a PCA do reduce dimensionality of datasets to perform faster and less computationally heavy integration:
for(dts_name in names(datasets)) {
  message(dts_name)
  datasets[[dts_name]] <- Seurat::CellCycleScoring(datasets[[dts_name]], g2m.features=g2m_genes, s.features=s_genes)
  datasets[[dts_name]]$CC.Difference <-  datasets[[dts_name]]$S.Score - datasets[[dts_name]]$G2M.Score
  datasets[[dts_name]] <- Seurat::ScaleData(datasets[[dts_name]], features = features, verbose = TRUE,
                                            vars.to.regress=c('CC.Difference', 'percent.mitochondrial_RNA'))
  datasets[[dts_name]] <- Seurat::RunPCA(datasets[[dts_name]], features = features, verbose = TRUE)
}
gc()
saveRDS(datasets, '/home/scardoso/Documents/PhD/CRC_ATLAS/datasets.Rdata')

# Find integration achors and integrate the data:
anchors <- Seurat::FindIntegrationAnchors(object.list=datasets, reference=c(2, 4), reduction = "rpca", dims = 1:50)
gc()
CRCatlas_integrated <- Seurat::IntegrateData(anchorset = anchors, dims = 1:50)
gc()

# Add a new metadata variable:
CRCatlas_integrated[['state']] = CRCatlas_integrated$sample_type
CRCatlas_integrated$state[CRCatlas_integrated$state=='Normal Tissue'] = 'Normal'
CRCatlas_integrated$state[CRCatlas_integrated$state=='Tumour Border'] = 'Tumor'
CRCatlas_integrated$state[CRCatlas_integrated$state=='Tumour Core'] = 'Tumor'

# Store integrated CRCatlas:
SeuratDisk::SaveH5Seurat(CRCatlas_integrated, '/home/scardoso/Documents/PhD/CRC_ATLAS/1_merge/CRC_integrated.h5Seurat')

