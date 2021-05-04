project_dir = '/home/scardoso/Documents/PhD/CRC_ATLAS'

CRC_Qian_data_dir = paste(project_dir, 'original_datasets/CRC_qian/data', sep='/')
CRC_Qian_metadata = paste(project_dir, 'original_datasets/CRC_qian/metadata.csv', sep='/')

# Cell-cycle markers:
cell_cycle_markers = paste(project_dir, 'utils/cycle_markers.rda', sep='/')
load(cell_cycle_markers)

# Read data to a seurat object:
seurat_data = Seurat::Read10X(data.dir=CRC_Qian_data_dir)
CRC_Qian = Seurat::CreateSeuratObject(counts=seurat_data, project='CRCQian')

# Read metadata file:
CRC_Qian_original_metadata = read.csv(CRC_Qian_metadata, row.names=1)

# Update metadata:
colnames(CRC_Qian@meta.data) = c('sample', 'nUMI', 'nGene')
CRC_Qian[['patient']] = CRC_Qian_original_metadata[rownames(CRC_Qian@meta.data), 'PatientNumber']
CRC_Qian[['sample_type']] = CRC_Qian_original_metadata[rownames(CRC_Qian@meta.data), 'TumorSite']
CRC_Qian[['sample_type']][CRC_Qian[['sample_type']]=='C'] = 'Tumour Core'
CRC_Qian[['sample_type']][CRC_Qian[['sample_type']]=='B'] = 'Tumour Border'
CRC_Qian[['sample_type']][CRC_Qian[['sample_type']]=='N'] = 'Normal Tissue'
CRC_Qian[['original_annotation']] = CRC_Qian_original_metadata[rownames(CRC_Qian@meta.data), 'CellType']

# Store original data as seurat file:
SeuratDisk::SaveH5Seurat(CRC_Qian, paste(project_dir, 'original_datasets/CRC_qian/original_data.h5Seurat', sep='/'))

# Percentage of mitochondrial RNA:
features = grep(pattern='^MT-', x=rownames(CRC_Qian@assays$RNA@counts), value = TRUE)
percent.featureset = colSums(as.matrix(Seurat::GetAssayData(CRC_Qian, assay='RNA', slot="counts")[features, , drop = FALSE]))/
  CRC_Qian$nUMI * 100
CRC_Qian[['percent.mitochondrial_RNA']] = percent.featureset

# Complexity
CRC_Qian[['log10.complexity']] = log10(CRC_Qian$nGene) / log10(CRC_Qian$nUMI)


# Quality control - Before
# Number of cell counts per sample
ggplot2::ggplot(CRC_Qian@meta.data, ggplot2::aes(x=sample, fill=sample)) + 
  ggplot2::geom_bar() + ggplot2::theme_classic() +
  ggplot2::theme(axis.text.x = ggplot2::element_text(angle = 45, vjust = 1, hjust=1)) +
  ggplot2::theme(plot.title = ggplot2::element_text(hjust=0.5, face="bold"), legend.position = "none") +
  ggplot2::ggtitle("NCells")
# UMI counts per cell
ggplot2::ggplot(CRC_Qian@meta.data, ggplot2::aes(x=nUMI, fill=sample)) + 
  ggplot2::geom_density(alpha = 0.2) + 
  ggplot2::scale_x_log10() + 
  ggplot2::theme_classic() +
  ggplot2::xlab('Nº UMIs') +
  ggplot2::ylab("Cell density") +
  ggplot2::geom_vline(xintercept = 1000) +
  ggplot2::ggtitle('UMI counts per cell')
# Number of detected genes per cell
ggplot2::ggplot(CRC_Qian@meta.data, ggplot2::aes(x=nGene, fill=sample)) + 
  ggplot2::geom_density(alpha = 0.2) + 
  ggplot2::theme_classic() + ggplot2::xlab('Nº Genes') + ggplot2::ylab('Density') +
  ggplot2::scale_x_log10() + 
  ggplot2::geom_vline(xintercept = 250) +
  ggplot2::ggtitle('Number of detected genes per cell')
# Distribution of genes detected per cell:
ggplot2::ggplot(CRC_Qian@meta.data, ggplot2::aes(y=log10(nGene), fill=sample)) + 
  ggplot2::geom_boxplot() + 
  ggplot2::theme_classic() +
  ggplot2::theme(plot.title = ggplot2::element_text(hjust=0.5, face="bold")) +
  ggplot2::ggtitle('Distribution of genes detected per cell:')
# UMIs vs Genes
ggplot2::ggplot(CRC_Qian@meta.data, ggplot2::aes(x=nUMI, y=nGene, color=percent.mitochondrial_RNA)) + 
  ggplot2::geom_point() + 
  ggplot2::scale_colour_gradient(low = "gray90", high = "black") +
  ggplot2::stat_smooth(method=lm) +
  ggplot2::scale_x_log10() + 
  ggplot2::scale_y_log10() + 
  ggplot2::theme_classic() +
  ggplot2::geom_vline(xintercept = 1000) +
  ggplot2::geom_hline(yintercept = 250) +
  ggplot2::xlab('Nº UMIs') + ggplot2::ylab('Nº Genes') +
  ggplot2::ggtitle('UMIs vs Genes')
# Percentage of mitochondrial gene counts per cell
ggplot2::ggplot(CRC_Qian@meta.data, ggplot2::aes(x=percent.mitochondrial_RNA, fill=sample)) + 
  ggplot2::geom_density(alpha = 0.2) + 
  ggplot2::scale_x_log10() + 
  ggplot2::theme_classic() +
  ggplot2::geom_vline(xintercept = 20) +
  ggplot2::ggtitle('Percentage of mitochondrial gene counts per cell')
# Complexity
ggplot2::ggplot(CRC_Qian@meta.data, ggplot2::aes(x=log10.complexity, fill=sample)) +
  ggplot2::geom_density(alpha = 0.2) +
  ggplot2::theme_classic() +
  ggplot2::geom_vline(xintercept = 0.8) +
  ggplot2::ggtitle('Complexity')

# Filter out low quality cells:
cells_keep = rownames(CRC_Qian@meta.data)[(CRC_Qian$nUMI >= 1000) &
                                            (CRC_Qian$nGene >= 250) &
                                            (CRC_Qian$nGene <= 6000) &
                                            (CRC_Qian$log10.complexity > 0.80) &
                                            (CRC_Qian$percent.mitochondrial_RNA <= 20)]
filtered_seurat = subset(CRC_Qian, cells=cells_keep)
# Filter out low-expressed genes:
counts = Seurat::GetAssayData(object=filtered_seurat, slot="counts")
nonzero <- counts > 0
keep_genes <- Matrix::rowSums(nonzero) >= 10
filtered_counts <- counts[keep_genes, ]
# Creation of filtered dataset:
CRC_Qian_filtered = Seurat::CreateSeuratObject(filtered_counts, meta.data=filtered_seurat@meta.data)

# Quality control - After
# Number of cell counts per sample
ggplot2::ggplot(CRC_Qian_filtered@meta.data, ggplot2::aes(x=sample, fill=sample)) + 
  ggplot2::geom_bar() + ggplot2::theme_classic() +
  ggplot2::theme(axis.text.x = ggplot2::element_text(angle = 45, vjust = 1, hjust=1)) +
  ggplot2::theme(plot.title = ggplot2::element_text(hjust=0.5, face="bold"), legend.position = "none") +
  ggplot2::ggtitle("NCells")
# UMI counts per cell
ggplot2::ggplot(CRC_Qian_filtered@meta.data, ggplot2::aes(x=nUMI, fill=sample)) + 
  ggplot2::geom_density(alpha = 0.2) + 
  ggplot2::scale_x_log10() + 
  ggplot2::theme_classic() +
  ggplot2::xlab('Nº UMIs') +
  ggplot2::ylab("Cell density") +
  ggplot2::geom_vline(xintercept = 1000) +
  ggplot2::ggtitle('UMI counts per cell')
# Number of detected genes per cell
ggplot2::ggplot(CRC_Qian_filtered@meta.data, ggplot2::aes(x=nGene, fill=sample)) + 
  ggplot2::geom_density(alpha = 0.2) + 
  ggplot2::theme_classic() + ggplot2::xlab('Nº Genes') + ggplot2::ylab('Density') +
  ggplot2::scale_x_log10() + 
  ggplot2::geom_vline(xintercept = 250) +
  ggplot2::ggtitle('Number of detected genes per cell')
# Distribution of genes detected per cell:
ggplot2::ggplot(CRC_Qian_filtered@meta.data, ggplot2::aes(y=log10(nGene), fill=sample)) + 
  ggplot2::geom_boxplot() + 
  ggplot2::theme_classic() +
  ggplot2::theme(plot.title = ggplot2::element_text(hjust=0.5, face="bold")) +
  ggplot2::ggtitle('Distribution of genes detected per cell:')
# UMIs vs Genes
ggplot2::ggplot(CRC_Qian_filtered@meta.data, ggplot2::aes(x=nUMI, y=nGene, color=percent.mitochondrial_RNA)) + 
  ggplot2::geom_point() + 
  ggplot2::scale_colour_gradient(low = "gray90", high = "black") +
  ggplot2::stat_smooth(method=lm) +
  ggplot2::scale_x_log10() + 
  ggplot2::scale_y_log10() + 
  ggplot2::theme_classic() +
  ggplot2::geom_vline(xintercept = 1000) +
  ggplot2::geom_hline(yintercept = 250) +
  ggplot2::xlab('Nº UMIs') + ggplot2::ylab('Nº Genes') +
  ggplot2::ggtitle('UMIs vs Genes')
# Percentage of mitochondrial gene counts per cell
ggplot2::ggplot(CRC_Qian_filtered@meta.data, ggplot2::aes(x=percent.mitochondrial_RNA, fill=sample)) + 
  ggplot2::geom_density(alpha = 0.2) + 
  ggplot2::scale_x_log10() + 
  ggplot2::theme_classic() +
  ggplot2::geom_vline(xintercept = 20) +
  ggplot2::ggtitle('Percentage of mitochondrial gene counts per cell')
# Complexity
ggplot2::ggplot(CRC_Qian_filtered@meta.data, ggplot2::aes(x=log10.complexity, fill=sample)) +
  ggplot2::geom_density(alpha = 0.2) +
  ggplot2::theme_classic() +
  ggplot2::geom_vline(xintercept = 0.8) +
  ggplot2::ggtitle('Complexity')


# Cell Cycle scoring
CRC_Qian_filtered = Seurat::NormalizeData(CRC_Qian_filtered, assay='RNA')
CRC_Qian_filtered = Seurat::CellCycleScoring(CRC_Qian_filtered, g2m.features=g2m_genes, s.features=s_genes, assay='RNA')
# Perform PCA to assess if cell cycle is a major source of variation
CRC_Qian_filtered = Seurat::FindVariableFeatures(CRC_Qian_filtered, selection.method="vst", nfeatures=2000, verbose=FALSE)
CRC_Qian_filtered = Seurat::ScaleData(CRC_Qian_filtered)
CRC_Qian_filtered = Seurat::RunPCA(CRC_Qian_filtered)
# Plot results:
Seurat::DimPlot(CRC_Qian_filtered, reduction="pca", group.by='Phase', pt.size=1)
Seurat::DimPlot(CRC_Qian_filtered, reduction="pca", group.by='Phase', split.by = 'Phase', pt.size=1)


# Visualize data in clusters:
CRC_Qian_filtered = Seurat::SCTransform(CRC_Qian_filtered, variable.features.n=2000, ncells=2000, conserve.memory=TRUE)
CRC_Qian_filtered = Seurat::RunPCA(CRC_Qian_filtered, verbose = FALSE)
CRC_Qian_filtered = Seurat::RunUMAP(CRC_Qian_filtered, dims = 1:30, verbose = FALSE)
oa = Seurat::DimPlot(CRC_Qian_filtered, group.by='original_annotation', label=TRUE, label.size = 3) * Seurat::NoLegend()
s = Seurat::DimPlot(CRC_Qian_filtered, group.by='sample', label.size = 3) +
  ggplot2::theme(legend.text = ggplot2::element_text(size = 8))
p = Seurat::DimPlot(CRC_Qian_filtered, group.by='patient', label.size = 3) +
  ggplot2::theme(legend.text = ggplot2::element_text(size = 8))
st = Seurat::DimPlot(CRC_Qian_filtered, group.by='sample_type', label.size = 3) +
  ggplot2::theme(legend.text = ggplot2::element_text(size = 8))
phase = Seurat::DimPlot(CRC_Qian_filtered, group.by='Phase', label.size = 3) +
  ggplot2::theme(legend.text = ggplot2::element_text(size = 8))
mito = Seurat::FeaturePlot(CRC_Qian_filtered, features=c('percent.mitochondrial_RNA'))

(oa | s) / (p | st)
(oa) / (phase | mito)



# Some cells in the samples from healthy tissue were classified by these authors as cancer ?!?
normal_tissue_cells = rownames(CRC_Qian_filtered@meta.data)[CRC_Qian_filtered@meta.data$sample_type=='Normal Tissue']
Seurat::DimPlot(subset(CRC_Qian_filtered, cells=normal_tissue_cells),  group.by='original_annotation', label=TRUE, label.size = 3)
Seurat::FeaturePlot(subset(CRC_Qian_filtered, cells=normal_tissue_cells), features=c('MYC', 'AXIN2', 'RNF43', 'AREG')) # oncogenes
Seurat::FeaturePlot(subset(CRC_Qian_filtered, cells=normal_tissue_cells), features=c('THY1', 'FAP')) # CAFs
Seurat::FeaturePlot(subset(CRC_Qian_filtered, cells=normal_tissue_cells), features=c('COL1A1')) # Fibroblasts

non_normal_tissue_cells = rownames(CRC_Qian_filtered@meta.data)[CRC_Qian_filtered@meta.data$sample_type!='Normal Tissue']
Seurat::DimPlot(subset(CRC_Qian_filtered, cells=non_normal_tissue_cells), group.by='original_annotation', label=TRUE, label.size = 3)
Seurat::FeaturePlot(subset(CRC_Qian_filtered, cells=non_normal_tissue_cells), features=c('MYC', 'AXIN2', 'RNF43', 'AREG')) # oncogenes
Seurat::FeaturePlot(subset(CRC_Qian_filtered, cells=non_normal_tissue_cells), features=c('THY1', 'FAP')) # CAFs
Seurat::FeaturePlot(subset(CRC_Qian_filtered, cells=non_normal_tissue_cells), features=c('COL1A1')) # Fibroblasts
