project_dir = '/home/scardoso/Documents/PhD/CRC_ATLAS'

GSE144735_data = paste(project_dir, 'original_datasets/GSE144735/data/data.txt', sep='/')
GSE144735_metadata = paste(project_dir, 'original_datasets/GSE144735/metadata.txt', sep='/')

# Cell-cycle markers:
cell_cycle_markers = paste(project_dir, 'utils/cycle_markers.rda', sep='/')
load(cell_cycle_markers)

# Read data to a seurat object:
data_matrix = bigreadr::big_fread2(GSE144735_data)
rownames(data_matrix) = data_matrix$Index
data_matrix = data_matrix[,-1]
GSE144735 = Seurat::CreateSeuratObject(counts=data_matrix, project='GSE144735')
remove(data_matrix)
gc()

# Read metadata file:
GSE144735_original_metadata = read.csv(GSE144735_metadata, row.names=1, header=TRUE)

# Update metadata:
colnames(GSE144735@meta.data) = c('sample', 'nUMI', 'nGene')
GSE144735[['patient']] = GSE144735_original_metadata[rownames(GSE144735@meta.data), 'Patient']
GSE144735[['sample_type']] = GSE144735_original_metadata[rownames(GSE144735@meta.data), 'Class']
GSE144735[['sample_type']][GSE144735[['sample_type']]=='Tumor'] = 'Tumour Core'
GSE144735[['sample_type']][GSE144735[['sample_type']]=='Border'] = 'Tumour Border'
GSE144735[['sample_type']][GSE144735[['sample_type']]=='Normal'] = 'Normal Tissue'
GSE144735[['original_annotation_1']] = GSE144735_original_metadata[rownames(GSE144735@meta.data), 'Cell_type']
GSE144735[['original_annotation_2']] = GSE144735_original_metadata[rownames(GSE144735@meta.data), 'Cell_subtype']

# Store original data as seurat file:
SeuratDisk::SaveH5Seurat(GSE144735, paste(project_dir, 'original_datasets/GSE144735/original_data.h5Seurat', sep='/'))

# Percentage of mitochondrial RNA:
features = grep(pattern='^MT-', x=rownames(GSE144735@assays$RNA@counts), value = TRUE)
percent.featureset = colSums(as.matrix(Seurat::GetAssayData(GSE144735, assay='RNA', slot="counts")[features, , drop = FALSE]))/
  GSE144735$nUMI * 100
GSE144735[['percent.mitochondrial_RNA']] = percent.featureset

# Complexity
GSE144735[['log10.complexity']] = log10(GSE144735$nGene) / log10(GSE144735$nUMI)


# Quality control - Before
# Number of cell counts per sample
ggplot2::ggplot(GSE144735@meta.data, ggplot2::aes(x=sample, fill=sample)) + 
  ggplot2::geom_bar() + ggplot2::theme_classic() +
  ggplot2::theme(axis.text.x = ggplot2::element_text(angle = 45, vjust = 1, hjust=1)) +
  ggplot2::theme(plot.title = ggplot2::element_text(hjust=0.5, face="bold"), legend.position = "none") +
  ggplot2::ggtitle("NCells")
# UMI counts per cell
ggplot2::ggplot(GSE144735@meta.data, ggplot2::aes(x=nUMI, fill=sample)) + 
  ggplot2::geom_density(alpha = 0.2) + 
  ggplot2::scale_x_log10() + 
  ggplot2::theme_classic() +
  ggplot2::xlab('Nº UMIs') +
  ggplot2::ylab("Cell density") +
  ggplot2::geom_vline(xintercept = 1000) +
  ggplot2::ggtitle('UMI counts per cell')
# Number of detected genes per cell
ggplot2::ggplot(GSE144735@meta.data, ggplot2::aes(x=nGene, fill=sample)) + 
  ggplot2::geom_density(alpha = 0.2) + 
  ggplot2::theme_classic() + ggplot2::xlab('Nº Genes') + ggplot2::ylab('Density') +
  ggplot2::scale_x_log10() + 
  ggplot2::geom_vline(xintercept = 250) +
  ggplot2::ggtitle('Number of detected genes per cell')
# Distribution of genes detected per cell:
ggplot2::ggplot(GSE144735@meta.data, ggplot2::aes(y=log10(nGene), fill=sample)) + 
  ggplot2::geom_boxplot() + 
  ggplot2::theme_classic() +
  ggplot2::theme(plot.title = ggplot2::element_text(hjust=0.5, face="bold")) +
  ggplot2::ggtitle('Distribution of genes detected per cell:')
# UMIs vs Genes
ggplot2::ggplot(GSE144735@meta.data, ggplot2::aes(x=nUMI, y=nGene, color=percent.mitochondrial_RNA)) + 
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
ggplot2::ggplot(GSE144735@meta.data, ggplot2::aes(x=percent.mitochondrial_RNA, fill=sample)) + 
  ggplot2::geom_density(alpha = 0.2) + 
  ggplot2::scale_x_log10() + 
  ggplot2::theme_classic() +
  ggplot2::geom_vline(xintercept = 20) +
  ggplot2::ggtitle('Percentage of mitochondrial gene counts per cell')
# Complexity
ggplot2::ggplot(GSE144735@meta.data, ggplot2::aes(x=log10.complexity, fill=sample)) +
  ggplot2::geom_density(alpha = 0.2) +
  ggplot2::theme_classic() +
  ggplot2::geom_vline(xintercept = 0.8) +
  ggplot2::ggtitle('Complexity')

# Filter out low quality cells:
cells_keep = rownames(GSE144735@meta.data)[(GSE144735$nUMI >= 1000) &
                                             (GSE144735$nGene >= 250) &
                                             (GSE144735$nGene <= 6000) &
                                             (GSE144735$log10.complexity > 0.80) &
                                             (GSE144735$percent.mitochondrial_RNA <= 20)]
filtered_seurat = subset(GSE144735, cells=cells_keep)
# Filter out low-expressed genes:
counts = Seurat::GetAssayData(object=filtered_seurat, slot="counts")
nonzero <- counts > 0
keep_genes <- Matrix::rowSums(nonzero) >= 10
filtered_counts <- counts[keep_genes, ]
# Creation of filtered dataset:
GSE144735_filtered = Seurat::CreateSeuratObject(filtered_counts, meta.data=filtered_seurat@meta.data)

# Quality control - After
# Number of cell counts per sample
ggplot2::ggplot(GSE144735_filtered@meta.data, ggplot2::aes(x=sample, fill=sample)) + 
  ggplot2::geom_bar() + ggplot2::theme_classic() +
  ggplot2::theme(axis.text.x = ggplot2::element_text(angle = 45, vjust = 1, hjust=1)) +
  ggplot2::theme(plot.title = ggplot2::element_text(hjust=0.5, face="bold"), legend.position = "none") +
  ggplot2::ggtitle("NCells")
# UMI counts per cell
ggplot2::ggplot(GSE144735_filtered@meta.data, ggplot2::aes(x=nUMI, fill=sample)) + 
  ggplot2::geom_density(alpha = 0.2) + 
  ggplot2::scale_x_log10() + 
  ggplot2::theme_classic() +
  ggplot2::xlab('Nº UMIs') +
  ggplot2::ylab("Cell density") +
  ggplot2::geom_vline(xintercept = 1000) +
  ggplot2::ggtitle('UMI counts per cell')
# Number of detected genes per cell
ggplot2::ggplot(GSE144735_filtered@meta.data, ggplot2::aes(x=nGene, fill=sample)) + 
  ggplot2::geom_density(alpha = 0.2) + 
  ggplot2::theme_classic() + ggplot2::xlab('Nº Genes') + ggplot2::ylab('Density') +
  ggplot2::scale_x_log10() + 
  ggplot2::geom_vline(xintercept = 250) +
  ggplot2::ggtitle('Number of detected genes per cell')
# Distribution of genes detected per cell:
ggplot2::ggplot(GSE144735_filtered@meta.data, ggplot2::aes(y=log10(nGene), fill=sample)) + 
  ggplot2::geom_boxplot() + 
  ggplot2::theme_classic() +
  ggplot2::theme(plot.title = ggplot2::element_text(hjust=0.5, face="bold")) +
  ggplot2::ggtitle('Distribution of genes detected per cell:')
# UMIs vs Genes
ggplot2::ggplot(GSE144735_filtered@meta.data, ggplot2::aes(x=nUMI, y=nGene, color=percent.mitochondrial_RNA)) + 
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
ggplot2::ggplot(GSE144735_filtered@meta.data, ggplot2::aes(x=percent.mitochondrial_RNA, fill=sample)) + 
  ggplot2::geom_density(alpha = 0.2) + 
  ggplot2::scale_x_log10() + 
  ggplot2::theme_classic() +
  ggplot2::geom_vline(xintercept = 20) +
  ggplot2::ggtitle('Percentage of mitochondrial gene counts per cell')
# Complexity
ggplot2::ggplot(GSE144735_filtered@meta.data, ggplot2::aes(x=log10.complexity, fill=sample)) +
  ggplot2::geom_density(alpha = 0.2) +
  ggplot2::theme_classic() +
  ggplot2::geom_vline(xintercept = 0.8) +
  ggplot2::ggtitle('Complexity')


# Cell Cycle scoring
GSE144735_filtered = Seurat::NormalizeData(GSE144735_filtered, assay='RNA')
GSE144735_filtered = Seurat::CellCycleScoring(GSE144735_filtered, g2m.features=g2m_genes, s.features=s_genes, assay='RNA')
# Perform PCA to assess if cell cycle is a major source of variation
GSE144735_filtered = Seurat::FindVariableFeatures(GSE144735_filtered, selection.method="vst", nfeatures=2000, verbose=FALSE)
GSE144735_filtered = Seurat::ScaleData(GSE144735_filtered)
GSE144735_filtered = Seurat::RunPCA(GSE144735_filtered)
# Plot results:
Seurat::DimPlot(GSE144735_filtered, reduction="pca", group.by='Phase', pt.size=1)
Seurat::DimPlot(GSE144735_filtered, reduction="pca", group.by='Phase', split.by = 'Phase', pt.size=1)


# Visualize data in clusters:
GSE144735_filtered = Seurat::SCTransform(GSE144735_filtered, variable.features.n=2000, ncells=2000, conserve.memory = TRUE)
GSE144735_filtered = Seurat::RunPCA(GSE144735_filtered, verbose = FALSE)
GSE144735_filtered = Seurat::RunUMAP(GSE144735_filtered, dims = 1:30, verbose = FALSE)
oa1 = Seurat::DimPlot(GSE144735_filtered, group.by='original_annotation_1', label=TRUE, label.size=2, repel=TRUE) +
  ggplot2::theme(legend.text = ggplot2::element_text(size = 8))
oa2 = Seurat::DimPlot(GSE144735_filtered, group.by='original_annotation_2', label=TRUE, label.size=2, repel=TRUE) +
  ggplot2::theme(legend.text = ggplot2::element_text(size = 8))
p = Seurat::DimPlot(GSE144735_filtered, group.by='patient', label.size = 3) +
  ggplot2::theme(legend.text = ggplot2::element_text(size = 8))
s = Seurat::DimPlot(GSE144735_filtered, group.by='sample', label.size = 3) +
  ggplot2::theme(legend.text = ggplot2::element_text(size = 8))
st = Seurat::DimPlot(GSE144735_filtered, group.by='sample_type', label.size = 3) +
  ggplot2::theme(legend.text = ggplot2::element_text(size = 8))
phase = Seurat::DimPlot(GSE144735_filtered, group.by='Phase', label.size = 3) +
  ggplot2::theme(legend.text = ggplot2::element_text(size = 8))
mito = Seurat::FeaturePlot(GSE144735_filtered, features=c('percent.mitochondrial_RNA'))

(oa1 | oa2) / (p | st)
(oa1 | s) / (phase | mito)

normal_tissue_cells = rownames(GSE144735_filtered@meta.data)[GSE144735_filtered@meta.data$sample_type=='Normal Tissue']
Seurat::DimPlot(subset(GSE144735_filtered, cells=normal_tissue_cells),  group.by='original_annotation_2',
                label=TRUE, label.size=3, repel=TRUE)

non_normal_tissue_cells = rownames(GSE144735_filtered@meta.data)[GSE144735_filtered@meta.data$sample_type!='Normal Tissue']
Seurat::DimPlot(subset(GSE144735_filtered, cells=non_normal_tissue_cells),  group.by='original_annotation_2',
                label=TRUE, label.size = 3, repel=TRUE)


