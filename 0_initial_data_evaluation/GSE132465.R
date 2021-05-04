project_dir = '/home/scardoso/Documents/PhD/CRC_ATLAS'

GSE132465_data = paste(project_dir, 'original_datasets/GSE132465/data/data.txt', sep='/')
GSE132465_metadata = paste(project_dir, 'original_datasets/GSE132465/metadata.txt', sep='/')

# Cell-cycle markers:
cell_cycle_markers = paste(project_dir, 'utils/cycle_markers.rda', sep='/')
load(cell_cycle_markers)

# Read data to a seurat object:
data_matrix = bigreadr::big_fread2(GSE132465_data)
rownames(data_matrix) = data_matrix$Index
data_matrix = data_matrix[,-1]
GSE132465 = Seurat::CreateSeuratObject(counts=data_matrix, project='GSE132465')
remove(data_matrix)
gc()

# Read metadata file:
GSE132465_original_metadata = read.csv(GSE132465_metadata, row.names=1, header=TRUE)

# Update metadata:
colnames(GSE132465@meta.data) = c('sample', 'nUMI', 'nGene')
GSE132465[['patient']] = GSE132465_original_metadata[rownames(GSE132465@meta.data), 'Patient']
GSE132465[['sample_type']] = GSE132465_original_metadata[rownames(GSE132465@meta.data), 'Class']
GSE132465[['original_annotation_1']] = GSE132465_original_metadata[rownames(GSE132465@meta.data), 'Cell_type']
GSE132465[['original_annotation_2']] = GSE132465_original_metadata[rownames(GSE132465@meta.data), 'Cell_subtype']

# Store original data as seurat file:
SeuratDisk::SaveH5Seurat(GSE132465, paste(project_dir, 'original_datasets/GSE132465/original_data.h5Seurat', sep='/'))

# Percentage of mitochondrial RNA:
features = grep(pattern='^MT-', x=rownames(GSE132465@assays$RNA@counts), value = TRUE)
percent.featureset = colSums(as.matrix(Seurat::GetAssayData(GSE132465, assay='RNA', slot="counts")[features, , drop = FALSE]))/
  GSE132465$nUMI * 100
GSE132465[['percent.mitochondrial_RNA']] = percent.featureset

# Complexity
GSE132465[['log10.complexity']] = log10(GSE132465$nGene) / log10(GSE132465$nUMI)


# Quality control - Before
# Number of cell counts per sample
ggplot2::ggplot(GSE132465@meta.data, ggplot2::aes(x=sample, fill=sample)) + 
  ggplot2::geom_bar() + ggplot2::theme_classic() +
  ggplot2::theme(axis.text.x = ggplot2::element_text(angle = 45, vjust = 1, hjust=1)) +
  ggplot2::theme(plot.title = ggplot2::element_text(hjust=0.5, face="bold"), legend.position = "none") +
  ggplot2::ggtitle("NCells")
# UMI counts per cell
ggplot2::ggplot(GSE132465@meta.data, ggplot2::aes(x=nUMI, fill=sample)) + 
  ggplot2::geom_density(alpha = 0.2) + 
  ggplot2::scale_x_log10() + 
  ggplot2::theme_classic() +
  ggplot2::xlab('Nº UMIs') +
  ggplot2::ylab("Cell density") +
  ggplot2::geom_vline(xintercept = 1000) +
  ggplot2::ggtitle('UMI counts per cell')
# Number of detected genes per cell
ggplot2::ggplot(GSE132465@meta.data, ggplot2::aes(x=nGene, fill=sample)) + 
  ggplot2::geom_density(alpha = 0.2) + 
  ggplot2::theme_classic() + ggplot2::xlab('Nº Genes') + ggplot2::ylab('Density') +
  ggplot2::scale_x_log10() + 
  ggplot2::geom_vline(xintercept = 250) +
  ggplot2::ggtitle('Number of detected genes per cell')
# Distribution of genes detected per cell:
ggplot2::ggplot(GSE132465@meta.data, ggplot2::aes(y=log10(nGene), fill=sample)) + 
  ggplot2::geom_boxplot() + 
  ggplot2::theme_classic() +
  ggplot2::theme(plot.title = ggplot2::element_text(hjust=0.5, face="bold")) +
  ggplot2::ggtitle('Distribution of genes detected per cell:')
# UMIs vs Genes
ggplot2::ggplot(GSE132465@meta.data, ggplot2::aes(x=nUMI, y=nGene, color=percent.mitochondrial_RNA)) + 
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
ggplot2::ggplot(GSE132465@meta.data, ggplot2::aes(x=percent.mitochondrial_RNA, fill=sample)) + 
  ggplot2::geom_density(alpha = 0.2) + 
  ggplot2::scale_x_log10() + 
  ggplot2::theme_classic() +
  ggplot2::geom_vline(xintercept = 20) +
  ggplot2::ggtitle('Percentage of mitochondrial gene counts per cell')
# Complexity
ggplot2::ggplot(GSE132465@meta.data, ggplot2::aes(x=log10.complexity, fill=sample)) +
  ggplot2::geom_density(alpha = 0.2) +
  ggplot2::theme_classic() +
  ggplot2::geom_vline(xintercept = 0.8) +
  ggplot2::ggtitle('Complexity')

# Filter out low quality cells:
cells_keep = rownames(GSE132465@meta.data)[(GSE132465$nUMI >= 1000) &
                                             (GSE132465$nGene >= 250) &
                                             (GSE132465$nGene <= 6000) &
                                             (GSE132465$log10.complexity > 0.80) &
                                             (GSE132465$percent.mitochondrial_RNA <= 20)]
filtered_seurat = subset(GSE132465, cells=cells_keep)
# Filter out low-expressed genes:
counts = Seurat::GetAssayData(object=filtered_seurat, slot="counts")
nonzero <- counts > 0
keep_genes <- Matrix::rowSums(nonzero) >= 10
filtered_counts <- counts[keep_genes, ]
# Creation of filtered dataset:
GSE132465_filtered = Seurat::CreateSeuratObject(filtered_counts, meta.data=filtered_seurat@meta.data)

# Quality control - After
# Number of cell counts per sample
ggplot2::ggplot(GSE132465_filtered@meta.data, ggplot2::aes(x=sample, fill=sample)) + 
  ggplot2::geom_bar() + ggplot2::theme_classic() +
  ggplot2::theme(axis.text.x = ggplot2::element_text(angle = 45, vjust = 1, hjust=1)) +
  ggplot2::theme(plot.title = ggplot2::element_text(hjust=0.5, face="bold"), legend.position = "none") +
  ggplot2::ggtitle("NCells")
# UMI counts per cell
ggplot2::ggplot(GSE132465_filtered@meta.data, ggplot2::aes(x=nUMI, fill=sample)) + 
  ggplot2::geom_density(alpha = 0.2) + 
  ggplot2::scale_x_log10() + 
  ggplot2::theme_classic() +
  ggplot2::xlab('Nº UMIs') +
  ggplot2::ylab("Cell density") +
  ggplot2::geom_vline(xintercept = 1000) +
  ggplot2::ggtitle('UMI counts per cell')
# Number of detected genes per cell
ggplot2::ggplot(GSE132465_filtered@meta.data, ggplot2::aes(x=nGene, fill=sample)) + 
  ggplot2::geom_density(alpha = 0.2) + 
  ggplot2::theme_classic() + ggplot2::xlab('Nº Genes') + ggplot2::ylab('Density') +
  ggplot2::scale_x_log10() + 
  ggplot2::geom_vline(xintercept = 250) +
  ggplot2::ggtitle('Number of detected genes per cell')
# Distribution of genes detected per cell:
ggplot2::ggplot(GSE132465_filtered@meta.data, ggplot2::aes(y=log10(nGene), fill=sample)) + 
  ggplot2::geom_boxplot() + 
  ggplot2::theme_classic() +
  ggplot2::theme(plot.title = ggplot2::element_text(hjust=0.5, face="bold")) +
  ggplot2::ggtitle('Distribution of genes detected per cell:')
# UMIs vs Genes
ggplot2::ggplot(GSE132465_filtered@meta.data, ggplot2::aes(x=nUMI, y=nGene, color=percent.mitochondrial_RNA)) + 
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
ggplot2::ggplot(GSE132465_filtered@meta.data, ggplot2::aes(x=percent.mitochondrial_RNA, fill=sample)) + 
  ggplot2::geom_density(alpha = 0.2) + 
  ggplot2::scale_x_log10() + 
  ggplot2::theme_classic() +
  ggplot2::geom_vline(xintercept = 20) +
  ggplot2::ggtitle('Percentage of mitochondrial gene counts per cell')
# Complexity
ggplot2::ggplot(GSE132465_filtered@meta.data, ggplot2::aes(x=log10.complexity, fill=sample)) +
  ggplot2::geom_density(alpha = 0.2) +
  ggplot2::theme_classic() +
  ggplot2::geom_vline(xintercept = 0.8) +
  ggplot2::ggtitle('Complexity')


# Cell Cycle scoring
GSE132465_filtered = Seurat::NormalizeData(GSE132465_filtered, assay='RNA')
GSE132465_filtered = Seurat::CellCycleScoring(GSE132465_filtered, g2m.features=g2m_genes, s.features=s_genes, assay='RNA')
# Perform PCA to assess if cell cycle is a major source of variation
GSE132465_filtered = Seurat::FindVariableFeatures(GSE132465_filtered, selection.method="vst", nfeatures=2000, verbose=FALSE)
GSE132465_filtered = Seurat::ScaleData(GSE132465_filtered)
GSE132465_filtered = Seurat::RunPCA(GSE132465_filtered)
# Plot results:
Seurat::DimPlot(GSE132465_filtered, reduction="pca", group.by='Phase', pt.size=1)
Seurat::DimPlot(GSE132465_filtered, reduction="pca", group.by='Phase', split.by = 'Phase', pt.size=1)


# Visualize data in clusters:
GSE132465_filtered = Seurat::SCTransform(GSE132465_filtered, variable.features.n=2000, ncells=2000, conserve.memory = TRUE)
GSE132465_filtered = Seurat::RunPCA(GSE132465_filtered, verbose = FALSE)
GSE132465_filtered = Seurat::RunUMAP(GSE132465_filtered, dims = 1:30, verbose = FALSE)
oa1 = Seurat::DimPlot(GSE132465_filtered, group.by='original_annotation_1', label=TRUE, label.size = 3) +
  ggplot2::theme(legend.text = ggplot2::element_text(size = 8))
oa2 = Seurat::DimPlot(GSE132465_filtered, group.by='original_annotation_2', label=TRUE, label.size = 2, repel=TRUE) +
  ggplot2::theme(legend.text = ggplot2::element_text(size = 8))
p = Seurat::DimPlot(GSE132465_filtered, group.by='patient', label.size = 3) +
  ggplot2::theme(legend.text = ggplot2::element_text(size = 8))
s = Seurat::DimPlot(GSE132465_filtered, group.by='sample', label.size = 3) +
  ggplot2::theme(legend.text = ggplot2::element_text(size = 8))
st = Seurat::DimPlot(GSE132465_filtered, group.by='sample_type', label.size = 3) +
  ggplot2::theme(legend.text = ggplot2::element_text(size = 8))
phase = Seurat::DimPlot(GSE132465_filtered, group.by='Phase', label.size = 3) +
  ggplot2::theme(legend.text = ggplot2::element_text(size = 8))
mito = Seurat::FeaturePlot(GSE132465_filtered, features=c('percent.mitochondrial_RNA'))

(oa1 | oa2) / (p | st)
(oa1 | s) / (phase | mito)

normal_tissue_cells = rownames(GSE132465_filtered@meta.data)[GSE132465_filtered@meta.data$sample_type=='Normal']
Seurat::DimPlot(subset(GSE132465_filtered, cells=normal_tissue_cells),  group.by='original_annotation_2',
                label=TRUE, label.size=3, repel=TRUE)

non_normal_tissue_cells = rownames(GSE132465_filtered@meta.data)[GSE132465_filtered@meta.data$sample_type!='Normal']
Seurat::DimPlot(subset(GSE132465_filtered, cells=non_normal_tissue_cells),  group.by='original_annotation_2',
                label=TRUE, label.size=3, repel=TRUE)

