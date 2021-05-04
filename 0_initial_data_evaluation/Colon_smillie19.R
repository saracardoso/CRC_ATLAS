project_dir = '/home/scardoso/Documents/PhD/CRC_ATLAS'

Colon_smillie19_data_dir = paste(project_dir, 'original_datasets/Colon_smillie19', sep='/')
Colon_smillie19_metadata = paste(project_dir, 'original_datasets/Colon_smillie19/metadata.txt', sep='/')

# Cell-cycle markers:
cell_cycle_markers = paste(project_dir, 'utils/cycle_markers.rda', sep='/')
load(cell_cycle_markers)

# Read data to a seurat object:
epithelial_seurat_data = Seurat::Read10X(data.dir=paste(Colon_smillie19_data_dir, 'epithelial', sep='/'), gene.column=1)
epithelial_Colon_smillie19 = Seurat::CreateSeuratObject(counts=epithelial_seurat_data, project='Colon_smillie19_epithelial')
gc()
stromal_seurat_data = Seurat::Read10X(data.dir=paste(Colon_smillie19_data_dir, 'stromal', sep='/'), gene.column=1)
stromal_Colon_smillie19 = Seurat::CreateSeuratObject(counts=stromal_seurat_data, project='Colon_smillie19_stromal')
gc()
immune_seurat_data = Seurat::Read10X(data.dir=paste(Colon_smillie19_data_dir, 'immune', sep='/'), gene.column=1)
immune_Colon_smillie19 = Seurat::CreateSeuratObject(counts=immune_seurat_data, project='Colon_smillie19_immune')
gc()
remove(epithelial_seurat_data, stromal_seurat_data, immune_seurat_data)
gc()

# Merge into one dataset:
Colon_smillie19 = merge(epithelial_Colon_smillie19, y=c(stromal_Colon_smillie19, immune_Colon_smillie19),
                        add.cell.ids = NULL, project = "Colon_smillie19")
remove(epithelial_Colon_smillie19, immune_Colon_smillie19, stromal_Colon_smillie19)
gc()

# Read metadata file:
Colon_smillie19_original_metadata = read.table(Colon_smillie19_metadata, header=TRUE, sep='\t')[-1,]
rownames(Colon_smillie19_original_metadata) = Colon_smillie19_original_metadata$NAME

# Get only healthy individuals (unhealthy individuals are not CRC patients):
healthy_cells = Colon_smillie19_original_metadata[Colon_smillie19_original_metadata$Health=='Healthy','NAME']
Colon_smillie19 = subset(Colon_smillie19, cells = healthy_cells)
gc()

# Update metadata:
Colon_smillie19@meta.data = Colon_smillie19@meta.data[,-1]
colnames(Colon_smillie19@meta.data) = c('nUMI', 'nGene')

Colon_smillie19[['sample']] = Colon_smillie19_original_metadata[rownames(Colon_smillie19@meta.data), 'Sample']
Colon_smillie19[['patient']] = Colon_smillie19_original_metadata[rownames(Colon_smillie19@meta.data), 'Subject']
Colon_smillie19[['sample_type']] = Colon_smillie19_original_metadata[rownames(Colon_smillie19@meta.data), 'Health']
Colon_smillie19[['original_annotation_1']] = Colon_smillie19_original_metadata[rownames(Colon_smillie19@meta.data), 'Cluster']
Colon_smillie19[['Location']] = Colon_smillie19_original_metadata[rownames(Colon_smillie19@meta.data), 'Location']

# Store original data as seurat file:
SeuratDisk::SaveH5Seurat(Colon_smillie19, paste(project_dir, 'original_datasets/Colon_smillie19/original_data.h5Seurat', sep='/'))

# Percentage of mitochondrial RNA:
features = grep(pattern='^MT-', x=rownames(Colon_smillie19@assays$RNA@counts), value = TRUE)
percent.featureset = colSums(as.matrix(Seurat::GetAssayData(Colon_smillie19, assay='RNA', slot="counts")[features, , drop = FALSE]))/
  Colon_smillie19$nUMI * 100
Colon_smillie19[['percent.mitochondrial_RNA']] = percent.featureset

# Complexity
Colon_smillie19[['log10.complexity']] = log10(Colon_smillie19$nGene) / log10(Colon_smillie19$nUMI)


# Quality control - Before
# Number of cell counts per sample
ggplot2::ggplot(Colon_smillie19@meta.data, ggplot2::aes(x=sample, fill=sample)) + 
  ggplot2::geom_bar() + ggplot2::theme_classic() +
  ggplot2::theme(axis.text.x = ggplot2::element_text(angle = 45, vjust = 1, hjust=1)) +
  ggplot2::theme(plot.title = ggplot2::element_text(hjust=0.5, face="bold"), legend.position = "none") +
  ggplot2::ggtitle("NCells")
# UMI counts per cell
ggplot2::ggplot(Colon_smillie19@meta.data, ggplot2::aes(x=nUMI, fill=sample)) + 
  ggplot2::geom_density(alpha = 0.2) + 
  ggplot2::scale_x_log10() + 
  ggplot2::theme_classic() +
  ggplot2::xlab('Nº UMIs') +
  ggplot2::ylab("Cell density") +
  ggplot2::geom_vline(xintercept = 1000) +
  ggplot2::ggtitle('UMI counts per cell')
# Number of detected genes per cell
ggplot2::ggplot(Colon_smillie19@meta.data, ggplot2::aes(x=nGene, fill=sample)) + 
  ggplot2::geom_density(alpha = 0.2) + 
  ggplot2::theme_classic() + ggplot2::xlab('Nº Genes') + ggplot2::ylab('Density') +
  ggplot2::scale_x_log10() + 
  ggplot2::geom_vline(xintercept = 250) +
  ggplot2::ggtitle('Number of detected genes per cell')
# Distribution of genes detected per cell:
ggplot2::ggplot(Colon_smillie19@meta.data, ggplot2::aes(y=log10(nGene), fill=sample)) + 
  ggplot2::geom_boxplot() + 
  ggplot2::theme_classic() +
  ggplot2::theme(plot.title = ggplot2::element_text(hjust=0.5, face="bold")) +
  ggplot2::ggtitle('Distribution of genes detected per cell:')
# UMIs vs Genes
ggplot2::ggplot(Colon_smillie19@meta.data, ggplot2::aes(x=nUMI, y=nGene, color=percent.mitochondrial_RNA)) + 
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
ggplot2::ggplot(Colon_smillie19@meta.data, ggplot2::aes(x=percent.mitochondrial_RNA, fill=sample)) + 
  ggplot2::geom_density(alpha = 0.2) + 
  ggplot2::scale_x_log10() + 
  ggplot2::theme_classic() +
  ggplot2::geom_vline(xintercept = 20) +
  ggplot2::ggtitle('Percentage of mitochondrial gene counts per cell')
# Complexity
ggplot2::ggplot(Colon_smillie19@meta.data, ggplot2::aes(x=log10.complexity, fill=sample)) +
  ggplot2::geom_density(alpha = 0.2) +
  ggplot2::theme_classic() +
  ggplot2::geom_vline(xintercept = 0.8) +
  ggplot2::ggtitle('Complexity')

# Filter out low quality cells:
cells_keep = rownames(Colon_smillie19@meta.data)[(Colon_smillie19$nUMI >= 1000) &
                                                   (Colon_smillie19$nGene >= 250) &
                                                   (Colon_smillie19$nGene <= 6000) &
                                                   (Colon_smillie19$log10.complexity > 0.80) &
                                                   (Colon_smillie19$percent.mitochondrial_RNA <= 20)]
filtered_seurat = subset(Colon_smillie19, cells=cells_keep)
# Filter out low-expressed genes:
counts = Seurat::GetAssayData(object=filtered_seurat, slot="counts")
nonzero <- counts > 0
keep_genes <- Matrix::rowSums(nonzero) >= 10
filtered_counts <- counts[keep_genes, ]
# Creation of filtered dataset:
Colon_smillie19_filtered = Seurat::CreateSeuratObject(filtered_counts, meta.data=filtered_seurat@meta.data)

# Quality control - After
# Number of cell counts per sample
ggplot2::ggplot(Colon_smillie19_filtered@meta.data, ggplot2::aes(x=sample, fill=sample)) + 
  ggplot2::geom_bar() + ggplot2::theme_classic() +
  ggplot2::theme(axis.text.x = ggplot2::element_text(angle = 45, vjust = 1, hjust=1)) +
  ggplot2::theme(plot.title = ggplot2::element_text(hjust=0.5, face="bold"), legend.position = "none") +
  ggplot2::ggtitle("NCells")
# UMI counts per cell
ggplot2::ggplot(Colon_smillie19_filtered@meta.data, ggplot2::aes(x=nUMI, fill=sample)) + 
  ggplot2::geom_density(alpha = 0.2) + 
  ggplot2::scale_x_log10() + 
  ggplot2::theme_classic() +
  ggplot2::xlab('Nº UMIs') +
  ggplot2::ylab("Cell density") +
  ggplot2::geom_vline(xintercept = 1000) +
  ggplot2::ggtitle('UMI counts per cell')
# Number of detected genes per cell
ggplot2::ggplot(Colon_smillie19_filtered@meta.data, ggplot2::aes(x=nGene, fill=sample)) + 
  ggplot2::geom_density(alpha = 0.2) + 
  ggplot2::theme_classic() + ggplot2::xlab('Nº Genes') + ggplot2::ylab('Density') +
  ggplot2::scale_x_log10() + 
  ggplot2::geom_vline(xintercept = 250) +
  ggplot2::ggtitle('Number of detected genes per cell')
# Distribution of genes detected per cell:
ggplot2::ggplot(Colon_smillie19_filtered@meta.data, ggplot2::aes(y=log10(nGene), fill=sample)) + 
  ggplot2::geom_boxplot() + 
  ggplot2::theme_classic() +
  ggplot2::theme(plot.title = ggplot2::element_text(hjust=0.5, face="bold")) +
  ggplot2::ggtitle('Distribution of genes detected per cell:')
# UMIs vs Genes
ggplot2::ggplot(Colon_smillie19_filtered@meta.data, ggplot2::aes(x=nUMI, y=nGene, color=percent.mitochondrial_RNA)) + 
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
ggplot2::ggplot(Colon_smillie19_filtered@meta.data, ggplot2::aes(x=percent.mitochondrial_RNA, fill=sample)) + 
  ggplot2::geom_density(alpha = 0.2) + 
  ggplot2::scale_x_log10() + 
  ggplot2::theme_classic() +
  ggplot2::geom_vline(xintercept = 20) +
  ggplot2::ggtitle('Percentage of mitochondrial gene counts per cell')
# Complexity
ggplot2::ggplot(Colon_smillie19_filtered@meta.data, ggplot2::aes(x=log10.complexity, fill=sample)) +
  ggplot2::geom_density(alpha = 0.2) +
  ggplot2::theme_classic() +
  ggplot2::geom_vline(xintercept = 0.8) +
  ggplot2::ggtitle('Complexity')


# Cell Cycle scoring
Colon_smillie19_filtered = Seurat::NormalizeData(Colon_smillie19_filtered, assay='RNA')
Colon_smillie19_filtered = Seurat::CellCycleScoring(Colon_smillie19_filtered, g2m.features=g2m_genes, s.features=s_genes, assay='RNA')
# Perform PCA to assess if cell cycle is a major source of variation
Colon_smillie19_filtered = Seurat::FindVariableFeatures(Colon_smillie19_filtered, selection.method="vst", nfeatures=2000, verbose=FALSE)
Colon_smillie19_filtered = Seurat::ScaleData(Colon_smillie19_filtered)
Colon_smillie19_filtered = Seurat::RunPCA(Colon_smillie19_filtered)
# Plot results:
Seurat::DimPlot(Colon_smillie19_filtered, reduction="pca", group.by='Phase', pt.size=1)
Seurat::DimPlot(Colon_smillie19_filtered, reduction="pca", group.by='Phase', split.by = 'Phase', pt.size=1)


# Visualize data in clusters:
Colon_smillie19_filtered = Seurat::SCTransform(Colon_smillie19_filtered, variable.features.n=2000, ncells=2000, conserve.memory = TRUE)
Colon_smillie19_filtered = Seurat::RunPCA(Colon_smillie19_filtered, verbose = FALSE)
Colon_smillie19_filtered = Seurat::RunUMAP(Colon_smillie19_filtered, dims = 1:30, verbose = FALSE)
oa = Seurat::DimPlot(Colon_smillie19_filtered, group.by='original_annotation', label=TRUE, label.size = 3, repel=TRUE) * Seurat::NoLegend()
s = Seurat::DimPlot(Colon_smillie19_filtered, group.by='sample', label.size = 3) +
  ggplot2::theme(legend.text = ggplot2::element_text(size = 8))
p = Seurat::DimPlot(Colon_smillie19_filtered, group.by='patient', label.size = 3) +
  ggplot2::theme(legend.text = ggplot2::element_text(size = 8))
phase = Seurat::DimPlot(Colon_smillie19_filtered, group.by='Phase', label.size = 3) +
  ggplot2::theme(legend.text = ggplot2::element_text(size = 8))
mito = Seurat::FeaturePlot(Colon_smillie19_filtered, features=c('percent.mitochondrial_RNA'))

oa / (p | s)
(oa | s) / (phase | mito)

# Gene expression:
Seurat::FeaturePlot(Colon_smillie19_filtered, features=c('CD3D', 'CD79A', 'CD14', 'COL1A1', 'ENG', 'EPCAM', 'TPSAB1', 'S100B'))
