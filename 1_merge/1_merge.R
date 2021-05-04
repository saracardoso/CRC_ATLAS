project_dir = '/home/scardoso/Documents/PhD/CRC_ATLAS'

# Load datasets separately:
CRC_Qian = SeuratDisk::LoadH5Seurat(paste(project_dir, 'original_datasets/CRC_qian/original_data.h5Seurat', sep='/'))
colnames(CRC_Qian@meta.data)[6] = 'original_annotation_1'
CRC_Qian[['dataset']] = rep('CRC_Qian', dim(CRC_Qian@meta.data)[1])
GSE132465 = SeuratDisk::LoadH5Seurat(paste(project_dir, 'original_datasets/GSE132465/original_data.h5Seurat', sep='/'))
GSE132465[['dataset']] = rep('GSE132465', dim(GSE132465@meta.data)[1])
GSE144735 = SeuratDisk::LoadH5Seurat(paste(project_dir, 'original_datasets/GSE144735/original_data.h5Seurat', sep='/'))
GSE144735[['dataset']] = rep('GSE144735', dim(GSE144735@meta.data)[1])
Colon_smillie19 = SeuratDisk::LoadH5Seurat(paste(project_dir, 'original_datasets/Colon_smillie19/original_data.h5Seurat', sep='/'))
Colon_smillie19[['dataset']] = rep('Colon_smillie19', dim(Colon_smillie19@meta.data)[1])

# Merge datasets in one:
CRCatlas = merge(CRC_Qian, y = c(GSE132465, GSE144735, Colon_smillie19),
                 add.cell.ids = c("CRC_Qian", "GSE132465", "GSE144735", "Colon_smillie19"), project = "CRCatlas")

# Percentage of mitochondrial RNA:
features = grep(pattern='^MT-', x=rownames(CRCatlas@assays$RNA@counts), value = TRUE)
percent.featureset = colSums(as.matrix(Seurat::GetAssayData(CRCatlas, assay='RNA', slot="counts")[features, , drop = FALSE]))/
  CRCatlas$nUMI * 100
CRCatlas[['percent.mitochondrial_RNA']] = percent.featureset

# Complexity
CRCatlas[['log10.complexity']] = log10(CRCatlas$nGene) / log10(CRCatlas$nUMI)


# Quality control - Before
# Number of cell counts per sample
ggplot2::ggplot(CRCatlas@meta.data, ggplot2::aes(x=sample, fill=dataset)) + 
  ggplot2::geom_bar() + ggplot2::theme_classic() +
  ggplot2::theme(axis.text.x = ggplot2::element_text(angle = 45, vjust = 1, hjust=1)) +
  ggplot2::theme(plot.title = ggplot2::element_text(hjust=0.5, face="bold")) +
  ggplot2::ggtitle("NCells")
# UMI counts per cell
ggplot2::ggplot(CRCatlas@meta.data, ggplot2::aes(x=nUMI,fill=sample)) + 
  ggplot2::geom_density(alpha = 0.2) +ggplot2::scale_x_log10() + ggplot2::theme_classic() +
  ggplot2::theme(legend.position="none") +
  ggplot2::xlab('Nº UMIs') +
  ggplot2::ylab("Cell density") +
  ggplot2::geom_vline(xintercept = 1000) +
  ggplot2::ggtitle('UMI counts per cell')
# Number of detected genes per cell
ggplot2::ggplot(CRCatlas@meta.data, ggplot2::aes(x=nGene, fill=sample)) + 
  ggplot2::geom_density(alpha = 0.2) + 
  ggplot2::theme_classic() + ggplot2::xlab('Nº Genes') + ggplot2::ylab('Density') +
  ggplot2::scale_x_log10() + 
  ggplot2::geom_vline(xintercept = 250) +
  ggplot2::theme(legend.position="none") +
  ggplot2::ggtitle('Number of detected genes per cell')
# Distribution of genes detected per cell:
ggplot2::ggplot(CRCatlas@meta.data, ggplot2::aes(y=log10(nGene), x=sample, fill=dataset)) + 
  ggplot2::geom_boxplot() +  ggplot2::theme_classic() +
  ggplot2::theme(axis.text.x = ggplot2::element_text(angle = 45, vjust = 1, hjust=1)) +
  ggplot2::theme(plot.title = ggplot2::element_text(hjust=0.5, face="bold")) +
  ggplot2::ggtitle('Distribution of genes detected per cell:')
# UMIs vs Genes
ggplot2::ggplot(CRCatlas@meta.data, ggplot2::aes(x=nUMI, y=nGene, color=percent.mitochondrial_RNA)) + 
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
ggplot2::ggplot(CRCatlas@meta.data, ggplot2::aes(x=percent.mitochondrial_RNA, fill=sample)) + 
  ggplot2::geom_density(alpha = 0.2) + ggplot2::scale_x_log10() + ggplot2::theme_classic() +
  ggplot2::geom_vline(xintercept = 20) +
  ggplot2::theme(legend.position="none") +
  ggplot2::ggtitle('Percentage of mitochondrial gene counts per cell')
# Complexity
ggplot2::ggplot(CRCatlas@meta.data, ggplot2::aes(x=log10.complexity, fill=sample)) +
  ggplot2::geom_density(alpha = 0.2) + ggplot2::theme_classic() +
  ggplot2::geom_vline(xintercept = 0.8) +
  ggplot2::theme(legend.position="none") +
  ggplot2::ggtitle('Complexity')

# Filter out low quality cells:
cells_keep = rownames(CRCatlas@meta.data)[(CRCatlas$nUMI >= 1000) &
                                            (CRCatlas$nGene >= 250) &
                                            (CRCatlas$nGene <= 6000) &
                                            (CRCatlas$log10.complexity > 0.80) &
                                            (CRCatlas$percent.mitochondrial_RNA <= 20)]
filtered_seurat = subset(CRCatlas, cells=cells_keep)
# Filter out low-expressed genes:
counts = Seurat::GetAssayData(object=filtered_seurat, slot="counts")
nonzero <- counts > 0
keep_genes <- Matrix::rowSums(nonzero) >= 10
filtered_counts <- counts[keep_genes, ]
# Creation of filtered dataset:
CRCatlas_filtered = Seurat::CreateSeuratObject(filtered_counts, meta.data=filtered_seurat@meta.data)

# Quality control - After
# Number of cell counts per sample
ggplot2::ggplot(CRCatlas_filtered@meta.data, ggplot2::aes(x=sample, fill=dataset)) + 
  ggplot2::geom_bar() + ggplot2::theme_classic() +
  ggplot2::theme(axis.text.x = ggplot2::element_text(angle = 45, vjust = 1, hjust=1)) +
  ggplot2::theme(plot.title = ggplot2::element_text(hjust=0.5, face="bold")) +
  ggplot2::ggtitle("NCells")
# UMI counts per cell
ggplot2::ggplot(CRCatlas_filtered@meta.data, ggplot2::aes(x=nUMI,fill=sample)) + 
  ggplot2::geom_density(alpha = 0.2) +ggplot2::scale_x_log10() + ggplot2::theme_classic() +
  ggplot2::theme(legend.position="none") +
  ggplot2::xlab('Nº UMIs') +
  ggplot2::ylab("Cell density") +
  ggplot2::geom_vline(xintercept = 1000) +
  ggplot2::ggtitle('UMI counts per cell')
# Number of detected genes per cell
ggplot2::ggplot(CRCatlas_filtered@meta.data, ggplot2::aes(x=nGene, fill=sample)) + 
  ggplot2::geom_density(alpha = 0.2) + 
  ggplot2::theme_classic() + ggplot2::xlab('Nº Genes') + ggplot2::ylab('Density') +
  ggplot2::scale_x_log10() + 
  ggplot2::geom_vline(xintercept = 250) +
  ggplot2::theme(legend.position="none") +
  ggplot2::ggtitle('Number of detected genes per cell')
# Distribution of genes detected per cell:
ggplot2::ggplot(CRCatlas_filtered@meta.data, ggplot2::aes(y=log10(nGene), x=sample, fill=dataset)) + 
  ggplot2::geom_boxplot() + ggplot2::theme_classic() +
  ggplot2::theme(axis.text.x = ggplot2::element_text(angle = 45, vjust = 1, hjust=1)) +
  ggplot2::theme(plot.title = ggplot2::element_text(hjust=0.5, face="bold")) +
  ggplot2::ggtitle('Distribution of genes detected per cell:')
# UMIs vs Genes
ggplot2::ggplot(CRCatlas_filtered@meta.data, ggplot2::aes(x=nUMI, y=nGene, color=percent.mitochondrial_RNA)) + 
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
ggplot2::ggplot(CRCatlas_filtered@meta.data, ggplot2::aes(x=percent.mitochondrial_RNA, fill=sample)) + 
  ggplot2::geom_density(alpha = 0.2) + ggplot2::scale_x_log10() + ggplot2::theme_classic() +
  ggplot2::geom_vline(xintercept = 20) +
  ggplot2::theme(legend.position="none") +
  ggplot2::ggtitle('Percentage of mitochondrial gene counts per cell')
# Complexity
ggplot2::ggplot(CRCatlas_filtered@meta.data, ggplot2::aes(x=log10.complexity, fill=sample)) +
  ggplot2::geom_density(alpha = 0.2) + ggplot2::theme_classic() +
  ggplot2::geom_vline(xintercept = 0.8) +
  ggplot2::theme(legend.position="none") +
  ggplot2::ggtitle('Complexity')

# Store initial merge of CRCatlas:
SeuratDisk::SaveH5Seurat(CRCatlas_filtered, paste(project_dir, '1_merge/CRC.h5Seurat', sep='/'))
