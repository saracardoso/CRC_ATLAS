project_dir = '/home/scardoso/Documents/PhD/CRC_ATLAS'


# ########################################### #
# #### INITTIAL VISUALIZATION OF DATASET #### #
# ########################################### #

# Load subset of epithelial cells:
Epithelial = SeuratDisk::LoadH5Seurat(paste(project_dir, '2_annotation/results_Epithelial/datasets/epithelial_cells.h5Seurat', sep='/'))
invisible(gc())





# -----
# - UMAP visualizations
# -----


# 1. Run UMAP:
Epithelial = Seurat::RunUMAP(Epithelial, dims=1:50, reduction='pca')#elbow, reduction='pca')


# 2. Save object:
SeuratDisk::SaveH5Seurat(Epithelial, paste(project_dir, '2_annotation/results_Epithelial/datasets/Epithelial2.h5Seurat', sep='/'))


# 3. Metadata:
state = Seurat::DimPlot(Epithelial, reduction="umap", group.by='state', label=FALSE) + ggplot2::theme_minimal()
patients = Seurat::DimPlot(Epithelial, reduction="umap", group.by='patient', label=TRUE, label.size=3, pt.size=.1) +
  ggplot2::theme_minimal() + Seurat::NoLegend()
datasets = Seurat::DimPlot(Epithelial, reduction="umap", group.by='dataset', label=FALSE, pt.size=.1) +
  ggplot2::theme_minimal()
(state / datasets) | patients

