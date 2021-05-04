project_dir = '/home/scardoso/Documents/PhD/CRC_ATLAS'

# Load integrated CRC atlas:
CRCatlas_integrated = SeuratDisk::LoadH5Seurat(paste(project_dir, '1_merge/CRC_integrated.h5Seurat', sep='/'))


# Integrated dataset visualization:
CRCatlas_integrated = Seurat::ScaleData(CRCatlas_integrated, verbose = FALSE)
gc()
CRCatlas_integrated = Seurat::RunPCA(CRCatlas_integrated, verbose = TRUE)
gc()
CRCatlas_integrated = Seurat::RunUMAP(CRCatlas_integrated, reduction = "pca", dims = 1:50)
gc()

# Vector with colors of cell-types:
colors_vec = c('#93AA00', '#00C19F', '#FF61C3', '#00B9E3', '#00B9E3', '#00C19F', '#F8766D', '#00B9E3', '#93AA00', '#93AA00', '#00B9E3',
               '#00C19F', '#FF61C3', '#F8766D', '#00C19F')
names(colors_vec) = c("Cancer", "Myeloid", "T_cell", "Fibroblast", "EC", "Mast_cell", "B_cell", "Enteric_glia", "Epithelial",
                      "Epithelial cells", "Stromal cells", "Myeloids", "T cells", "B cells", "Mast cells")

# Visualize data in clusters:
dt =  Seurat::DimPlot(CRCatlas_integrated, group.by='dataset', label=TRUE, label.size = 3) +
  ggplot2::theme(legend.text = ggplot2::element_text(size = 8))
oa1 = Seurat::DimPlot(CRCatlas_integrated, group.by='original_annotation_1', label=TRUE, label.size = 3) +
  ggplot2::theme(legend.text = ggplot2::element_text(size = 8))
p = Seurat::DimPlot(CRCatlas_integrated, group.by='patient', label.size = 3) +
  ggplot2::theme(legend.text = ggplot2::element_text(size = 8))
st = Seurat::DimPlot(CRCatlas_integrated, group.by='sample_type', label.size = 3) +
  ggplot2::theme(legend.text = ggplot2::element_text(size = 8))
(dt | oa1) / (p | st)

Seurat::DimPlot(CRCatlas_integrated, group.by='sample_type', split.by='dataset',
                label=TRUE, label.size=3, repel=TRUE) +
  Seurat::NoLegend() + ggplot2::theme(legend.text=ggplot2::element_text(size=5))

# Some gene signatures will be looked into:
dim_dt_ann = Seurat::DimPlot(CRCatlas_integrated, group.by='original_annotation_1', split.by='dataset',
                             label=TRUE, label.size=3, repel=TRUE) +
  Seurat::NoLegend() + ggplot2::theme(legend.text=ggplot2::element_text(size=5))
CRCatlas_integrated = Seurat::SetIdent(CRCatlas_integrated, value='original_annotation_1')
Qian_cells = rownames(CRCatlas_integrated@meta.data)[CRCatlas_integrated@meta.data$dataset=='CRC_Qian']
GSE132465_cells = rownames(CRCatlas_integrated@meta.data)[CRCatlas_integrated@meta.data$dataset=='GSE132465']
GSE144735_cells = rownames(CRCatlas_integrated@meta.data)[CRCatlas_integrated@meta.data$dataset=='GSE144735']

tcell_1 = Seurat::FeaturePlot(CRCatlas_integrated, features=c('rna_CD3D'), cells=Qian_cells) +
  Seurat::NoLegend() + ggplot2::xlab('') +  ggplot2::ggtitle('')
tcell_2 = Seurat::FeaturePlot(CRCatlas_integrated, features=c('rna_CD3D'), cells=GSE132465_cells) +
  Seurat::NoLegend() + ggplot2::ylab('') +  ggplot2::ggtitle('') +
  ggplot2::theme(axis.ticks.y=ggplot2::element_blank(), axis.text.y=ggplot2::element_blank(), axis.line.y=ggplot2::element_blank())
tcell_3 = Seurat::FeaturePlot(CRCatlas_integrated, features=c('rna_CD3D'), cells=GSE144735_cells) +
  Seurat::NoLegend() + ggplot2::xlab('') + ggplot2::ylab('') +  ggplot2::ggtitle('') +
  ggplot2::theme(axis.ticks.y=ggplot2::element_blank(), axis.text.y=ggplot2::element_blank(), axis.line.y=ggplot2::element_blank())
(dim_dt_ann + ggplot2::ggtitle('CD3D')) / (tcell_1 | tcell_2 | tcell_3)

bcell_1 = Seurat::FeaturePlot(CRCatlas_integrated, features=c('rna_MS4A1'), cells=Qian_cells) +
  Seurat::NoLegend() + ggplot2::xlab('') +  ggplot2::ggtitle('')
bcell_2 = Seurat::FeaturePlot(CRCatlas_integrated, features=c('rna_MS4A1'), cells=GSE132465_cells) +
  Seurat::NoLegend() + ggplot2::ylab('') +  ggplot2::ggtitle('') +
  ggplot2::theme(axis.ticks.y=ggplot2::element_blank(), axis.text.y=ggplot2::element_blank(), axis.line.y=ggplot2::element_blank())
bcell_3 = Seurat::FeaturePlot(CRCatlas_integrated, features=c('rna_MS4A1'), cells=GSE144735_cells) +
  Seurat::NoLegend() + ggplot2::xlab('') + ggplot2::ylab('') +  ggplot2::ggtitle('') +
  ggplot2::theme(axis.ticks.y=ggplot2::element_blank(), axis.text.y=ggplot2::element_blank(), axis.line.y=ggplot2::element_blank())
(dim_dt_ann + ggplot2::ggtitle('MS4A1')) / (bcell_1 | bcell_2 | bcell_3)

plasmacell_1 = Seurat::FeaturePlot(CRCatlas_integrated, features=c('rna_CD79A'), cells=Qian_cells) +
  Seurat::NoLegend() + ggplot2::xlab('') +  ggplot2::ggtitle('')
plasmacell_2 = Seurat::FeaturePlot(CRCatlas_integrated, features=c('rna_CD79A'), cells=GSE132465_cells) +
  Seurat::NoLegend() + ggplot2::ylab('') +  ggplot2::ggtitle('') +
  ggplot2::theme(axis.ticks.y=ggplot2::element_blank(), axis.text.y=ggplot2::element_blank(), axis.line.y=ggplot2::element_blank())
plasmacell_3 = Seurat::FeaturePlot(CRCatlas_integrated, features=c('rna_CD79A'), cells=GSE144735_cells) +
  Seurat::NoLegend() + ggplot2::xlab('') + ggplot2::ylab('') +  ggplot2::ggtitle('') +
  ggplot2::theme(axis.ticks.y=ggplot2::element_blank(), axis.text.y=ggplot2::element_blank(), axis.line.y=ggplot2::element_blank())
(dim_dt_ann + ggplot2::ggtitle('CD79A')) / (plasmacell_1 | plasmacell_2 | plasmacell_3)

mono_1 = Seurat::FeaturePlot(CRCatlas_integrated, features=c('rna_CD14'), cells=Qian_cells) +
  Seurat::NoLegend() + ggplot2::xlab('') +  ggplot2::ggtitle('')
mono_2 = Seurat::FeaturePlot(CRCatlas_integrated, features=c('rna_CD14'), cells=GSE132465_cells) +
  Seurat::NoLegend() + ggplot2::ylab('') +  ggplot2::ggtitle('') +
  ggplot2::theme(axis.ticks.y=ggplot2::element_blank(), axis.text.y=ggplot2::element_blank(), axis.line.y=ggplot2::element_blank())
mono_3 = Seurat::FeaturePlot(CRCatlas_integrated, features=c('rna_CD14'), cells=GSE144735_cells) +
  Seurat::NoLegend() + ggplot2::xlab('') + ggplot2::ylab('') +  ggplot2::ggtitle('') +
  ggplot2::theme(axis.ticks.y=ggplot2::element_blank(), axis.text.y=ggplot2::element_blank(), axis.line.y=ggplot2::element_blank())
(dim_dt_ann + ggplot2::ggtitle('CD14')) / (mono_1 | mono_2 | mono_3)

epithelial_1 = Seurat::FeaturePlot(CRCatlas_integrated, features=c('rna_EPCAM'), cells=Qian_cells) +
  Seurat::NoLegend() + ggplot2::xlab('') +  ggplot2::ggtitle('')
epithelial_2 = Seurat::FeaturePlot(CRCatlas_integrated, features=c('rna_EPCAM'), cells=GSE132465_cells) +
  Seurat::NoLegend() + ggplot2::ylab('') +  ggplot2::ggtitle('') +
  ggplot2::theme(axis.ticks.y=ggplot2::element_blank(), axis.text.y=ggplot2::element_blank(), axis.line.y=ggplot2::element_blank())
epithelial_3 = Seurat::FeaturePlot(CRCatlas_integrated, features=c('rna_EPCAM'), cells=GSE144735_cells) +
  Seurat::NoLegend() + ggplot2::xlab('') + ggplot2::ylab('') +  ggplot2::ggtitle('') +
  ggplot2::theme(axis.ticks.y=ggplot2::element_blank(), axis.text.y=ggplot2::element_blank(), axis.line.y=ggplot2::element_blank())
(dim_dt_ann + ggplot2::ggtitle('EPCAM')) / (epithelial_1 | epithelial_2 | epithelial_3)

stromal_1 = Seurat::FeaturePlot(CRCatlas_integrated, features=c('rna_COL1A1'), cells=Qian_cells) +
  Seurat::NoLegend() + ggplot2::xlab('') +  ggplot2::ggtitle('')
stromal_2 = Seurat::FeaturePlot(CRCatlas_integrated, features=c('rna_COL1A1'), cells=GSE132465_cells) +
  Seurat::NoLegend() + ggplot2::ylab('') +  ggplot2::ggtitle('') +
  ggplot2::theme(axis.ticks.y=ggplot2::element_blank(), axis.text.y=ggplot2::element_blank(), axis.line.y=ggplot2::element_blank())
stromal_3 = Seurat::FeaturePlot(CRCatlas_integrated, features=c('rna_COL1A1'), cells=GSE144735_cells) +
  Seurat::NoLegend() + ggplot2::xlab('') + ggplot2::ylab('') +  ggplot2::ggtitle('') +
  ggplot2::theme(axis.ticks.y=ggplot2::element_blank(), axis.text.y=ggplot2::element_blank(), axis.line.y=ggplot2::element_blank())
(dim_dt_ann + ggplot2::ggtitle('COL1A1')) / (stromal_1 | stromal_2 | stromal_3)

endothelial_1 = Seurat::FeaturePlot(CRCatlas_integrated, features=c('rna_ENG'), cells=Qian_cells) +
  Seurat::NoLegend() + ggplot2::xlab('') +  ggplot2::ggtitle('')
endothelial_2 = Seurat::FeaturePlot(CRCatlas_integrated, features=c('rna_ENG'), cells=GSE132465_cells) +
  Seurat::NoLegend() + ggplot2::ylab('') +  ggplot2::ggtitle('') +
  ggplot2::theme(axis.ticks.y=ggplot2::element_blank(), axis.text.y=ggplot2::element_blank(), axis.line.y=ggplot2::element_blank())
endothelial_3 = Seurat::FeaturePlot(CRCatlas_integrated, features=c('rna_ENG'), cells=GSE144735_cells) +
  Seurat::NoLegend() + ggplot2::xlab('') + ggplot2::ylab('') +  ggplot2::ggtitle('') +
  ggplot2::theme(axis.ticks.y=ggplot2::element_blank(), axis.text.y=ggplot2::element_blank(), axis.line.y=ggplot2::element_blank())
(dim_dt_ann + ggplot2::ggtitle('ENG')) / (endothelial_1 | endothelial_2 | endothelial_3)

mast_1 = Seurat::FeaturePlot(CRCatlas_integrated, features=c('rna_KIT'), cells=Qian_cells) +
  Seurat::NoLegend() + ggplot2::xlab('') +  ggplot2::ggtitle('')
mast_2 = Seurat::FeaturePlot(CRCatlas_integrated, features=c('rna_KIT'), cells=GSE132465_cells) +
  Seurat::NoLegend() + ggplot2::ylab('') +  ggplot2::ggtitle('') +
  ggplot2::theme(axis.ticks.y=ggplot2::element_blank(), axis.text.y=ggplot2::element_blank(), axis.line.y=ggplot2::element_blank())
mast_3 = Seurat::FeaturePlot(CRCatlas_integrated, features=c('rna_KIT'), cells=GSE144735_cells) +
  Seurat::NoLegend() + ggplot2::xlab('') + ggplot2::ylab('') +  ggplot2::ggtitle('') +
  ggplot2::theme(axis.ticks.y=ggplot2::element_blank(), axis.text.y=ggplot2::element_blank(), axis.line.y=ggplot2::element_blank())
(dim_dt_ann + ggplot2::ggtitle('KIT')) / (mast_1 | mast_2 | mast_3)

entericglia_1 = Seurat::FeaturePlot(CRCatlas_integrated, features=c('rna_S100B'), cells=Qian_cells) +
  Seurat::NoLegend() + ggplot2::xlab('') +  ggplot2::ggtitle('')
entericglia_2 = Seurat::FeaturePlot(CRCatlas_integrated, features=c('rna_S100B'), cells=GSE132465_cells) +
  Seurat::NoLegend() + ggplot2::ylab('') +  ggplot2::ggtitle('') +
  ggplot2::theme(axis.ticks.y=ggplot2::element_blank(), axis.text.y=ggplot2::element_blank(), axis.line.y=ggplot2::element_blank())
entericglia_3 = Seurat::FeaturePlot(CRCatlas_integrated, features=c('rna_S100B'), cells=GSE144735_cells) +
  Seurat::NoLegend() + ggplot2::xlab('') + ggplot2::ylab('') +  ggplot2::ggtitle('') +
  ggplot2::theme(axis.ticks.y=ggplot2::element_blank(), axis.text.y=ggplot2::element_blank(), axis.line.y=ggplot2::element_blank())
(dim_dt_ann + ggplot2::ggtitle('S100B')) / (entericglia_1 | entericglia_2 | entericglia_3)

mdc_1 = Seurat::FeaturePlot(CRCatlas_integrated, features=c('rna_FCER1A'), cells=Qian_cells) +
  Seurat::NoLegend() + ggplot2::xlab('') +  ggplot2::ggtitle('')
mdc_2 = Seurat::FeaturePlot(CRCatlas_integrated, features=c('rna_FCER1A'), cells=GSE132465_cells) +
  Seurat::NoLegend() + ggplot2::ylab('') +  ggplot2::ggtitle('') +
  ggplot2::theme(axis.ticks.y=ggplot2::element_blank(), axis.text.y=ggplot2::element_blank(), axis.line.y=ggplot2::element_blank())
mdc_3 = Seurat::FeaturePlot(CRCatlas_integrated, features=c('rna_FCER1A'), cells=GSE144735_cells) +
  Seurat::NoLegend() + ggplot2::xlab('') + ggplot2::ylab('') +  ggplot2::ggtitle('') +
  ggplot2::theme(axis.ticks.y=ggplot2::element_blank(), axis.text.y=ggplot2::element_blank(), axis.line.y=ggplot2::element_blank())
(dim_dt_ann + ggplot2::ggtitle('FCER1A')) / (mdc_1 | mdc_2 | mdc_3)

pdc_1 = Seurat::FeaturePlot(CRCatlas_integrated, features=c('rna_GZMB'), cells=Qian_cells) +
  Seurat::NoLegend() + ggplot2::xlab('') +  ggplot2::ggtitle('')
pdc_2 = Seurat::FeaturePlot(CRCatlas_integrated, features=c('rna_GZMB'), cells=GSE132465_cells) +
  Seurat::NoLegend() + ggplot2::ylab('') +  ggplot2::ggtitle('') +
  ggplot2::theme(axis.ticks.y=ggplot2::element_blank(), axis.text.y=ggplot2::element_blank(), axis.line.y=ggplot2::element_blank())
pdc_3 = Seurat::FeaturePlot(CRCatlas_integrated, features=c('rna_GZMB'), cells=GSE144735_cells) +
  Seurat::NoLegend() + ggplot2::xlab('') + ggplot2::ylab('') +  ggplot2::ggtitle('') +
  ggplot2::theme(axis.ticks.y=ggplot2::element_blank(), axis.text.y=ggplot2::element_blank(), axis.line.y=ggplot2::element_blank())
(dim_dt_ann + ggplot2::ggtitle('GZMB')) / (pdc_1 | pdc_2 | pdc_3)

cd8_1 = Seurat::FeaturePlot(CRCatlas_integrated, features=c('rna_CD8A'), cells=Qian_cells) +
  Seurat::NoLegend() + ggplot2::xlab('') +  ggplot2::ggtitle('')
cd8_2 = Seurat::FeaturePlot(CRCatlas_integrated, features=c('rna_CD8A'), cells=GSE132465_cells) +
  Seurat::NoLegend() + ggplot2::ylab('') +  ggplot2::ggtitle('') +
  ggplot2::theme(axis.ticks.y=ggplot2::element_blank(), axis.text.y=ggplot2::element_blank(), axis.line.y=ggplot2::element_blank())
cd8_3 = Seurat::FeaturePlot(CRCatlas_integrated, features=c('rna_CD8A'), cells=GSE144735_cells) +
  Seurat::NoLegend() + ggplot2::xlab('') + ggplot2::ylab('') +  ggplot2::ggtitle('') +
  ggplot2::theme(axis.ticks.y=ggplot2::element_blank(), axis.text.y=ggplot2::element_blank(), axis.line.y=ggplot2::element_blank())
(dim_dt_ann + ggplot2::ggtitle('CD8A')) / (cd8_1 | cd8_2 | cd8_3)

cd4_1 = Seurat::FeaturePlot(CRCatlas_integrated, features=c('rna_IL7R'), cells=Qian_cells) +
  Seurat::NoLegend() + ggplot2::xlab('') +  ggplot2::ggtitle('')
cd4_2 = Seurat::FeaturePlot(CRCatlas_integrated, features=c('rna_IL7R'), cells=GSE132465_cells) +
  Seurat::NoLegend() + ggplot2::ylab('') +  ggplot2::ggtitle('') +
  ggplot2::theme(axis.ticks.y=ggplot2::element_blank(), axis.text.y=ggplot2::element_blank(), axis.line.y=ggplot2::element_blank())
cd4_3 = Seurat::FeaturePlot(CRCatlas_integrated, features=c('rna_IL7R'), cells=GSE144735_cells) +
  Seurat::NoLegend() + ggplot2::xlab('') + ggplot2::ylab('') +  ggplot2::ggtitle('') +
  ggplot2::theme(axis.ticks.y=ggplot2::element_blank(), axis.text.y=ggplot2::element_blank(), axis.line.y=ggplot2::element_blank())
(dim_dt_ann + ggplot2::ggtitle('IL7R')) / (cd4_1 | cd4_2 | cd4_3)

macro_1 = Seurat::FeaturePlot(CRCatlas_integrated, features=c('rna_MARCO'), cells=Qian_cells) +
  Seurat::NoLegend() + ggplot2::xlab('') +  ggplot2::ggtitle('')
macro_2 = Seurat::FeaturePlot(CRCatlas_integrated, features=c('rna_MARCO'), cells=GSE132465_cells) +
  Seurat::NoLegend() + ggplot2::ylab('') +  ggplot2::ggtitle('') +
  ggplot2::theme(axis.ticks.y=ggplot2::element_blank(), axis.text.y=ggplot2::element_blank(), axis.line.y=ggplot2::element_blank())
macro_3 = Seurat::FeaturePlot(CRCatlas_integrated, features=c('rna_MARCO'), cells=GSE144735_cells) +
  Seurat::NoLegend() + ggplot2::xlab('') + ggplot2::ylab('') +  ggplot2::ggtitle('') +
  ggplot2::theme(axis.ticks.y=ggplot2::element_blank(), axis.text.y=ggplot2::element_blank(), axis.line.y=ggplot2::element_blank())
(dim_dt_ann + ggplot2::ggtitle('MARCO')) / (macro_1 | macro_2 | macro_3)

nk_1 = Seurat::FeaturePlot(CRCatlas_integrated, features=c('rna_GNLY'), cells=Qian_cells) +
  Seurat::NoLegend() + ggplot2::xlab('') +  ggplot2::ggtitle('')
nk_2 = Seurat::FeaturePlot(CRCatlas_integrated, features=c('rna_GNLY'), cells=GSE132465_cells) +
  Seurat::NoLegend() + ggplot2::ylab('') +  ggplot2::ggtitle('') +
  ggplot2::theme(axis.ticks.y=ggplot2::element_blank(), axis.text.y=ggplot2::element_blank(), axis.line.y=ggplot2::element_blank())
nk_3 = Seurat::FeaturePlot(CRCatlas_integrated, features=c('rna_GNLY'), cells=GSE144735_cells) +
  Seurat::NoLegend() + ggplot2::xlab('') + ggplot2::ylab('') +  ggplot2::ggtitle('') +
  ggplot2::theme(axis.ticks.y=ggplot2::element_blank(), axis.text.y=ggplot2::element_blank(), axis.line.y=ggplot2::element_blank())
(dim_dt_ann + ggplot2::ggtitle('GNLY')) / (nk_1 | nk_2 | nk_3)

cd31cd38_cells = colnames(CRCatlas_integrated@assays$RNA@counts)[CRCatlas_integrated@assays$RNA@counts['CD38',]>0 &
                                                                   CRCatlas_integrated@assays$RNA@counts['PECAM1',]>0]
cd31cd38_1 = Seurat::DimPlot(CRCatlas_integrated, cells.highlight=cd31cd38_cells, cells=Qian_cells, sizes.highlight=0.3, pt.size=.3) +
  Seurat::NoLegend() + ggplot2::xlab('') +  ggplot2::ggtitle('')
cd31cd38_2 = Seurat::DimPlot(CRCatlas_integrated, cells.highlight=cd31cd38_cells, cells=GSE132465_cells, sizes.highlight=0.3, pt.size=.3) +
  Seurat::NoLegend() + ggplot2::ylab('') +  ggplot2::ggtitle('') +
  ggplot2::theme(axis.ticks.y=ggplot2::element_blank(), axis.text.y=ggplot2::element_blank(), axis.line.y=ggplot2::element_blank())
cd31cd38_3 = Seurat::DimPlot(CRCatlas_integrated, cells.highlight=cd31cd38_cells, cells=GSE144735_cells, sizes.highlight=0.3, pt.size=.3) +
  Seurat::NoLegend() + ggplot2::xlab('') + ggplot2::ylab('') +  ggplot2::ggtitle('') +
  ggplot2::theme(axis.ticks.y=ggplot2::element_blank(), axis.text.y=ggplot2::element_blank(), axis.line.y=ggplot2::element_blank())
(dim_dt_ann + ggplot2::ggtitle('CD38+CD31+ Cells')) / (cd31cd38_1 | cd31cd38_2 | cd31cd38_3)

cd31cd38_cells1 = colnames(CRCatlas_integrated@assays$RNA@counts)[CRCatlas_integrated@assays$RNA@counts['CD38',]>1 &
                                                                    CRCatlas_integrated@assays$RNA@counts['PECAM1',]>1]
cd31cd38_11 = Seurat::DimPlot(CRCatlas_integrated, cells.highlight=cd31cd38_cells1, cells=Qian_cells, sizes.highlight=0.3, pt.size=.3) +
  Seurat::NoLegend() + ggplot2::xlab('') +  ggplot2::ggtitle('')
cd31cd38_21 = Seurat::DimPlot(CRCatlas_integrated, cells.highlight=cd31cd38_cells1, cells=GSE132465_cells, sizes.highlight=0.3, pt.size=.3) +
  Seurat::NoLegend() + ggplot2::ylab('') +  ggplot2::ggtitle('') +
  ggplot2::theme(axis.ticks.y=ggplot2::element_blank(), axis.text.y=ggplot2::element_blank(), axis.line.y=ggplot2::element_blank())
cd31cd38_31 = Seurat::DimPlot(CRCatlas_integrated, cells.highlight=cd31cd38_cells1, cells=GSE144735_cells, sizes.highlight=0.3, pt.size=.3) +
  Seurat::NoLegend() + ggplot2::xlab('') + ggplot2::ylab('') +  ggplot2::ggtitle('') +
  ggplot2::theme(axis.ticks.y=ggplot2::element_blank(), axis.text.y=ggplot2::element_blank(), axis.line.y=ggplot2::element_blank())
(dim_dt_ann + ggplot2::ggtitle('CD38+CD31+ Cells')) / (cd31cd38_11 | cd31cd38_21 | cd31cd38_31)

