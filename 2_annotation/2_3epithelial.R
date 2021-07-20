project_dir = '/home/scardoso/Documents/PhD/CRC_ATLAS'
source(paste(project_dir, 'utils/modified_plots.R', sep='/'))

# ############################################### #
# #### CHECK GENE MARKERS OF TUMOUR VS NORMAL ### #
# ############################################### #

# Load dataset of epithelial cells:
Epithelial = SeuratDisk::LoadH5Seurat(paste(project_dir, '2_annotation/results_Epithelial/datasets/Epithelial3.h5Seurat', sep='/'))
invisible(gc())





# -----
# - Check CNV results
# -----


# 1. Plots colored by state and patient:
state = Seurat::DimPlot(Epithelial, reduction="umap", group.by='state', label=FALSE, pt.size=.3) + ggplot2::theme_minimal()
patients = Seurat::DimPlot(Epithelial, reduction="umap", group.by='patient', label=TRUE, label.size=3, pt.size=.3) +
  ggplot2::theme_minimal() + Seurat::NoLegend()


# 2. Color by CNA annotation:
cna_annotation = Seurat::DimPlot(Epithelial, reduction="umap", group.by='CNA_annotation', label=FALSE, pt.size=.3) +
  ggplot2::theme_minimal() + ggplot2::theme(legend.position = 'bottom')
(state + ggplot2::theme(legend.position = 'bottom')) | cna_annotation





# -----
# - Check genes
# -----

genes_to_check = c('rna_MYC', 'rna_RNF43', 'rna_AXIN2', 'rna_CTNNB1', 'rna_CD44', 'rna_MLH1')

patients_putNorm = c('SMC03', 'SMC10', 'SMC19', 'SMC24', 'SMC07', 'KUL01', '31')
Epithelial_noPutNorm = subset(Epithelial, cells=rownames(Epithelial[[]])[!Epithelial$patient%in%patients_putNorm])
invisible(gc())

# 1. Feature and violin plots only with patients with no putative normal cells
feature_plots(Epithelial_noPutNorm, genes_to_check, ncol=3, with_dimplot=FALSE, point_size=.4)
violin_plots(Epithelial_noPutNorm, genes_to_check, ncol=3, with_dimplot=FALSE, group.by='CNA_annotation')
# Difference in expression of these genes between normal and tumour is as expected (in literature)


# 3. Patient SMC03
cells_use_SMC03 = rownames(Epithelial[[]])[Epithelial$patient=='SMC03' | Epithelial$state=='Normal']
SMC03_epithelial = subset(Epithelial, cells=cells_use_SMC03)
SMC03_df = Seurat::FetchData(SMC03_epithelial, c(genes_to_check, 'CNA_annotation'))
# 3.1. Boxplots
SMC03_boxplots = list()
for(gene in  genes_to_check){
  SMC03_boxplots[[gene]] = ggplot2::ggplot(SMC03_df, ggplot2::aes_string('CNA_annotation', gene, colour='CNA_annotation')) +
    ggplot2::geom_boxplot() + ggplot2::theme_minimal() + ggplot2::ylab(gene) + ggplot2::theme(legend.position = 'none')
}
patchwork::wrap_plots(SMC03_boxplots, ncol=3)
# 3.2. We will focus on genes CD44 and CTNNB1
CD44_cutoff_SMC03 = max(summary(SMC03_df$rna_CD44[SMC03_df$CNA_annotation=='Normal'])[2:5])
CTNNB1_cutoff_SMC03 = max(summary(SMC03_df$rna_CTNNB1[SMC03_df$CNA_annotation=='Normal'])[2:5])
SMC03_df_new = SMC03_df
SMC03_df_new$CNA_annotation[(SMC03_df_new$rna_CD44 > CD44_cutoff_SMC03 | SMC03_df_new$rna_CTNNB1 > CTNNB1_cutoff_SMC03)&
                              SMC03_df_new$CNA_annotation=='Putative Normal'] = 'New Putative Tumour'
SMC03_df_new2 = SMC03_df_new
SMC03_df_new2$CNA_annotation[SMC03_df_new$CNA_annotation=='New Putative Tumour'] = 'Putative Tumour'
table(SMC03_df_new$CNA_annotation)
# 3.3. Check boxplots with new annotations
CD44_new = ggplot2::ggplot(SMC03_df_new, ggplot2::aes_string('CNA_annotation', 'rna_CD44', colour='CNA_annotation')) +
  ggplot2::geom_boxplot() + ggplot2::theme_minimal() + ggplot2::ylab('rna_CD44') + ggplot2::theme(legend.position = 'none')
CTNNB1_new = ggplot2::ggplot(SMC03_df_new, ggplot2::aes_string('CNA_annotation', 'rna_CTNNB1', colour='CNA_annotation')) +
  ggplot2::geom_boxplot() + ggplot2::theme_minimal() + ggplot2::ylab('rna_CTNNB1') + ggplot2::theme(legend.position = 'none')
CD44_new | CTNNB1_new
# 3.4. Check boxplots normal vs putative normal vs putative tumour before and after new assignments
boxplots_all_genes_before = list()
for(gene in  genes_to_check){
  boxplots_all_genes_before[[gene]] = ggplot2::ggplot(SMC03_df, ggplot2::aes_string('CNA_annotation', gene, colour='CNA_annotation')) +
    ggplot2::geom_boxplot() + ggplot2::theme_minimal() + ggplot2::ylab(gene) + ggplot2::theme(legend.position = 'none')
}
patchwork::wrap_plots(boxplots_all_genes_before, ncol=3)
boxplots_all_genes_after = list()
for(gene in  genes_to_check){
  boxplots_all_genes_after[[gene]] = ggplot2::ggplot(SMC03_df_new2, ggplot2::aes_string('CNA_annotation', gene, colour='CNA_annotation')) +
    ggplot2::geom_boxplot() + ggplot2::theme_minimal() + ggplot2::ylab(gene) + ggplot2::theme(legend.position = 'none')
}
patchwork::wrap_plots(boxplots_all_genes_after, ncol=3)
# 3.5. Check UMAP with new assignments
normcells_before = rownames(SMC03_df)[SMC03_df$CNA_annotation=='Putative Normal']
normacells_after = rownames(SMC03_df_new2)[SMC03_df_new2$CNA_annotation=='Putative Normal']
badNormSMC03_highlighted = Seurat::DimPlot(Epithelial, cells.highlight = normcells_before, sizes.highlight = .2, pt.size=.1) +
  ggplot2::ggtitle('Putative normal - Before') + ggplot2::theme_minimal() + Seurat::NoLegend()
goodNormSMC03_highlighted = Seurat::DimPlot(Epithelial, cells.highlight = normacells_after, sizes.highlight = .2, pt.size=.1) +
  ggplot2::ggtitle('Putative normal - After') + ggplot2::theme_minimal() + Seurat::NoLegend()
(badNormSMC03_highlighted / goodNormSMC03_highlighted) | state


# 4. Patient SMC10
cells_use_SMC10 = rownames(Epithelial[[]])[Epithelial$patient=='SMC10' | Epithelial$state=='Normal']
SMC10_epithelial = subset(Epithelial, cells=cells_use_SMC10)
SMC10_df = Seurat::FetchData(SMC10_epithelial, c(genes_to_check, 'CNA_annotation'))
# 4.1. Boxplots
SMC10_boxplots = list()
for(gene in  genes_to_check){
  SMC10_boxplots[[gene]] = ggplot2::ggplot(SMC10_df, ggplot2::aes_string('CNA_annotation', gene, colour='CNA_annotation')) +
    ggplot2::geom_boxplot() + ggplot2::theme_minimal() + ggplot2::ylab(gene) + ggplot2::theme(legend.position = 'none')
}
patchwork::wrap_plots(SMC10_boxplots, ncol=3)
# 4.2. We will focus on genes RNF43, CD44 and CTNNB1
RNF43_cutoff_SMC10 = max(summary(SMC10_df$rna_RNF43[SMC10_df$CNA_annotation=='Normal'])[2:5])
CD44_cutoff_SMC10 = max(summary(SMC10_df$rna_CD44[SMC10_df$CNA_annotation=='Normal'])[2:5])
CTNNB1_cutoff_SMC10 = max(summary(SMC10_df$rna_CTNNB1[SMC10_df$CNA_annotation=='Normal'])[2:5])
SMC10_df_new = SMC10_df
SMC10_df_new$CNA_annotation[(SMC10_df_new$rna_RNF43 > RNF43_cutoff_SMC10 | SMC10_df_new$rna_CD44 > CD44_cutoff_SMC10 | 
                               SMC10_df_new$rna_CTNNB1 > CTNNB1_cutoff_SMC10) &
                              SMC10_df_new$CNA_annotation=='Putative Normal'] = 'New Putative Tumour'
SMC10_df_new2 = SMC10_df_new
SMC10_df_new2$CNA_annotation[SMC10_df_new$CNA_annotation=='New Putative Tumour'] = 'Putative Tumour'
table(SMC10_df_new$CNA_annotation)
# 4.3. Check boxplots with new annotations
RNF43_new = ggplot2::ggplot(SMC10_df_new, ggplot2::aes_string('CNA_annotation', 'rna_RNF43', colour='CNA_annotation')) +
  ggplot2::geom_boxplot() + ggplot2::theme_minimal() + ggplot2::ylab('rna_RNF43') + ggplot2::theme(legend.position = 'none')
CD44_new = ggplot2::ggplot(SMC10_df_new, ggplot2::aes_string('CNA_annotation', 'rna_CD44', colour='CNA_annotation')) +
  ggplot2::geom_boxplot() + ggplot2::theme_minimal() + ggplot2::ylab('rna_CD44') + ggplot2::theme(legend.position = 'none')
CTNNB1_new = ggplot2::ggplot(SMC10_df_new, ggplot2::aes_string('CNA_annotation', 'rna_CTNNB1', colour='CNA_annotation')) +
  ggplot2::geom_boxplot() + ggplot2::theme_minimal() + ggplot2::ylab('rna_CTNNB1') + ggplot2::theme(legend.position = 'none')
RNF43_new | CD44_new | CTNNB1_new
# 4.4. Check boxplots normal vs putative normal vs putative tumour before and after new assignments
boxplots_all_genes_before = list()
for(gene in  genes_to_check){
  boxplots_all_genes_before[[gene]] = ggplot2::ggplot(SMC10_df, ggplot2::aes_string('CNA_annotation', gene, colour='CNA_annotation')) +
    ggplot2::geom_boxplot() + ggplot2::theme_minimal() + ggplot2::ylab(gene) + ggplot2::theme(legend.position = 'none')
}
patchwork::wrap_plots(boxplots_all_genes_before, ncol=3)
boxplots_all_genes_after = list()
for(gene in  genes_to_check){
  boxplots_all_genes_after[[gene]] = ggplot2::ggplot(SMC10_df_new2, ggplot2::aes_string('CNA_annotation', gene, colour='CNA_annotation')) +
    ggplot2::geom_boxplot() + ggplot2::theme_minimal() + ggplot2::ylab(gene) + ggplot2::theme(legend.position = 'none')
}
patchwork::wrap_plots(boxplots_all_genes_after, ncol=3)
# 4.5. Check UMAP with new assignments
normcells_before = rownames(SMC10_df)[SMC10_df$CNA_annotation=='Putative Normal']
normacells_after = rownames(SMC10_df_new2)[SMC10_df_new2$CNA_annotation=='Putative Normal']
badNormSMC10_highlighted = Seurat::DimPlot(Epithelial, cells.highlight = normcells_before, sizes.highlight = .2, pt.size=.1) +
  ggplot2::ggtitle('Putative normal - Before') + ggplot2::theme_minimal() + Seurat::NoLegend()
goodNormSMC10_highlighted = Seurat::DimPlot(Epithelial, cells.highlight = normacells_after, sizes.highlight = .2, pt.size=.1) +
  ggplot2::ggtitle('Putative normal - After') + ggplot2::theme_minimal() + Seurat::NoLegend()
(badNormSMC10_highlighted / goodNormSMC10_highlighted) | state


# 5. Patient SMC19
cells_use_SMC19 = rownames(Epithelial[[]])[Epithelial$patient=='SMC19' | Epithelial$state=='Normal']
SMC19_epithelial = subset(Epithelial, cells=cells_use_SMC19)
SMC19_df = Seurat::FetchData(SMC19_epithelial, c(genes_to_check, 'CNA_annotation'))
# 5.1. Boxplots
SMC19_boxplots = list()
for(gene in  genes_to_check){
  SMC19_boxplots[[gene]] = ggplot2::ggplot(SMC19_df, ggplot2::aes_string('CNA_annotation', gene, colour='CNA_annotation')) +
    ggplot2::geom_boxplot() + ggplot2::theme_minimal() + ggplot2::ylab(gene) + ggplot2::theme(legend.position = 'none')
}
patchwork::wrap_plots(SMC19_boxplots, ncol=3)
# 5.2. We will focus on genes MYC, RNF43 and CD44
MYC_cutoff_SMC19 = max(summary(SMC19_df$rna_MYC[SMC19_df$CNA_annotation=='Normal'])[2:5])
RNF43_cutoff_SMC19 = max(summary(SMC19_df$rna_RNF43[SMC19_df$CNA_annotation=='Normal'])[2:5])
CD44_cutoff_SMC19 = max(summary(SMC19_df$rna_CD44[SMC19_df$CNA_annotation=='Normal'])[2:5])
SMC19_df_new = SMC19_df
SMC19_df_new$CNA_annotation[(SMC19_df_new$rna_MYC > MYC_cutoff_SMC19 | SMC19_df_new$rna_RNF43 > RNF43_cutoff_SMC19 |
                              SMC19_df_new$rna_CD44 > CD44_cutoff_SMC19) &
                              SMC19_df_new$CNA_annotation=='Putative Normal'] = 'New Putative Tumour'
SMC19_df_new2 = SMC19_df_new
SMC19_df_new2$CNA_annotation[SMC19_df_new$CNA_annotation=='New Putative Tumour'] = 'Putative Tumour'
table(SMC19_df_new$CNA_annotation)
# 5.3. Check boxplots with new annotations
MYC_new = ggplot2::ggplot(SMC10_df_new, ggplot2::aes_string('CNA_annotation', 'rna_MYC', colour='CNA_annotation')) +
  ggplot2::geom_boxplot() + ggplot2::theme_minimal() + ggplot2::ylab('rna_MYC') + ggplot2::theme(legend.position = 'none')
RNF43_new = ggplot2::ggplot(SMC10_df_new, ggplot2::aes_string('CNA_annotation', 'rna_RNF43', colour='CNA_annotation')) +
  ggplot2::geom_boxplot() + ggplot2::theme_minimal() + ggplot2::ylab('rna_RNF43') + ggplot2::theme(legend.position = 'none')
CD44_new = ggplot2::ggplot(SMC10_df_new, ggplot2::aes_string('CNA_annotation', 'rna_CD44', colour='CNA_annotation')) +
  ggplot2::geom_boxplot() + ggplot2::theme_minimal() + ggplot2::ylab('rna_CD44') + ggplot2::theme(legend.position = 'none')
MYC_new | RNF43_new | CD44_new
# 5.4. Check boxplots normal vs putative normal vs putative tumour before and after new assignments
boxplots_all_genes_before = list()
for(gene in  genes_to_check){
  boxplots_all_genes_before[[gene]] = ggplot2::ggplot(SMC19_df, ggplot2::aes_string('CNA_annotation', gene, colour='CNA_annotation')) +
    ggplot2::geom_boxplot() + ggplot2::theme_minimal() + ggplot2::ylab(gene) + ggplot2::theme(legend.position = 'none')
}
patchwork::wrap_plots(boxplots_all_genes_before, ncol=3)
boxplots_all_genes_after = list()
for(gene in  genes_to_check){
  boxplots_all_genes_after[[gene]] = ggplot2::ggplot(SMC19_df_new2, ggplot2::aes_string('CNA_annotation', gene, colour='CNA_annotation')) +
    ggplot2::geom_boxplot() + ggplot2::theme_minimal() + ggplot2::ylab(gene) + ggplot2::theme(legend.position = 'none')
}
patchwork::wrap_plots(boxplots_all_genes_after, ncol=3)
# 5.5. Check UMAP with new assignments
normcells_before = rownames(SMC19_df)[SMC19_df$CNA_annotation=='Putative Normal']
normacells_after = rownames(SMC19_df_new2)[SMC19_df_new2$CNA_annotation=='Putative Normal']
badNormSMC19_highlighted = Seurat::DimPlot(Epithelial, cells.highlight = normcells_before, sizes.highlight = .2, pt.size=.1) +
  ggplot2::ggtitle('Putative normal - Before') + ggplot2::theme_minimal() + Seurat::NoLegend()
goodNormSMC19_highlighted = Seurat::DimPlot(Epithelial, cells.highlight = normacells_after, sizes.highlight = .2, pt.size=.1) +
  ggplot2::ggtitle('Putative normal - After') + ggplot2::theme_minimal() + Seurat::NoLegend()
(badNormSMC19_highlighted / goodNormSMC19_highlighted) | state


# 6. Patient SMC24
cells_use_SMC24 = rownames(Epithelial[[]])[Epithelial$patient=='SMC24' | Epithelial$state=='Normal']
SMC24_epithelial = subset(Epithelial, cells=cells_use_SMC24)
SMC24_df = Seurat::FetchData(SMC24_epithelial, c(genes_to_check, 'CNA_annotation'))
# 6.1. Boxplots
SMC24_boxplots = list()
for(gene in  genes_to_check){
  SMC24_boxplots[[gene]] = ggplot2::ggplot(SMC24_df, ggplot2::aes_string('CNA_annotation', gene, colour='CNA_annotation')) +
    ggplot2::geom_boxplot() + ggplot2::theme_minimal() + ggplot2::ylab(gene) + ggplot2::theme(legend.position = 'none')
}
patchwork::wrap_plots(SMC24_boxplots, ncol=3)
# 6.2. We will focus on genes MYC, RNF43, CTNNB1 and CD44
MYC_cutoff_SMC24 = max(summary(SMC24_df$rna_MYC[SMC24_df$CNA_annotation=='Normal'])[2:5])
RNF43_cutoff_SMC24 = max(summary(SMC24_df$rna_RNF43[SMC24_df$CNA_annotation=='Normal'])[2:5])
CTNNB1_cutoff_SMC24 = max(summary(SMC24_df$rna_CTNNB1[SMC24_df$CNA_annotation=='Normal'])[2:5])
CD44_cutoff_SMC24 = max(summary(SMC24_df$rna_CD44[SMC24_df$CNA_annotation=='Normal'])[2:5])
SMC24_df_new = SMC24_df
SMC24_df_new$CNA_annotation[(SMC24_df_new$rna_MYC > MYC_cutoff_SMC24 | SMC24_df_new$rna_RNF43 > RNF43_cutoff_SMC24 |
                               SMC24_df_new$rna_CTNNB1 > CTNNB1_cutoff_SMC24 | SMC24_df_new$rna_CD44 > CD44_cutoff_SMC24) &
                              SMC24_df_new$CNA_annotation=='Putative Normal'] = 'New Putative Tumour'
SMC24_df_new2 = SMC24_df_new
SMC24_df_new2$CNA_annotation[SMC24_df_new$CNA_annotation=='New Putative Tumour'] = 'Putative Tumour'
table(SMC24_df_new$CNA_annotation)
# 6.3. Check boxplots with new annotations
MYC_new = ggplot2::ggplot(SMC24_df_new, ggplot2::aes_string('CNA_annotation', 'rna_MYC', colour='CNA_annotation')) +
  ggplot2::geom_boxplot() + ggplot2::theme_minimal() + ggplot2::ylab('rna_MYC') + ggplot2::theme(legend.position = 'none')
RNF43_new = ggplot2::ggplot(SMC24_df_new, ggplot2::aes_string('CNA_annotation', 'rna_RNF43', colour='CNA_annotation')) +
  ggplot2::geom_boxplot() + ggplot2::theme_minimal() + ggplot2::ylab('rna_RNF43') + ggplot2::theme(legend.position = 'none')
CTNNB1_new = ggplot2::ggplot(SMC24_df_new, ggplot2::aes_string('CNA_annotation', 'rna_CTNNB1', colour='CNA_annotation')) +
  ggplot2::geom_boxplot() + ggplot2::theme_minimal() + ggplot2::ylab('rna_CTNNB1') + ggplot2::theme(legend.position = 'none')
CD44_new = ggplot2::ggplot(SMC24_df_new, ggplot2::aes_string('CNA_annotation', 'rna_CD44', colour='CNA_annotation')) +
  ggplot2::geom_boxplot() + ggplot2::theme_minimal() + ggplot2::ylab('rna_CD44') + ggplot2::theme(legend.position = 'none')
(MYC_new | RNF43_new) / (CTNNB1_new | CD44_new)
# 6.4. Check boxplots normal vs putative normal vs putative tumour before and after new assignments
boxplots_all_genes_before = list()
for(gene in  genes_to_check){
  boxplots_all_genes_before[[gene]] = ggplot2::ggplot(SMC24_df, ggplot2::aes_string('CNA_annotation', gene, colour='CNA_annotation')) +
    ggplot2::geom_boxplot() + ggplot2::theme_minimal() + ggplot2::ylab(gene) + ggplot2::theme(legend.position = 'none')
}
patchwork::wrap_plots(boxplots_all_genes_before, ncol=3)
boxplots_all_genes_after = list()
for(gene in  genes_to_check){
  boxplots_all_genes_after[[gene]] = ggplot2::ggplot(SMC24_df_new2, ggplot2::aes_string('CNA_annotation', gene, colour='CNA_annotation')) +
    ggplot2::geom_boxplot() + ggplot2::theme_minimal() + ggplot2::ylab(gene) + ggplot2::theme(legend.position = 'none')
}
patchwork::wrap_plots(boxplots_all_genes_after, ncol=3)
# 6.5. Check UMAP with new assignments
normcells_before = rownames(SMC24_df)[SMC24_df$CNA_annotation=='Putative Normal']
normacells_after = rownames(SMC24_df_new2)[SMC24_df_new2$CNA_annotation=='Putative Normal']
badNormSMC24_highlighted = Seurat::DimPlot(Epithelial, cells.highlight = normcells_before, sizes.highlight = .2, pt.size=.1) +
  ggplot2::ggtitle('Putative normal - Before') + ggplot2::theme_minimal() + Seurat::NoLegend()
goodNormSMC24_highlighted = Seurat::DimPlot(Epithelial, cells.highlight = normacells_after, sizes.highlight = .2, pt.size=.1) +
  ggplot2::ggtitle('Putative normal - After') + ggplot2::theme_minimal() + Seurat::NoLegend()
(badNormSMC24_highlighted / goodNormSMC24_highlighted) | state


# 7. Patient KUL01
cells_use_KUL01 = rownames(Epithelial[[]])[Epithelial$patient=='KUL01' | Epithelial$state=='Normal']
KUL01_epithelial = subset(Epithelial, cells=cells_use_KUL01)
KUL01_df = Seurat::FetchData(KUL01_epithelial, c(genes_to_check, 'CNA_annotation'))
# 7.1. Boxplots
KUL01_boxplots = list()
for(gene in  genes_to_check){
  KUL01_boxplots[[gene]] = ggplot2::ggplot(KUL01_df, ggplot2::aes_string('CNA_annotation', gene, colour='CNA_annotation')) +
    ggplot2::geom_boxplot() + ggplot2::theme_minimal() + ggplot2::ylab(gene) + ggplot2::theme(legend.position = 'none')
}
patchwork::wrap_plots(KUL01_boxplots, ncol=3)
# 7.2. No changes will be made to the CNA assignments


# 8. Patient 31
cells_use_31 = rownames(Epithelial[[]])[Epithelial$patient=='31' | Epithelial$state=='Normal']
p31_epithelial = subset(Epithelial, cells=cells_use_31)
p31_df = Seurat::FetchData(p31_epithelial, c(genes_to_check, 'CNA_annotation'))
# 8.1. Boxplots
p31_boxplots = list()
for(gene in  genes_to_check){
  p31_boxplots[[gene]] = ggplot2::ggplot(p31_df, ggplot2::aes_string('CNA_annotation', gene, colour='CNA_annotation')) +
    ggplot2::geom_boxplot() + ggplot2::theme_minimal() + ggplot2::ylab(gene) + ggplot2::theme(legend.position = 'none')
}
patchwork::wrap_plots(p31_boxplots, ncol=3)
# 8.2. No changes will be made to the CNA assignments


# 9. Patient SMC07
cells_use_SMC07 = rownames(Epithelial[[]])[Epithelial$patient=='SMC07' | Epithelial$state=='Normal']
SMC07_epithelial = subset(Epithelial, cells=cells_use_SMC07)
SMC07_df = Seurat::FetchData(SMC07_epithelial, c(genes_to_check, 'CNA_annotation'))
# 9.1. Boxplots
SMC07_boxplots = list()
for(gene in  genes_to_check){
  SMC07_boxplots[[gene]] = ggplot2::ggplot(SMC07_df, ggplot2::aes_string('CNA_annotation', gene, colour='CNA_annotation')) +
    ggplot2::geom_boxplot() + ggplot2::theme_minimal() + ggplot2::ylab(gene) + ggplot2::theme(legend.position = 'none')
}
patchwork::wrap_plots(SMC07_boxplots, ncol=3)
# 9.2. We will focus on gene CD44
CD44_cutoff_SMC07 = max(summary(SMC07_df$rna_CD44[SMC07_df$CNA_annotation=='Normal'])[2:5])
SMC07_df_new = SMC07_df
SMC07_df_new$CNA_annotation[(SMC07_df_new$rna_CD44 > CD44_cutoff_SMC07) &
                              SMC07_df_new$CNA_annotation=='Putative Normal'] = 'New Putative Tumour'
SMC07_df_new2 = SMC07_df_new
SMC07_df_new2$CNA_annotation[SMC07_df_new$CNA_annotation=='New Putative Tumour'] = 'Putative Tumour'
table(SMC07_df_new$CNA_annotation)
# 9.3. Check boxplots with new annotations
CD44_new = ggplot2::ggplot(SMC07_df_new, ggplot2::aes_string('CNA_annotation', 'rna_CD44', colour='CNA_annotation')) +
  ggplot2::geom_boxplot() + ggplot2::theme_minimal() + ggplot2::ylab('rna_CD44') + ggplot2::theme(legend.position = 'none')
CD44_new
# 9.4. Check boxplots normal vs putative normal vs putative tumour before and after new assignments
boxplots_all_genes_before = list()
for(gene in  genes_to_check){
  boxplots_all_genes_before[[gene]] = ggplot2::ggplot(SMC07_df, ggplot2::aes_string('CNA_annotation', gene, colour='CNA_annotation')) +
    ggplot2::geom_boxplot() + ggplot2::theme_minimal() + ggplot2::ylab(gene) + ggplot2::theme(legend.position = 'none')
}
patchwork::wrap_plots(boxplots_all_genes_before, ncol=3)
boxplots_all_genes_after = list()
for(gene in  genes_to_check){
  boxplots_all_genes_after[[gene]] = ggplot2::ggplot(SMC07_df_new2, ggplot2::aes_string('CNA_annotation', gene, colour='CNA_annotation')) +
    ggplot2::geom_boxplot() + ggplot2::theme_minimal() + ggplot2::ylab(gene) + ggplot2::theme(legend.position = 'none')
}
patchwork::wrap_plots(boxplots_all_genes_after, ncol=3)
# 9.5. Check UMAP with new assignments
normcells_before = rownames(SMC07_df)[SMC07_df$CNA_annotation=='Putative Normal']
normacells_after = rownames(SMC07_df_new2)[SMC07_df_new2$CNA_annotation=='Putative Normal']
badNormSMC07_highlighted = Seurat::DimPlot(Epithelial, cells.highlight = normcells_before, sizes.highlight = .2, pt.size=.1) +
  ggplot2::ggtitle('Putative normal - Before') + ggplot2::theme_minimal() + Seurat::NoLegend()
goodNormSMC07_highlighted = Seurat::DimPlot(Epithelial, cells.highlight = normacells_after, sizes.highlight = .2, pt.size=.1) +
  ggplot2::ggtitle('Putative normal - After') + ggplot2::theme_minimal() + Seurat::NoLegend()
(badNormSMC07_highlighted / goodNormSMC07_highlighted) | state


# 10. Set new assignments:
new_tumour = c(rownames(SMC03_df_new)[SMC03_df_new$CNA_annotation=='New Putative Tumour'],
               rownames(SMC10_df_new)[SMC10_df_new$CNA_annotation=='New Putative Tumour'],
               rownames(SMC19_df_new)[SMC19_df_new$CNA_annotation=='New Putative Tumour'],
               rownames(SMC24_df_new)[SMC24_df_new$CNA_annotation=='New Putative Tumour'],
               rownames(SMC07_df_new)[SMC07_df_new$CNA_annotation=='New Putative Tumour'])
Epithelial[['NormTum_assignment']] = Epithelial$CNA_annotation
Epithelial@meta.data[new_tumour, 'NormTum_assignment'] = 'Putative Tumour'
table(Epithelial$NormTum_assignment) # We now have 478 Putative Normal cells


# 11. How do the normal epithelial cells (normal + putative) group together?
normal_epithelial = subset(Epithelial, cells=rownames(Epithelial[[]])[Epithelial$NormTum_assignment!='Putative Tumour'])
normal_epithelial = Seurat::RunPCA(normal_epithelial, assay='integrated')
pct = normal_epithelial[["pca"]]@stdev /sum(normal_epithelial[["pca"]]@stdev) * 100
cumu = cumsum(pct)
co1 = which(cumu > 90 & pct < 5)[1]
co2 = sort(which((pct[1:length(pct) - 1] - pct[2:length(pct)]) > 0.1), decreasing = T)[1] + 1
elbow = min(co1, co2) # 14
normal_epithelial = Seurat::RunUMAP(normal_epithelial, dims=1:elbow)
Seurat::DimPlot(normal_epithelial, group.by = 'NormTum_assignment', pt.size=.5) + ggplot2::theme_minimal() |
  Seurat::DimPlot(normal_epithelial, group.by = 'state', pt.size=.5) + ggplot2::theme_minimal()


# 12. Separate normal from tumour epithelial cells and store datasets.
# 12.1. Normal
normal_epithelial = subset(Epithelial, cells=rownames(Epithelial[[]])[Epithelial$NormTum_assignment!='Putative Tumour'])
SeuratDisk::SaveH5Seurat(normal_epithelial, paste(project_dir, '2_annotation/results_Epithelial/datasets/Epithelial_normal.h5Seurat', sep='/'))
invisible(gc())

# 12.2. Tumour
tumour_epithelial = subset(Epithelial, cells=rownames(Epithelial[[]])[Epithelial$NormTum_assignment=='Putative Tumour'])
SeuratDisk::SaveH5Seurat(tumour_epithelial, paste(project_dir, '2_annotation/results_Epithelial/datasets/Epithelial_tumour.h5Seurat', sep='/'))
invisible(gc())
