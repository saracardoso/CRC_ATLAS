feature_plots = function(seurat_object, features, ncol=2, with_dimplot=T, dimplot_group='seurat_clusters', point_size=0.2){
  feature_plots = list()
  for(feature in features){
    feature_plots[[feature]] = Seurat::FeaturePlot(seurat_object, features=feature, pt.size=point_size) +
      ggplot2::scale_color_gradientn(colours=c("navy", "yellow"), na.value=grDevices::rgb(0.75, 0.75, 0.75, alpha=0.2), limits=c(1e-5,NA)) +
      ggplot2::theme_minimal()
  }
  if(with_dimplot) feature_plots[['dimplot']] = Seurat::DimPlot(seurat_object, reduction="umap", group.by=dimplot_group,
                                                                label=TRUE, label.size=3, pt.size=point_size) + ggplot2::theme_minimal() + 
                                                        Seurat::NoLegend()
  patchwork::wrap_plots(feature_plots, ncol=ncol)
}


violin_plots = function(seurat_object, features, ncol=2, with_dimplot=T, group.by='seurat_clusters'){
  violin_plots = list()
  for(feature in features){
    violin_plots[[feature]] = Seurat::VlnPlot(seurat_object, features=feature, group.by=group.by, pt.size=0) +
      ggplot2::theme_minimal() + Seurat::NoLegend() + Seurat::RotatedAxis()
  }
  if(with_dimplot) violin_plots[['dimplot']] = Seurat::DimPlot(seurat_object, reduction="umap", group.by=group.by,
                                                                label=TRUE, label.size=3, pt.size=.2) + ggplot2::theme_minimal() +
                                                      Seurat::NoLegend()
  patchwork::wrap_plots(violin_plots, ncol=ncol)
}
