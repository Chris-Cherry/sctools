#' Make violin and feature plots for a given set of features
#'
#' This function accepts a Seurat object and a vector of features to create
#' violin and feature plots. It chunks the vector into peices containing up to
#' 9 features each. The output folder will contain a pdf of violin plots, each
#' page containing one chunk worth of features, and two png files for each
#' chunk - one feature plot of the chunk with no thresholding and one with
#' thresholding on the 10% and 90% quantile.
#' 
#' @param ser           Seurat file to use for plots
#' @param features      Vector of features to make plots for
#' @param out_dir       Directory to write plots
#' @param use_phate     Boolean indicating whether to use phate dimensional reductions
#' @import Seurat
#' @importFrom grDevices png
#' @importFrom grDevices pdf
#' @importFrom grDevices dev.off
#' @export

make_tar_feat_plots <- function(ser, features, out_dir = '2_de', 
	prefix = 'tar_features', use_phate = TRUE){
    dir.create(out_dir)
    prefix = paste(out_dir, prefix, sep = '/')
    tmp = ser
    DefaultAssay(ser) = 'RNA'

    fp_cols = c("#2c7bb6", "#00a6ca", "#00ccbc", "#90eb9d", "#ffff8c", "#f9d057",
        "#f29e2e", "#e76818", "#d7191c")

    sets = split(features, ceiling(seq_along(features)/12))
    pdf(paste0(prefix, '_vlns.pdf'))
    for(set in sets){
        suppressWarnings(print(Seurat::VlnPlot(ser, features = set, pt.size = 0, 
            ncol = 3)))
    }
    dev.off()

    for(i in 1:length(sets)){
        set = sets[[i]]

        png(paste0(prefix, '_feature_plot_', i, '.png'), height = 2000, 
            width = 2000)
        suppressWarnings(print(Seurat::FeaturePlot(ser, features = set, 
            ncol = 3, cols = fp_cols, reduction = 'umap')))
        dev.off()

        png(paste0(prefix, '_feature_plot_', i, '_q10.png'), height = 2000, 
            width = 2000)
        suppressWarnings(print(Seurat::FeaturePlot(ser, features = set, ncol = 3, 
            min.cutoff = 'q10', max.cutoff = 'q90', cols = fp_cols, 
            reduction = 'umap')))
        dev.off()

        if(use_phate){
            png(paste0(prefix, '_feature_plot_', i, '_phate.png'), height = 2000, 
                width = 2000)
            suppressWarnings(print(Seurat::FeaturePlot(ser, features = set, 
                ncol = 3, cols = fp_cols, reduction = 'phate')))
            dev.off()

            png(paste0(prefix, '_feature_plot_', i, '_q10_phate.png'), height = 2000, 
                width = 2000)
            suppressWarnings(print(Seurat::FeaturePlot(ser, features = set, ncol = 3, 
                min.cutoff = 'q10', max.cutoff = 'q90', cols = fp_cols, 
                reduction = 'phate')))
            dev.off()
        }
    }
}
