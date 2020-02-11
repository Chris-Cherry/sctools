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
#' @export

make_tar_feat_plots <- function(ser, features, out_dir = '2_de', 
	prefix = 'tar_features'){
    dir.create(out_dir)
    prefix = paste(out_dir, prefix, sep = '/')
    tmp = ser
    DefaultAssay(ser) = 'RNA'

    sets = split(features, ceiling(seq_along(features)/12))
    pdf(paste0(prefix, '_vlns.pdf'))
    for(set in sets){
        suppressWarnings(print(Seurat::VlnPlot(ser, features = set, pt.size = 0, 
            ncol = 3)))
    }
    dev.off()

    for(i in 1:length(sets)){
        set = sets[[i]]

        png(paste0(prefix, '_feature_plot_', i, '.png'), height = 1000, 
            width = 1000)
        suppressWarnings(print(Seurat::FeaturePlot(ser, features = set, 
            ncol = 3)))
        dev.off()

        png(paste0(prefix, '_feature_plot_', i, '_q10.png'), height = 1000, 
            width = 1000)
        suppressWarnings(print(Seurat::FeaturePlot(ser, features = set, ncol = 3, 
            min.cutoff = 'q10', max.cutoff = 'q90')))
        dev.off()
    }
}
