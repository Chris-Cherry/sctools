#' Make violin and feature plots for a given set of genes
#' Parameters:
#' @param ser           Seurat file to use for plots
#' @param genes         Vector of genes to make plots for
#' @param out_dir       Directory to write plots
#' @export

make_tar_gene_plots <- function(ser, genes, out_dir = '2_de', 
	prefix = '2_de/tar_genes'){
    
    tmp = ser
    DefaultAssay(ser) = 'RNA'

    sets = split(genes, ceiling(seq_along(genes)/12))
    pdf(paste0(prefix, '_vlns.pdf'))
    for(set in sets){
        print(Seurat::VlnPlot(ser, features = set, pt.size = 0, assay = 'RNA', ncol = 3))
    }
    dev.off()

    for(i in 1:length(sets)){
        set = sets[[i]]

        png(paste0(prefix, '_feature_plot_', i, '.png'), height = 1000, width = 1000)
        print(Seurat::FeaturePlot(ser, features = set, ncol = 3))
        dev.off()

        png(paste0(prefix, '_feature_plot_', i, '_q10.png'), height = 1000, width = 1000)
        print(Seurat::FeaturePlot(ser, features = set, ncol = 3, 
            min.cutoff = 'q10', max.cutoff = 'q90'))
        dev.off()
    }
}
