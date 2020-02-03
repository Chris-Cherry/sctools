#' Make default plot array associated with processing step.
#' Parameters:
#' @param ser           Seurat object to process.
#' @param out_dir       Output directory
#' @param colors        Custom color scheme
#'   
#' Creates a number of plots and csv files relating to cluster composition
#' by all metadata in the ser object. This includes dim plots, feature plots,
#' and pie plots for cluster composition. Point size is automatically scaled
#' based on total number of cells and 10%/90% quantiles are used for all
#' feature plots. Reductions used are umap and phate.
#'
#' @import grDevices
#' @export

make_processing_plots <- function(ser, out_dir = '1_process/', colors = NULL){

    dir.create(out_dir)
    if(ncol(ser) < 5000){
        pt.size = 2
    } else if(ncol(ser) < 10000){
        pt.size = 1.5
    } else {
        pt.size = 1
    }

    for(meta in colnames(ser[[]])){
        print(meta)
        if(class(ser[[meta]][,1]) == 'factor'){
            if(meta == 'seurat_clusters'){
                cols = colors
            } else {
                cols = NULL
            }
            # Color dim plots by meta
            png(paste0(out_dir, '/', meta, '_umap_nolabel.png'), 
                height = 1000, width = 1000)
            print(Seurat::DimPlot(ser, group.by = meta, cols = cols, pt.size = pt.size))
            dev.off()

            png(paste0(out_dir, '/', meta, '_phate_nolabel.png'), 
                height = 1000, width = 1000)
            print(Seurat::DimPlot(ser, group.by = meta, cols = cols, pt.size = pt.size, reduction = 'phate'))
            dev.off()

            png(paste0(out_dir, '/', meta, '_umap.png'), 
                height = 1000, width = 1000)
            print(Seurat::DimPlot(ser, group.by = meta, cols = cols, pt.size = pt.size, label = TRUE, label.size = 8))
            dev.off()

            png(paste0(out_dir, '/', meta, '_phate.png'), 
                height = 1000, width = 1000)
            print(Seurat::DimPlot(ser, group.by = meta, cols = cols, pt.size = pt.size, label = TRUE, label.size = 8, reduction = 'phate'))
            dev.off()

            # Split dim plots by meta
            n_factor = length(levels(ser[[meta]][,1]))
            ncol = ceiling(sqrt(n_factor))
            nrow = ceiling(n_factor/ncol)
            h = 500*nrow
            w = 500*ncol
            png(paste0(out_dir, '/', meta, '_umap_nolabel_split.png'), 
                height = h, width = w)
            print(Seurat::DimPlot(ser, split.by = meta, cols = colors, pt.size = pt.size, ncol = ncol))
            dev.off()

            png(paste0(out_dir, '/', meta, '_phate_nolabel_split.png'), 
                height = h, width = w)
            print(Seurat::DimPlot(ser, split.by = meta, cols = colors, pt.size = pt.size, ncol = ncol, reduction = 'phate'))
            dev.off()

            png(paste0(out_dir, '/', meta, '_umap_split.png'), 
                height = h, width = w)
            print(Seurat::DimPlot(ser, split.by = meta, cols = colors, pt.size = pt.size, ncol = ncol, label.size = 8, label = TRUE))
            dev.off()

            png(paste0(out_dir, '/', meta, '_phate_split.png'), 
                height = h, width = w)
            print(Seurat::DimPlot(ser, split.by = meta, cols = colors, pt.size = pt.size, ncol = ncol, label.size = 8, label = TRUE, reduction = 'phate'))
            dev.off() 

            clust_proportions(ser, meta, 
                paste0(out_dir, '/', meta, '_proportions.pdf'),
                paste0(out_dir, '/', meta, '_proportions.csv'))           
        }
        if(class(ser[[meta]][,1]) == 'numeric'){
            png(paste0(out_dir, '/', meta, '_umap_feature.png'),
                height = 1000, width = 1000)
            print(Seurat::FeaturePlot(ser, features = meta, pt.size = pt.size))
            dev.off()

            png(paste0(out_dir, '/', meta, '_phate_feature.png'),
                height = 1000, width = 1000)
            print(Seurat::FeaturePlot(ser, features = meta, pt.size = pt.size, reduction = 'phate'))
            dev.off()

            png(paste0(out_dir, '/', meta, '_umap_feature_q10.png'),
                height = 1000, width = 1000)
            print(Seurat::FeaturePlot(ser, features = meta,  pt.size = pt.size, min.cutoff = 'q10', max.cutoff = 'q90'))
            dev.off()

            png(paste0(out_dir, '/', meta, '_phate_feature_q10.png'),
                height = 1000, width = 1000)
            print(Seurat::FeaturePlot(ser, features = meta,  pt.size = pt.size, reduction = 'phate', min.cutoff = 'q10', max.cutoff = 'q90'))
            dev.off()
        }
    }
}
