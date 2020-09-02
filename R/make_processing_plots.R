#' Make default plot array associated with processing step.
#'
#' Creates a number of plots and csv files relating to cluster composition
#' by all metadata in the ser object. This includes dim plots, feature plots,
#' and pie plots for cluster composition. Point size is automatically scaled
#' based on total number of cells and 10%/90% quantiles are used for all
#' feature plots. Reductions used are umap and phate.
#'
#' @param ser           Seurat object to process.
#' @param out_dir       Output directory
#' @param colors        Optional color schemes for metadata. Each scheme (for a given metadata) must be a vector of 
#'                      colors named with the levels of the metadata. The colors object is a named list of these named vectors. 
#'                      The name of the scheme must match the metadata name in the seurat object.
#' @param use_phate     Boolean indicating whether to use phate dimensional reduction for plots.  
#' @param cluster_enrichment_pvals  Boolean indicating whether to calculate P values for clusters by metadata.
#' @return Makes plots in user inputed directory
#' @import grDevices
#' @export

make_processing_plots <- function(ser, out_dir = '1_process/', colors = NULL, 
    use_phate = TRUE, cluster_enrichment_pvals = FALSE){

    dir.create(out_dir)
    if(ncol(ser) < 5000){
        pt.size = 2
    } else if(ncol(ser) < 10000){
        pt.size = 1.5
    } else {
        pt.size = 1
    }

    fp_cols = c("#2c7bb6", "#00a6ca", "#00ccbc", "#90eb9d", "#ffff8c", "#f9d057",
        "#f29e2e", "#e76818", "#d7191c")
    cl_cols = colors[['seurat_clusters']]

    for(meta in colnames(ser[[]])){
        print(meta)
        if(class(ser[[meta]][,1]) == 'factor'){
            # Pull out colors for metadata if it exists.
            id = which(names(colors) == meta)
            if(length(id)==1){
                cols = colors[[meta]]
            } else {
                cols = NULL
            }
            if(length(id)>1){
                warning('Found multiple color schemes for a metadata - ignoring colors.')
            }

            # Color dim plots by meta
            png(paste0(out_dir, '/', meta, '_umap_nolabel.png'), 
                height = 1000, width = 1000)
            print(Seurat::DimPlot(ser, group.by = meta, cols = cols, pt.size = pt.size))
            dev.off()

            png(paste0(out_dir, '/', meta, '_umap.png'), 
                height = 1000, width = 1000)
            print(Seurat::DimPlot(ser, group.by = meta, cols = cols, pt.size = pt.size, label = TRUE, label.size = 8))
            dev.off()

            if(use_phate){
                png(paste0(out_dir, '/', meta, '_phate_nolabel.png'), 
                    height = 1000, width = 1000)
                print(Seurat::DimPlot(ser, group.by = meta, cols = cols, pt.size = pt.size, reduction = 'phate'))
                dev.off()

                png(paste0(out_dir, '/', meta, '_phate.png'), 
                    height = 1000, width = 1000)
                print(Seurat::DimPlot(ser, group.by = meta, cols = cols, pt.size = pt.size, label = TRUE, label.size = 8, reduction = 'phate'))
                dev.off()
            }

            # Split dim plots by meta
            n_factor = length(levels(ser[[meta]][,1]))
            ncol = ceiling(sqrt(n_factor))
            nrow = ceiling(n_factor/ncol)
            h = 500*nrow
            w = 500*ncol
            png(paste0(out_dir, '/', meta, '_umap_nolabel_split.png'), 
                height = h, width = w)
            print(Seurat::DimPlot(ser, split.by = meta, cols = cl_cols, pt.size = pt.size, ncol = ncol))
            dev.off()

            png(paste0(out_dir, '/', meta, '_umap_split.png'), 
                height = h, width = w)
            print(Seurat::DimPlot(ser, split.by = meta, cols = cl_cols, pt.size = pt.size, ncol = ncol, label.size = 8, label = TRUE))
            dev.off()

            if(use_phate){
                png(paste0(out_dir, '/', meta, '_phate_nolabel_split.png'), 
                    height = h, width = w)
                print(Seurat::DimPlot(ser, split.by = meta, cols = cl_cols, pt.size = pt.size, ncol = ncol, reduction = 'phate'))
                dev.off()

                png(paste0(out_dir, '/', meta, '_phate_split.png'), 
                    height = h, width = w)
                print(Seurat::DimPlot(ser, split.by = meta, cols = cl_cols, pt.size = pt.size, ncol = ncol, label.size = 8, label = TRUE, reduction = 'phate'))
                dev.off()
            } 

            clust_proportions(ser, meta, 
                paste0(out_dir, '/', meta, '_proportions.pdf'),
                paste0(out_dir, '/', meta, '_proportions.csv'),
                cols, cluster_enrichment_pvals)           
        }
        if(class(ser[[meta]][,1]) == 'numeric'){
            png(paste0(out_dir, '/', meta, '_umap_feature.png'),
                height = 1000, width = 1000)
            print(Seurat::FeaturePlot(ser, features = meta, pt.size = pt.size, 
                cols = fp_cols, order = TRUE))
            dev.off()

            png(paste0(out_dir, '/', meta, '_umap_feature_q10.png'),
                height = 1000, width = 1000)
            print(Seurat::FeaturePlot(ser, features = meta,  pt.size = pt.size, 
                cols = fp_cols, order = TRUE, min.cutoff = 'q10',
                max.cutoff = 'q90'))
            dev.off()

            if(use_phate){
                png(paste0(out_dir, '/', meta, '_phate_feature.png'),
                    height = 1000, width = 1000)
                print(Seurat::FeaturePlot(ser, features = meta, pt.size = pt.size, 
                    cols = fp_cols, order = TRUE, reduction = 'phate'))
                dev.off()

                png(paste0(out_dir, '/', meta, '_phate_feature_q10.png'),
                    height = 1000, width = 1000)
                print(Seurat::FeaturePlot(ser, features = meta,  pt.size = pt.size, 
                    cols = fp_cols, order = TRUE, reduction = 'phate', 
                    min.cutoff = 'q10', max.cutoff = 'q90'))
                dev.off()
            }
        }
    }
}
