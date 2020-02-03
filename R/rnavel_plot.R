#' Creates an RNA velocity projection
#' 
#' The function takes a loom file and Seurat object and creates a projection of
#' RNA velocity onto a dimensional reduction of your choice. The dimensional
#' reduction should be present in the Seurat object.
#' 
#' @param loom List of loom files as read in from process_loom.R.
#' @param ser Seurat object containing clusters and dimensional reduction.
#' @param out_file File name of png file to save plot.
#' @param dr Name of dimensional reduction slot to use from Seurat object. Defaults to phate.
#' @param cols Named vector of colors for cell clusters. If NULL then ggplot colors are generated.
#' @param n_core Number of cores to use for velocyto.R functions. If NULL then future::availableCores will be used.
#' @param plot_info If RNA velocity was run previously you can provide the output file to skip much of the processing time. This can be very helpful when tweaking plotting parameters.
#' @param ... Parameters to pass to velocyto.R::show.velocity.on.embedding.cor
#'
rnavel_plot <- function(loom, ser, out_file, dr = 'phate', cols = NULL, n_core = NULL, plot_info = NULL, ...){
    if(is.null(n_core)){
        n_core = future::availableCores()
    }
    if(is.null(cols)){
        cols = ggplot_col_gen(length(levels(Seurat::Idents(ser))))
        names(cols) = levels(Seurat::Idents(ser))
    }
    cell.colors = cols[Seurat::Idents(ser)]
    names(cell.colors) = names(Seurat::Idents(ser))

    if(is.null(plot_info)){
        tar_cells = intersect(colnames(loom$spliced), 
            colnames(ser@assays$RNA@counts))
        spliced_counts = loom$spliced[, tar_cells]
        unspliced_counts = loom$unspliced[, tar_cells]

        idents = Seurat::Idents(ser)[tar_cells]
        pcs = ser@reductions$pca@cell.embeddings[tar_cells,]
        cell_dist <- as.dist(1-velocyto.R::armaCor(t(pcs)))
        emb = ser[[dr]]@cell.embeddings[tar_cells, 1:2]

        spliced = velocyto.R::filter.genes.by.cluster.expression(spliced_counts, 
            idents, min.max.cluster.average = 0.2)
        unspliced = velocyto.R::filter.genes.by.cluster.expression(
            unspliced_counts, idents, min.max.cluster.average = 0.05)

        vel_estimates = velocyto.R::gene.relative.velocity.estimates(spliced, 
            unspliced, n.cores = n_core, kCells = 40, 
            fit.quantile = 0.01, cell.dist = cell_dist)

        png(out_file, height = 1000, width = 1000)
        plot_out = velocyto.R::show.velocity.on.embedding.cor(emb, vel_estimates, 
            cell.colors = cell.colors, n.cores = n_core, return.details = TRUE, ...)
        dev.off()

        return(list(emb, vel_estimates, plot_out))
    } else {
        emb = plot_info[[1]]
        vel_estimates = plot_info[[2]]
        cc = plot_info[[3]]$cc
        png(out_file, height = 1000, width = 1000)
        velocyto.R::show.velocity.on.embedding.cor(emb, vel_estimates, 
            cell.colors = cell.colors, n.cores = n_core, 
            return.details = FALSE, cc = cc, ...)
        dev.off()
        return()
    }
}

