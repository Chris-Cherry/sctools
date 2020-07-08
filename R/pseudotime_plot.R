#' Creates pseudotime
#' 
#' @param ser           Seurat object to process
#' @param out_dir       Output directory
#' @param subset_gse    User defined subset of gene set (Optional) 
#' @param cluster_start Indicate what cluster to start with
#' @param cluster_end   Indicate what cluster to end with
#' @param reduction     Reduction method
#' @param lineage_ls    User customized list of cluster that want to be perform slingshot on 
#' @import Seurat
#' @import slingshot
#' @import RColorBrewer
#' 
#' @export

pseudotime_plot <- function(ser, out_dir = '4_pseudotime', cluster_start = NULL, cluster_end = NULL, 
    cal_reduction = 'phate', plot_reduction = 'phate', group_by = NULL, lineage_ls = NULL){
    
    dir.create(out_dir)

    reads <- Embeddings(ser, cal_reduction)
    clust <- Idents(object = ser)  
    
    if (!is.null(group_by)){
        subser = subset(ser, idents = group_by)
        clust <- Idents(object = subser)
        reads <- Embeddings(subser, cal_reduction)
    }
    else{
        subser = ser
    }

    tmp = list()

    if (!is.null(lineage_ls)){
        for (num_ls in 1:length(lineage_ls)){
            subser = subset(ser, idents = lineage_ls[[num_ls]])
            clust <- Idents(object = subser)
            reads <- Embeddings(subser, cal_reduction)
            lin = getLineages(reads, clust, reducedDim = cal_reduction, start.clus = cluster_start, end.clus = cluster_end)
        
            reads <- Embeddings(subser, plot_reduction)

            # Generate slingshot lineages
            lin_plot = getLineages(reads, clust, reducedDim = plot_reduction, start.clus = cluster_start, end.clus = cluster_end) 
            png(paste0(out_dir, '/slingshot_lineages', num_ls,'.png'), height = 2000, width = 2000)
            colcode = brewer.pal(9,'Set1')[Idents(subser)]
            lg_colcode = brewer.pal(9,'Set1')[unique(Idents(subser))]
            plot(reads, col = colcode)
            lines(lin_plot, type ='l', lwd = 3, col = 'black')
            legend("right", legend = unique(clust), text.font = 20, fill = lg_colcode)
            dev.off()
            
            # Generate slingshot curves
            curv = getCurves(lin_plot)
            # Plot curve slinghot in 2D phate format
            png(paste0(out_dir, '/slingshot_curves', num_ls,'.png'), height = 2000, width = 2000)
            plot(reads, col = colcode)
            lines(curv, lwd = 3, col = 'black')
            legend("right", legend = unique(clust), text.font = 20, fill = lg_colcode)
            dev.off()

            # How to make orthogonal projection onto the 2D plot

            # Create input for 3D plot
            reads <- Embeddings(subser, "phate3d")
            lin_3dplot = getLineages(reads, clust, reducedDim = plot_reduction, start.clus = cluster_start, end.clus = cluster_end) 
            curv_3d = getCurves(lin_3dplot)

            saveRDS(curv_3d@curves, paste0(out_dir, '/curves_lin_', num_ls,'.RDS'))
            
            spt = slingPseudotime(curv_3d)
            tmp = spt[,]
            spt_tmp = list()
            if (is.null(ncol(tmp))){
                spt_tmp = list(tmp)
            }
            else{
                for (j in 1:ncol(tmp)){
                    spt_tmp[j] <- list(spt[,j])
                }
            }

            saveRDS(spt_tmp, paste0(out_dir, '/pseudotime_lin_', num_ls, '.RDS'))
        }
    }
    else{

        lin = getLineages(reads, clust, reducedDim = cal_reduction, start.clus = cluster_start, end.clus = cluster_end)
   
        reads <- Embeddings(subser, plot_reduction)
        # Generate slingshot lineages
        lin_plot = getLineages(reads, clust, reducedDim = plot_reduction, start.clus = cluster_start, end.clus = cluster_end) 
        png(paste0(out_dir, '/slingshot_lineages_all.png'), height = 2000, width = 2000)
        colcode = colorRampPalette(brewer.pal(9,'Set1'))(20)[Idents(subser)]
        lg_colcode = colorRampPalette(brewer.pal(9,'Set1'))(20)[unique(Idents(subser))]
        plot(reads, col = colcode)
        lines(lin_plot, type ='l', lwd = 3, col = 'black')
        legend("right", legend = unique(clust), text.font = 20, fill = lg_colcode)
        dev.off()
        
        # Generate slingshot curves
        curv = getCurves(lin_plot)
        # Plot curve slinghot in 2D phate format
        png(paste0(out_dir, '/slingshot_curves_all.png'), height = 2000, width = 2000)
        plot(reads, col = colcode)
        lines(curv, lwd = 3, col = 'black')
        legend("right", legend = unique(clust), text.font = 20, fill = lg_colcode)
        dev.off()

        # Create inputs for 3D plot
        reads <- Embeddings(ser, "phate3d")
        lin_3dplot = getLineages(reads, clust, reducedDim = plot_reduction, start.clus = cluster_start, end.clus = cluster_end) 
        curv_3d = getCurves(lin_3dplot)
        saveRDS(curv_3d@curves, paste0(out_dir, '/curves_lin_all.RDS'))
        
        spt = slingPseudotime(curv_3d)
        tmp = spt[,]
        spt_tmp = list()
        if (is.null(ncol(tmp))){
            spt_tmp = list(spt)
        }
        else{
            for (j in 1:ncol(tmp)){
                spt_tmp[j] <- list(spt[,j])
            }
        }
        saveRDS(spt_tmp, paste0(out_dir, '/pseudotime_lin_all.RDS'))

    }

}
