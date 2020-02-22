#' Runs differential expression and creates relevant plots for a ser object
#' Parameters:
#' @param ser           The Seurat object to use
#' @param list_ser      List of DE seurat object 
#' @param feats         A subset of featurs to run DE (I.E. only surface markers)
#'                      If NULL then will run on all genes
#' @param out_dir       Directory to write plots and save markers.
#' @param cols          Name vector of colors by cluster
#' @param feature_plots Whether or not to include feature plots in the pdf output.
#' @param meta          subset that you want to graph
#' @import dplyr
#' @import future
#' @importFrom dplyr %>%
#' @export

make_de_plots <- function(ser = NULL, list_ser = NULL, feats = NULL, out_dir = '2_de/', cols = NULL,
    feature_plots = TRUE, meta = NULL){

    future::plan(future::multiprocess)
    options(globals.future.maxSize = Inf)

    # Process each cluster
    if (!is.null(list_ser)){
        DefaultAssay(list_ser[[1]]) = 'RNA'
        for (i in meta){
            out_dir = paste0('2_de/cluster_', i, '/')
            submarks = readRDS(paste0(out_dir, '/submarks_', i, '.RDS'))
            top10 = submarks %>% group_by(cluster) %>% top_n(10, avg_logFC)
            if (any(submarks$avg_logFC > 1) == TRUE)
                top_all = submarks[which(submarks$avg_logFC > 1),]
            else {
                top_all = top10
            }       
            png(paste0(out_dir, '/heatmap_all_', i, '.png'), height = 2000, width = 2000)
            print(Seurat::DoHeatmap(list_ser[[i]], top_all$gene, assay = 'RNA', raster = FALSE, 
                draw.lines = FALSE))
            dev.off()

            png(paste0(out_dir, '/heatmap_trim_gene_', i, '.png'), height = 2000, width = 2000)
            print(Seurat::DoHeatmap(list_ser[[i]], top10$gene, assay = 'RNA', raster = FALSE, 
                draw.lines = FALSE))
            dev.off()

            smallest = min(table(Idents(list_ser[[i]])))
            subser = subset(list_ser[[i]], downsample = smallest*3)

            png(paste0(out_dir, '/heatmap_trim_cell_', i, '.png'), height = 2000, width = 2000)
            print(Seurat::DoHeatmap(subser, top_all$gene, assay = 'RNA', raster = FALSE, 
                draw.lines = FALSE))
            dev.off()

            png(paste0(out_dir, '/heatmap_trim_both_', i, '.png'), height = 2000, width = 2000)
            print(Seurat::DoHeatmap(subser, top10$gene, assay = 'RNA', raster = FALSE, 
                draw.lines = FALSE))
            dev.off()
            
            pdf(paste0(out_dir, '/clust_vlns_cluster_', i, '.pdf'), height = 20, width = 20)
            
            for(clust in levels(Seurat::Idents(list_ser[[i]]))){
                top = top10$gene[which(top10$cluster == clust)]
                print(Seurat::VlnPlot(list_ser[[i]], cols = cols, pt.size = 0, features = top, ncol = 3, 
                    assay = 'RNA'))
                if(feature_plots){
                print(Seurat::FeaturePlot(list_ser[[i]], features = top, ncol = 3))
                print(Seurat::FeaturePlot(list_ser[[i]], features = top, ncol = 3, min.cutoff = 'q10', max.cutoff = 'q90'))
                }
            }
        }
    }
    else if(!is.null(ser)){
    # Process the whole dataset
        DefaultAssay(ser) = 'RNA'
        marks = readRDS(paste0(out_dir, '/marks.RDS'))
        top10 = marks %>% group_by(cluster) %>% top_n(10, avg_logFC)
        top_all = marks[which(marks$avg_logFC > 1),]
        
        png(paste0(out_dir, '/heatmap_all.png'), height = 2000, width = 2000)
        print(Seurat::DoHeatmap(ser, top_all$gene, assay = 'RNA', raster = FALSE, 
            draw.lines = FALSE))
        dev.off()

        png(paste0(out_dir, '/heatmap_trim_gene.png'), height = 2000, width = 2000)
        print(Seurat::DoHeatmap(ser, top10$gene, assay = 'RNA', raster = FALSE, 
            draw.lines = FALSE))
        dev.off()

        smallest = min(table(Idents(ser)))
        subser = subset(ser, downsample = smallest*3)

        png(paste0(out_dir, '/heatmap_trim_cell.png'), height = 2000, width = 2000)
        print(Seurat::DoHeatmap(subser, top_all$gene, assay = 'RNA', raster = FALSE, 
            draw.lines = FALSE))
        dev.off()

        png(paste0(out_dir, '/heatmap_trim_both.png'), height = 2000, width = 2000)
        print(Seurat::DoHeatmap(subser, top10$gene, assay = 'RNA', raster = FALSE, 
            draw.lines = FALSE))
        dev.off()
        
        pdf(paste0(out_dir, '/clust_vlns.pdf'), height = 20, width = 20)
        
        for(clust in levels(Seurat::Idents(ser))){
            top = top10$gene[which(top10$cluster == clust)]
            print(Seurat::VlnPlot(ser, cols = cols, pt.size = 0, features = top, ncol = 3, 
                assay = 'RNA'))
            if(feature_plots){
            print(Seurat::FeaturePlot(ser, features = top, ncol = 3))
            print(Seurat::FeaturePlot(ser, features = top, ncol = 3, min.cutoff = 'q10', max.cutoff = 'q90'))
            }
        }
        dev.off()
    }
    else{
        print("Please input the seurat dataset, either input whole dataset or a list of clusters")
    }
}
