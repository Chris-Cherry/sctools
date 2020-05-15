#' Creates cluster proportions and plots from metadata.
#'
#' Creates a pdf of pie charts (one for each cluster with proportions from 
#' the given metadata) and a csv of raw cell numbers by cluster and metadata.
#' It also uses a permutation based test to determine statistical enrichment
#' of clusters by the given metadata. The stats are returned in the subtitles
#' of the pie charts.
#' @param ser           Seurat object with clusters and given metadata
#' @param meta          Name of metadata in ser object to get cluster proportions
#' @param csv           Name of output csv to save cell numbers by cluster/meta
#' @param pdf           Name of output pdf file
#' @param cluster_enrichment_pvals Boolean for whether to calculate pvals
#' @import grDevices
#' @importFrom utils write.table

clust_proportions = function(ser, meta, pdf, csv, cols, cluster_enrichment_pvals){

    # Start by generating null distribution for idents
    n = length(Seurat::Idents(ser))
    perms = 1000
    perm_idents = matrix(0, nrow = perms, ncol = n)
    colnames(perm_idents) = names(Seurat::Idents(ser))
    for(i in 1:perms){
        perm_idents[i,] = sample(Seurat::Idents(ser), n)
    }
    perm_idents = perm_idents - 1
    clust_comps = c()
    meta_totals = table(ser[[meta]])
    meta_totals[which(meta_totals == 0)] = 1
    # Origins of pooled clusters
    pdf(pdf, width = 18)
    for(clust in levels(Seurat::Idents(ser))){
        cell_ids = which(Seurat::Idents(ser) == clust)
        comps = table(ser[[meta]][cell_ids,])
        clust_comps = rbind(clust_comps, comps)

        # While we're here lets make pie charts
        # Normalize data to cell totals
        comps = comps/meta_totals
        temp = data.frame(comps)
        plot = ggplot2::ggplot(temp, ggplot2::aes(x="", y=Freq, fill=Var1)) +
            ggplot2::geom_bar(width=1, stat="identity") +
            ggplot2::labs(title = clust) +
            ggplot2::theme_void()
        if(cluster_enrichment_pvals){
            # Calculate P vals for the clust comps
            tables = matrix(0, ncol = length(levels(ser[[meta]][,1])), nrow = perms)
            for(i in 1:perms){
                cells = colnames(perm_idents)[which(perm_idents[i,] == clust)]
                table = table(ser[[meta]][cells,])
                tables[i,] = table
            }
            colnames(tables) = names(table)
            pvals = c()
            for(col in colnames(tables)){
                ngreater = sum(tables[,col] > comps[col])
                nless = sum(tables[,col] < comps[col])
                if(nless < ngreater){
                    p = nless*2/perms
                }else{
                    p = ngreater*2/perms
                }
                pvals[col] = p
            }
            pval_names = paste(names(pvals), pvals, sep = ':')
            plot = plot + ggplot2::labs(subtitle = paste(pval_names, collapse = ' '))
        }

        if(!is.null(cols)){
            plot = plot + ggplot2::scale_fill_manual(values = cols[names(comps)])
        } 

        # Pie chart:
        plot = plot + ggplot2::coord_polar('y', start=0)
        print(plot)
        rm(list = c('plot','temp','comps'))
    }
    dev.off()

    rownames(clust_comps) = levels(Seurat::Idents(ser))
    write.table(clust_comps, sep=',', col.names = NA,
                file=csv)
    rm(clust_comps)
}
