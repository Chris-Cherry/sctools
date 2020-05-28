#' This function will perform fast preranked gene set encrichment analysis
#' 
#' @param directory     Directory to get gmt file
#' @param rank_data     This can be a Directory to get RDS file or directly passing in dataframe
#' @param out_dir       Output directory
#' @param rank_by       User defined data ranking by logFC, Pvalue, or sign_Pvalue (Default setting at logFC)
#' @param from_gene     "ENSG" or "ENSMUSG" or "HGNC" or "MGI"             
#' @param cluster_name  Cluster name that user want to do the geneset analysis on
#' @param to_gene       "MGI" or HGNC
#' 
#' 
#' @return  Outputs dataframe from fgsea function
#' @export

rank_gse <- function(directory, rank_data, out_dir = '3_gse', from_gene, to_gene, 
    cluster_name, rank_by = 'logFC'){

    if (typeof(rank_data) == 'character' || typeof(rank_data) == 'list'){
        
        if (typeof(rank_data) == 'character'){
            gse_rank = readRDS(paste0(rank_data))
        }

        gse_sub = gse_rank[which(gse_rank$cluster == cluster_name),]

        if (rank_by == 'Pvalue'){
            rank = gse_sub$p_value
        }
        else if(rank_by == 'sign_Pvalue'){
            rank = gse_sub$p_value * gse_sub$avg_logFC
        }
        else{
            rank = gse_sub$avg_logFC
        }

        fgsea_pathway = fgsea::gmtPathways(directory)
        for (i in 1:length(fgsea_pathway)){
            # Convert genes
            if (from_gene != to_gene){
            tmp = convert_genes(fgsea_pathway[[i]], from = from_gene, to = to_gene)
            tmp = tmp[,2]
            } else {
                tmp = fgsea_pathway[[i]]
            }
            # Remove empty gene
            tmp = tmp[tmp != ""]
            # Find all the unique gene
            tmp = unique(tmp)
            fgsea_pathway[[i]] = tmp
        }

        names(rank) <- gse_sub$gene
        
        # Eliminate the Inf and -Inf in the rank
        if(max(rank) == Inf){
            tmp = which(rank == Inf)
            rank = rank[-tmp]
        }
         if(min(rank) == -Inf){
            tmp = which(rank == -Inf)
            rank = rank[-tmp]
        }

        fgseaResult = fgsea::fgsea(pathways = fgsea_pathway, stats = rank, nproc = 12, maxSize = 500, nperm = 10000)
    }
    else{
        print("Please input the correct inputs: either directory of RDS file or dataframe")
    }

    return(fgseaResult)
}
