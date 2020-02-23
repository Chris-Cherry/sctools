#' Parameters:
#' 
#' @param directory     Directory to get gmt file
#' @param rank_data     This can be a Directory to get RDS file or directly passing in dataframe
#' @param out_dir       Output directory
#' @param rank_by       User defined data ranking by logFC, Pvalue, or sign_Pvalue (Default setting at logFC)
#' @param from_gene     "ENSG" or "ENSMUSG" or "HGNC" or "MGI"             
#' @param cluster_name  Cluster name that user want to do the geneset analysis on
#' @param to_gene       "MGI" or HGNC
#' 
#' This function will perform fast preranked gene set encrichment analysis
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
            tmp = convert_genes(fgsea_pathway[[i]], from = from_gene, to = to_gene)
            tmp = tmp[,2]
            # Remove empty gene
            tmp = tmp[tmp != ""]
            # Find all the unique gene
            tmp = unique(tmp)
            fgsea_pathway[[i]] = tmp
        }

        names(rank) <- gse_sub$gene
        rank <- rank %>% mutate(rank = rank(avg_logFC, ties.method = "random"))
        print("Starting geneset analysis")
        fgseaResult = fgsea::fgsea(pathways = fgsea_pathway, stats = rank, nproc = 12, nperm = 10000)
    }
    else{
        print("Please input the correct inputs: either directory of RDS file or dataframe")
    }

    return(fgseaResult)
}
