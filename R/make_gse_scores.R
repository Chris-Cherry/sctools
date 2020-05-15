#' Makes gene set enrichment scores by cell
#'
#' This function will read in the gene sets from an existing list, gmt file, 
#' or csv file and score cells based on scaled gene expression. The results 
#' will be appended to the set_scores assay in the provided Seurat object or 
#' the assay will be created.
#'
#' @param directory     Directory to get gmt file
#' @param gene_sets    User defined subset of gene set (Optional) 
#' @param ser           Seurat object to process
#' @param csv_dir       Directory to get csv file (Optional)
#' @param type          Output score in Absolute or Real values
#' @param from_gene     "ENSG" or "ENSMUSG" or "HGNC" or "MGI"             
#' @param to_gene       "MGI" or "HGNC"
#' @importFrom fgsea gmtPathways
#' @importFrom dplyr %>%
#' @importFrom Matrix colMeans
#' @importFrom utils read.table
#' @importFrom stats na.omit
#' @return  Outputs seurat object
#' @export

make_gse_scores <- function(ser, directory = NULL, from_gene = 'HGNC', 
    to_gene = 'MGI', gene_sets = NULL, csv_dir = NULL, type = 'Real'){
    
    gse = list()
    # Read in gene sets
    if (!is.null(gene_sets)){
        gse = c(gse, gene_sets)
    } 
    if (!is.null(directory)){
        gse = fgsea::gmtPathways(directory)
    } 
    if (!is.null(csv_dir)){
        table = read.table(csv_dir, sep = ',', header = TRUE, stringsAsFactors = FALSE)
        for(col in 1:ncol(table)){
            gse[[table[1,col]]] = table[-1, col]
        }
    }

    # Get scale data and convert genes if necessary
    scale_data = Seurat::GetAssayData(ser, slot = 'scale.data')
    if(to_gene != from_gene){
        gene_conv = convert_genes(rownames(scale_data), to_gene, from_gene)
        for(set in names(gse)){
            g = gse[[set]]
            mid = match(g, gene_conv[,2])
            mid = na.omit(mid)
            gse[[set]] = gene_conv[mid, 1]
        }
    } else {
        for(set in names(gse)){
            g = gse[[set]]
            mid = match(g, rownames(scale_data))
            mid = na.omit(mid)
            gse[[set]] = rownames(scale_data)[mid]
        }
    }

    # Score gene sets
    if (length(which(duplicated(gse))) > 0){
        dup = which(duplicated(gse))
        gse = gse[-dup]
    }
    set_scores = matrix(0, ncol = ncol(scale_data), nrow = length(gse))
    rownames(set_scores) = names(gse)
    colnames(set_scores) = colnames(scale_data) 
    for(set in names(gse)){
        if(length(gse[[set]]) == 1){next}
        m = scale_data[gse[[set]],]
        scores = Matrix::colMeans(m)
        set_scores[set, names(scores)] = scores
    }

    if (type != "Real" && type == "Abs")
        {
            set_scores = abs(set_scores)
        }
        else{
            print("Default data type is Real values")
        }

    # Add scores to Seurat object
    if(sum(names(ser@assays) == 'set_scores')){
        set_scores = rbind(ser@assays$set_scores@counts, set_scores)
    }
    ser[['set_scores']] = Seurat::CreateAssayObject(set_scores)
    return(ser)
}  
