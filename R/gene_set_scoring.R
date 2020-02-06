#' Scores cell expression of genes in set
#'
#' Takes in a gene set and uses normalized gene expression data to calculate a score describing extent to which
#' cells are expressing genes in the set.
#'
#' Parameters:
#' @param geneset   A list of genes and their fold changes. Should be organized such that "genes" are genes and "FC" is (log)foldchange. All genes will be assumed to be significant.
#' @param ser       A Seurat object to be scored by the gene set. Must contain scaled data
#' @return          Outputs a vector of score values named as cell names
#' @export

gene_set_scoring <- function(ser, geneset){

    # Separate genes into positive and negative changes    
    gene_list = list()
    gene_list[["all"]] = geneset$gene
    gene_list[["positive"]] = geneset$gene[which(geneset$FC>0)]
    gene_list[["negative"]] = geneset$gene[which(geneset$FC<0)]

    # Calculate summed z-scores, taking into account directionality of fold change
    pos = gene_list[["positive"]]
    neg = gene_list[["negative"]]
    pos_ind = match(pos, rownames(ser@assays$RNA@scale.data))
    pos_ind = pos_ind[-which(is.na(pos_ind))]
    pos_gene_subset = ser@assays$RNA@scale.data[pos_ind,]
    pos_gene_neg_score = pos_gene_subset
    pos_gene_neg_score[which(pos_gene_subset > 0)] = 0
    pos_gene_pos_score = pos_gene_subset
    pos_gene_pos_score[which(pos_gene_subset < 0)] = 0
    neg_ind = match(neg, rownames(ser@assays$RNA@scale.data))
    neg_ind = neg_ind[-which(is.na(neg_ind))]
    neg_gene_subset = ser@assays$RNA@scale.data[neg_ind,]
    neg_gene_neg_score = neg_gene_subset
    neg_gene_neg_score[which(neg_gene_subset > 0)] = 0
    neg_gene_pos_score = neg_gene_subset
    neg_gene_pos_score[which(neg_gene_subset < 0)] = 0
    scores = abs(colSums(neg_gene_neg_score))+abs(colSums(pos_gene_pos_score))-abs(colSums(neg_gene_pos_score))-abs(colSums(pos_gene_neg_score))
    names(scores) = colnames(ser)
    scores = scores/length(gene_list[['all']])
    return(scores)
}

