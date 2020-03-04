#' Scores cell expression of genes in set
#'
#' Takes in a gene set and uses normalized gene expression data to calculate a score describing extent to which
#' cells are expressing genes in the set.
#'
#' Parameters:
#' @param geneset   A list of genes and their fold changes. Should be organized such that "genes" are genes and "FC" is (log)foldchange. All genes will be assumed to be significant.
#' @param ser       A Seurat object to be scored by the gene set. Must contain scaled data
#' @param scaled    Boolean to determine whether to scale score by the number of genes in the set
#' @return          Outputs a vector of score values named as cell names
#' @export

gene_set_scoring <- function(ser, geneset, scaled = TRUE){

    # Separate genes into positive and negative changes    
    gene_list = list()
    gene_list[["all"]] = geneset$genes
    gene_list[["positive"]] = geneset$genes[which(geneset$FC>0)]
    gene_list[["negative"]] = geneset$genes[which(geneset$FC<0)]

    # Calculate summed z-scores, taking into account directionality of fold change
    pos = gene_list[["positive"]]
    neg = gene_list[["negative"]]
    pos_ind = match(pos, rownames(ser@assays$RNA@scale.data))
    if (!identical(which(is.na(pos_ind)), integer(0))){
        pos_ind = pos_ind[-which(is.na(pos_ind))]}
    pos_gene_subset = ser@assays$RNA@scale.data[pos_ind,]
    neg_ind = match(neg, rownames(ser@assays$RNA@scale.data))
    if (!identical(which(is.na(neg_ind)), integer(0))){
        neg_ind = neg_ind[-which(is.na(neg_ind))]}
    neg_gene_subset = ser@assays$RNA@scale.data[neg_ind,]
    scores = colSums(neg_gene_subset*-1)+colSums(pos_gene_subset)
    names(scores) = colnames(ser)
    if (scaled) {scores = scores/length(gene_list[['all']])}
    return(scores)
}

