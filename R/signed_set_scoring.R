#' Scores cell expression of genes in set taking into account sign of fold change
#'
#' Takes in a gene set and uses normalized scaled gene expression data to calculate a score describing extent to which
#' cells are expressing genes in the set. Seurat object must contain scale_data assay. Gene set must be organized
#' such that "genes" are genes and "FC" is log(foldchange). All genes will be assumed to be significant.
#'
#' @param geneset   A list of genes and their fold changes. 
#' @param ser       A Seurat object to be scored by the gene set. Must contain scaled data
#' @param from_gene Original gene type
#' @param to_gene   Gene type for conversion
#' @param scaled    Boolean to determine whether to scale score by the number of genes in the set
#' @importFrom stats na.omit
#' @return          Outputs a vector of score values named as cell names
#' @export

signed_set_scoring <- function(ser, geneset, from_gene = "MGI", to_gene = "MGI", scaled = TRUE){
    
    # Gene conversion if required
    if(to_gene != from_gene){
        gene_conv = convert_genes(geneset$genes, from_gene, to_gene)
        mid = match(geneset$genes, gene_conv[,1])
        geneset$genes = gene_conv[mid, 2]
        geneset = na.omit(geneset)
        }
    
    # Calculate summed z-scores, taking into account directionality of fold change
    pos = geneset$genes[which(geneset$FC>0)]
    neg = geneset$genes[which(geneset$FC<0)]
    scale_dat = GetAssayData(object = ser, slot = "scale.data")
    pos_ind = match(pos, rownames(scale_dat))
    if (!identical(which(is.na(pos_ind)), integer(0))){
        pos_ind = pos_ind[-which(is.na(pos_ind))]}
    pos_gene_subset = scale_dat[pos_ind,]
    neg_ind = match(neg, rownames(scale_dat))
    if (!identical(which(is.na(neg_ind)), integer(0))){
        neg_ind = neg_ind[-which(is.na(neg_ind))]}
    neg_gene_subset = scale_dat[neg_ind,]
    if (is.null(dim(neg_gene_subset))){
        neg_score = neg_gene_subset * -1
    } else {
        neg_score = colSums(neg_gene_subset*-1)
    }
    if (is.null(dim(pos_gene_subset))){
        pos_score = pos_gene_subset
    } else {
        pos_score = colSums(pos_gene_subset)
    }
    scores = neg_score + pos_score
    names(scores) = colnames(ser)
    
    # Scale by number of genes
    if (scaled) {scores = scores/length(geneset$genes)}
    
    return(scores)
}
