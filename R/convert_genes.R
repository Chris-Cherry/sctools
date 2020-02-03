convert_genes <- function(genes, from, to){
    if (from == 'ENSMUSG'){
        srcMart = biomaRt::useMart("ensembl", dataset = "mmusculus_gene_ensembl")
        sourceAtts =  'ensembl_gene_id'   
    }
    if (from == 'ENSG'){
        srcMart = biomaRt::useMart("ensembl", dataset = "hsapiens_gene_ensembl")
        sourceAtts = 'ensembl_gene_id'
    }
    if (from == 'MGI'){
        srcMart = biomaRt::useMart("ensembl", dataset = "mmusculus_gene_ensembl")
        sourceAtts = 'mgi_symbol'    
    }
    if (from == 'HGNC'){
        srcMart = biomaRt::useMart("ensembl", dataset = "hsapiens_gene_ensembl")
        sourceAtts = 'hgnc_symbol'
    }
    if (to == 'MGI'){
        tarMart = biomaRt::useMart("ensembl", dataset = "mmusculus_gene_ensembl")
        tarAtts = 'mgi_symbol'
    }
    if (to == 'HGNC'){
        tarMart = biomaRt::useMart("ensembl", dataset = "hsapiens_gene_ensembl")
        tarAtts = 'hgnc_symbol'
    }
    genesV2 = biomaRt::getLDS(attributes = sourceAtts, filters = sourceAtts,
                     values = genes, mart = srcMart, 
                     attributesL = tarAtts, martL = tarMart,
                     uniqueRows = F)
    return(genesV2)
}
