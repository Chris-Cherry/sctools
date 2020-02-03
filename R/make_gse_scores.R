#' Parameters:
#' @param directory     Directory to get gmt file
#' @param out_dir       Output directory
#' @param subset_gse    User defined subset of gene set (Optional) 
#' @param ser           Seurat object to process
#' @param csv_dir       Directory to get csv file (Optional)
#' @param type          Output score in Absolute or Real values
#' @param from_gene     "ENSG" or "ENSMUSG" or "HGNC" or "MGI"             
#' @param to_gene       "MGI" or HGNC
#' 
#' This function will read in the genesets, either default or customized by users
#' and make feature plots from those genesets
#' 
#' @import dplyr
#' @importFrom dplyr %>%
#' @return  Outputs seurat object
#' @export

make_gse_scores <- function(ser, directory = NULL, out_dir = '3_gse', from_gene = 'HGNC', to_gene = 'MGI', 
    subset_gse = NULL, csv_dir = NULL, type = 'Real'){
    
    dir.create(out_dir)

    if (!is.null(subset_gse) || !is.null(csv_dir) || !is.null(directory)){
        if (!is.null(csv_dir)){
            # Red in an excel file
            gse_conv = read.table(csv_dir, sep = ',')
            for (i in 1:length(gse_conv)){
                tmp = gse_conv[i]
                tmp = tmp[tmp != ""]
                if (i == 1){
                    gse <- list(t(tmp))
                }
                else{
                    gse[i] <- list(t(tmp))
                }
            }
            
            # Extract the Seurat data
            gse_ser = ser@assays$RNA@scale.data
            gse_ser_row = rownames(gse_ser)
            gse_ser_col = colnames(gse_ser)

            for (i in 1:length(gse)){
                # Find unique genes
                gse = unique(gse)
                # Convert genes
                tmp = convert_genes(gse[[i]], from = from_gene, to = to_gene)
                tmp = tmp[,2]

                # Find index of each gene set
                idx_gse = match(tmp, gse_ser_row)

                idx_gse = na.omit(idx_gse)

                gse_ser_sub = gse_ser[idx_gse,]
                gse_ser_tmp <- Matrix::colSums(gse_ser_sub)
                if (i == 1){
                    gse_ser_set = gse_ser_tmp
                    set_name = paste0("GeneSet", i)
                }
                else{
                    gse_ser_set = rbind(gse_ser_set, gse_ser_tmp)
                    set_name[i] = paste0("GeneSet", i)
                }
            }
        }
        else if (!is.null(directory)){
            # Read in geneset with gmt pathway
            gse = fgsea::gmtPathways(directory)

            # Extract the Seurat data
            gse_ser = ser@assays$RNA@scale.data
            gse_ser_row = rownames(gse_ser)
            gse_ser_col = colnames(gse_ser)

            for (i in 1:length(gse)){
                # Remove empty cells
                gse = gse[gse != ""]
                # Find unique genes
                gse = unique(gse)
                # Convert genes
                tmp = convert_genes(gse[[i]], from = from_gene, to = to_gene)
                tmp = tmp[,2]

                # Find index of each gene set
                idx_gse = match(tmp, gse_ser_row)

                idx_gse = na.omit(idx_gse)

                gse_ser_sub = gse_ser[idx_gse,]
                gse_ser_tmp <- Matrix::colSums(gse_ser_sub)
                if (i == 1){
                    gse_ser_set = gse_ser_tmp
                    set_name = paste0("GeneSet", i)
                }
                else{
                    gse_ser_set = rbind(gse_ser_set, gse_ser_tmp)
                    set_name[i] = paste0("GeneSet", i)
                }
            }
        }
        else{
            if (typeof(subset_gse) =='list'){
            # Extract the Seurat data
            gse = subset_gse
            gse_ser = ser@assays$RNA@scale.data
            gse_ser_row = rownames(gse_ser)
            gse_ser_col = colnames(gse_ser)

            for (i in 1:length(gse)){
                # Remove empty cells
                gse = gse[gse != ""]
                # Find unique genes
                gse = unique(gse)
                # Convert genes
                tmp = convert_genes(gse[[i]], from = from_gene, to = to_gene)
                tmp = tmp[,2]

                # Find index of each gene set
                idx_gse = match(tmp, gse_ser_row)

                idx_gse = na.omit(idx_gse)

                gse_ser_sub = gse_ser[idx_gse,]
                gse_ser_tmp <- Matrix::colSums(gse_ser_sub)
                if (i == 1){
                    gse_ser_set = gse_ser_tmp
                    set_name = paste0("GeneSet", i)
                }
                else{
                    gse_ser_set = rbind(gse_ser_set, gse_ser_tmp)
                    set_name[i] = paste0("GeneSet", i)
                }
                }
            }
            else{
                print("Please insert the right data type (only except list data type)")
            }
        }

        rownames(gse_ser_set) = set_name

        if (type != "Real" && type == "Abs")
        {
            gse_ser_set = abs(gse_ser_set)
        }
        else{
            print("Default data type is Real values")
        }

        ser[['GeneSet']] = Seurat::CreateAssayObject(gse_ser_set)

        #Make feature plots
        tmp = ser
        DefaultAssay(tmp) = 'GeneSet'

        pdf(paste0(out_dir, '/FeaturePlot_GeneSet.pdf'), height = 20, width = 20)
        print(Seurat::FeaturePlot(tmp, features = rownames(tmp), ncol = 3))
    }

    else{
        print("Please specify directory of the geneset, or list of gene set that you want to analyse")
    }

    dev.off()

    return(ser)
}
