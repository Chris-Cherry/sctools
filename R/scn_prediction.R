#' singleCellNet training and prediction
#'
#' Trains a singleCellNet classifier and then predicts classification of query data
#' 
#' @param query     The data which will be classified by the random forest classifier
#' @param stRef     Reference sample data. Needs to have a column named "cell" with cell names and one 
#'                  named "newAnn" with annotations/desired clasifications. 
#' @param expRef    Reference experiment data. Should contain genes as row names, cell names as column
#'                  names and contain gene expression levels in the matrix
#'
#' @return          Returns a list with the trained classifier, predictions for the queried data, 
#'                  and assessment information from classifier training.
#' @export
#'

scn_prediction <- function(query, stRef, expRef){
    # Extract query data
    seuratfile = singleCellNet::extractSeurat(query, exp_slot_name = "counts")
    sampTab = seuratfile$sampTab
    expDat = query@assays$RNA@counts
    stQuery = sampTab
    expQuery = expDat
    genesQuery = rownames(expQuery)
    
    # Limit analysis to common genes
    commonGenes = intersect(rownames(expRef), genesQuery)
    expRef = expRef[commonGenes,]

    # Split for training and testing
    stList = singleCellNet::splitCommon(sampTab = stRef, ncells = 100, dLevel = "newAnn")
    stTrain = stList[[1]]
    expTrain = expRef[,rownames(stTrain)]

    # Train the classifier  
    topgenes = 10
    toppairs = 25
    system.time(class_info <- singleCellNet::scn_train(stTrain = stTrain, expTrain = expTrain,
        nTopGenes = topgenes, nRand = 70, nTrees = 1000, nTopGenePairs = toppairs,
        dLevel = "newAnn", colName_samp = "cell"))

    # Test
    stTestList = singleCellNet::splitCommon(sampTab = stList[[2]], ncells = 100, dLevel = "newAnn")
    stTest = stTestList[[1]]
    expTest = expRef[commonGenes, rownames(stTest)]
    classRes_val_all = singleCellNet::scn_predict(cnProc=class_info[['cnProc']], expDat=expTest, nrand = 50)

    # Assess the classifier
    tm_heldoutassessment = singleCellNet::assess_comm(ct_scores = classRes_val_all, stTrain = stTrain, 
        stQuery = stTest, dLevelSID = "cell", classTrain = "newAnn", classQuery = "newAnn", nRand = 50)

    # Query our data
    nqRand = 50 #Number of random profiles generated for evaluation
    system.time(cnPred<-singleCellNet::scn_predict(class_info[['cnProc']], expQuery, nrand=nqRand))

    things_we_want = list()
    things_we_want[["classifier"]] = class_info
    things_we_want[["predictions"]] = cnPred
    things_we_want[["assessment"]] = tm_heldoutassessment
    return(things_we_want)
}
