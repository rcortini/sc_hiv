# load the WGCNA library and allow multithreading
suppressMessages(library(WGCNA))
suppressMessages(allowWGCNAThreads())
suppressMessages(library(biomaRt))

PrepareDatExpr <- function (exprMatrix, sampleSheet,
                          ngenes,
                          cut) {
    # select the group of genes from the untreated J-Lat cells
    jlat.untreated <- exprMatrix[, sampleSheet$label == 'J-Lat+DMSO']
    
    # establish which are the most highly varying genes, based on a simple
    # criterion of maximum variance/mean.
    gene.variances <- apply(jlat.untreated, 1, var)
    gene.means <- apply(jlat.untreated, 1, mean)
    gene.variability <- gene.variances/gene.means
    
    # get the names of the genes that have the greatest biological variation, 
    # excluding the FILIONG01 gene (not really necessary)
    selected <- order(gene.variability, decreasing = TRUE)[1:ngenes]
    most.variable.genes <- rownames(jlat.untreated[selected, ])
    most.variable.genes <- most.variable.genes[most.variable.genes != 'FILIONG01']
    
    # extract a data frame with the values of the expressions for each of the genes
    # with the highest biological variation
    datExpr0 <- as.data.frame(t(jlat.untreated[most.variable.genes, ]))
    
    # do quality control
    gsg <- goodSamplesGenes(datExpr0, verbose = 3);
    if (!gsg$allOK) {
        stop("Do proper quality control on genes!") 
    }
    
    # extract the hierarchical clustering tree of the samples
    sampleTree <- hclust(dist(datExpr0), method = "average");

    # cut the tree according to the user-supplied `cut` parameter, and then 
    # clust 1 will contains the samples we want to keep.
    clust <- cutreeStatic(sampleTree, cutHeight = cut, minSize = 10)
    keepSamples <- (clust == 1)

    # if everything is okay, define a new data expression data.frame
    datExpr0[keepSamples, ]
}

ClusterGenes <- function (datExpr, softThresholdPower) {
    # reconstruct the network
    net <- blockwiseModules(datExpr,
			    power             = softThresholdPower,
			    TOMType           = "unsigned",
			    inModuleSize      = 30,
			    reassignThreshold = 0,
			    mergeCutHeight    = 0.25,
			    numericLabels     = TRUE,
			    pamRespectsDendro = FALSE,
			    verbose           = 0)
    net
}

AssociateClustersToHIV <- function (datExpr, exprMatrix, sampleSheet, net,
				    aliveThreshold) {
    # get the module labels, transform them into colors
    moduleColors <- labels2colors(net$colors)

    # get the names of the genes we selected from the original ones
    myGenes <- colnames(datExpr)

    # select only the genes that we selected before
    myExprMatrix <- exprMatrix[myGenes, ]

    # select only J-Lat treated cells
    myExprMatrix <- myExprMatrix[, sampleSheet$label == "J-Lat+SAHA"]

    # select only alive cells
    myExprMatrix <- myExprMatrix[, colSums(myExprMatrix) > aliveThreshold]

    # finally, transpose to be interfaced to WGCNA
    myExprMatrix <- t(myExprMatrix)

    # get the module eigengenes of the *new* data set: that is, we assign the
    # expression profiles of the treated data set based on the gene modules of the
    # untreated cells
    MEs <- moduleEigengenes(myExprMatrix, moduleColors)$eigengenes
    MEs <- orderMEs(MEs)

    # get the names of the cells that we have selected, and extract the HIV profile
    # of those cells
    myCells <- rownames(myExprMatrix)
    hiv <- t(exprMatrix["FILIONG01", myCells])

    # parameters of our data set
    nGenes <- ncol(myExprMatrix)
    nSamples <- nrow(myExprMatrix)

    # correlate the module eigengenes to the HIV expression patterns, and 
    # calculate the corresponding p value
    moduleHivCor <- cor(MEs, hiv, use = "p")
    moduleHivPvalue <- corPvalueStudent(moduleHivCor, nSamples)

    # prepare the return data structure `module`
    modules <- list()

    # add the information on the module eigengenes, together with the hiv
    # expression associated to each cell
    modules[["MEs"]] <- data.frame(MEs)
    rownames(modules[["MEs"]]) <- myCells
    modules[["MEs"]]$hiv <- hiv

    # add the statistics associated to the module eigengenes - to - HIV
    # correlation
    modules[["stats"]] <- data.frame(cor = moduleHivCor, p = moduleHivPvalue)
    names(modules[["stats"]]) <- c("cor", "p")
    modules
}
