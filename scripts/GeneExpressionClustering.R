# load the WGCNA library and allow multithreading
suppressMessages(library(WGCNA))
suppressMessages(allowWGCNAThreads())
suppressMessages(library(biomaRt))

PrepareDataForClustering <- function(X, sampleSheet, ngenes, cut, plot = TRUE) {
    # select the group of genes from the untreated J-Lat cells
    jlat.untreated <- X[, sampleSheet$status == 'nontreated']

    # establish which are the most highly varying genes, based on a simple
    # criterion of maximum variance/mean.
    gene.variances <- apply(jlat.untreated, 1, var)
    gene.means <- apply(jlat.untreated, 1, mean)
    gene.variability <- gene.variances/gene.means

    # get the names of the genes that have the greatest biological variation, 
    selected <- order(gene.variability, decreasing = TRUE)[1:ngenes]
    most.variable.genes <- rownames(jlat.untreated[selected, ])

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

    # plot results if user wishes
    if (plot) {
      # plot size
      options(repr.plot.width = 10, repr.plot.height = 6)

      # detect outliers
      par(cex = 0.6);
      par(mar = c(0,4,2,0))
      plot(sampleTree,
	   main     = "Sample clustering to detect outliers",
	   sub      = "",
	   xlab     = "",
	   cex.lab  = 1.5,
	   cex.axis = 1.5,
	   cex.main = 2)

      # Plot a line to show the cut
      abline(h = cut, col = "red");
    }

    # cut the tree according to the user-supplied `cut` parameter, and then 
    # clust 1 will contains the samples we want to keep.
    clust <- cutreeStatic(sampleTree, cutHeight = cut, minSize = 10)
    keepSamples <- (clust == 1)

    # if everything is okay, define a new data expression data.frame
    datExpr0[keepSamples, ]
}

# this function outputs a plot that allows to choose the best value of the
# soft thresholding power
PrepareClustering <- function (datExpr) {
    # Choose a set of soft-thresholding powers
    powers <- c(c(1:10), seq(from = 12, to=20, by=2))

    # Call the network topology analysis function
    sft <- pickSoftThreshold(datExpr, powerVector = powers, verbose = 5)
    
    # number of genes and number of samples
    nGenes <- ncol(datExpr)
    nSamples <- nrow(datExpr)

    # Plot the results:
    par(mfrow = c(1,2))
    cex1 = 0.9
    options(repr.plot.width = 10, repr.plot.height = 6)

    # Scale-free topology fit index as a function of the soft-thresholding power
    plot(sft$fitIndices[,1], -sign(sft$fitIndices[,3])*sft$fitIndices[,2],
         xlab = "Soft Threshold (power)",
         ylab = "Scale Free Topology Model Fit,signed R^2",
         type = "n",
         main = paste("Scale independence"))

    text(sft$fitIndices[,1], -sign(sft$fitIndices[,3])*sft$fitIndices[,2],
         labels = powers,
         cex    = cex1,
         col    = "red");

    # this line corresponds to using an R^2 cut-off of h
    abline(h = 0.90, col = "red")

    # Mean connectivity as a function of the soft-thresholding power
    plot(sft$fitIndices[,1], sft$fitIndices[,5],
         xlab = "Soft Threshold (power)",
         ylab = "Mean Connectivity",
         type = "n",
         main = paste("Mean connectivity"))

    text(sft$fitIndices[,1], sft$fitIndices[,5],
         labels = powers,
         cex    = cex1,
         col    = "red")
}

VisualizeClustering <- function (net) {
    # plot size
    options(repr.plot.width = 10, repr.plot.height = 6)

    # Convert labels to colors for plotting
    mergedColors <- labels2colors(net$colors)

    # Plot the dendrogram and the module colors underneath
    plotDendroAndColors(net$dendrograms[[1]],
                        mergedColors[net$blockGenes[[1]]],
                        "Module colors",
                        dendroLabels = FALSE,
                        hang = 0.03,
                        addGuide = TRUE,
                        guideHang = 0.05)
}

GeneColors <- function (datExpr, net) {
    # get the module labels, transform them into colors
    colors <- labels2colors(net$colors)

    # get the names of the genes
    genes <- colnames(datExpr)
    
    # put the things together
    C <- data.frame(color = colors)
    rownames(C) <- genes
    
    # return
    C
}

ModuleAnalysis <- function (colors, exprMatrix, sampleSheet) {

    # here, "colors" is a data frame that has as row names the names of 
    # the genes that were use in the identification of the modules. Then
    # there is a column that is called "color" that represents what module that
    # particular gene belongs to.

    # select only the genes that we selected before, of the treated cells,
    # and then transpose the matrix to be given to WGCNA
    myExprMatrix <- t(exprMatrix[rownames(colors), sampleSheet$status == "treated"])

    # get the module eigengenes of the *new* data set: that is, we assign the
    # expression profiles of the treated data set based on the gene modules of the
    # untreated cells
    MEs <- moduleEigengenes(myExprMatrix, colors$color)$eigengenes
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

    # extract the names from the MEs (because they are MEgrey...)
    modNames <- substring(names(MEs), 3)

    # evaluate gene module membership, with associated p-values, and gene
    # to HIV correlations, together with p-values. The following two data frames are
    # full matrices: the row is the gene, the column is the module membership score,
    # and in the second one it is the p-value associated to belonging to that module.
    geneModuleMembership  <- cor(myExprMatrix, MEs, use = "p")
    MMPvalue              <- corPvalueStudent(as.matrix(geneModuleMembership), nSamples)
    colnames(geneModuleMembership) <- modNames
    colnames(MMPvalue) <- modNames
    modules[["MM"]] <- as.data.frame(geneModuleMembership)
    modules[["MMP"]] <- as.data.frame(MMPvalue)

    # calculate gene correlation to HIV expression, along with its p-value
    geneTraitSignificance <- cor(myExprMatrix, hiv, use = "p")
    GSPvalue              <- corPvalueStudent(as.matrix(geneTraitSignificance), nSamples)
    GS <- data.frame(TS = geneTraitSignificance, TSP = GSPvalue)
    colnames(GS) <- c("GS", "GSP")
    modules[["GS"]] <- GS
    
    # return
    modules
}
