process.sc_hiv <- function  () {
    # load libraries needed for the analysis
    suppressMessages(library(scran))
    suppressMessages(library(scater))
    suppressMessages(library(SingleCellExperiment))

    # file names
    matrices_dir <- "/home/rcortini/work/CRG/projects/sc_hiv/data/matrices"

    # init variables
    sample_names <- c("P2449", "P2458")
    normalized <- list()
    dec <- list()

    # gene annotations file
    gene_annotations <- sprintf("%s/gene_annotations.tsv", matrices_dir)
    gene_data <- read.delim(gene_annotations, header = TRUE, row.names = 1, sep = "\t")

    # load sample sheets and expression matrices
    for (sample_name in sample_names) {
	
	# build file names
	matrix_fname <- sprintf("%s/%s.tsv.gz", matrices_dir, sample_name)
	sample_sheet_fname <- sprintf("%s/monocle/%s.pd.tsv", matrices_dir, sample_name)
	
	# parse files
	expr_matrix  <- read.table(matrix_fname, header = TRUE, row.names = 1,
				   sep = "\t", check.names = FALSE)
	sample_sheet <- read.delim(sample_sheet_fname,
				   header = TRUE, row.names = 1)
	
	# build the SingleCellExperiment object
	sce <- SingleCellExperiment(list(counts=as.matrix(expr_matrix)),
				    rowData = DataFrame(gene_data),
				    colData = DataFrame(sample_sheet))
	
	# do quality control
	sce <- calculateQCMetrics(sce, compact=TRUE)
	QC <- sce$scater_qc
	low.lib <- isOutlier(QC$all$log10_total_counts,
			     type="lower", nmad=3)
	low.genes <- isOutlier(QC$all$log10_total_features_by_counts,
			       type="lower", nmad=3)
	discard <- low.lib | low.genes
	sce <- sce[,!discard]
	
	# estimate size factors
	sizes <- c(sum(sce$label == "Jurkat"),
		   sum(sce$label == "J-Lat+DMSO"),
		   sum(sce$label == "J-Lat+SAHA"))
	sce <- computeSumFactors(sce, sizes=sizes)
	
	# normalize
	normalized[[sample_name]] <- normalize(sce)
	
	# now normalize based on the sample labels to obtain a list of genes
	# that have maximal biological variation
	collected <- list()
	for (x in unique(sce$label)) {
	    current <- sce[, sce$label==x]
	    current <- normalize(current)
	    fit <- trendVar(current, parametric=TRUE, use.spikes=FALSE) 
	    dec[[sample_name]] <- decomposeVar(current, fit)
	    collected[[x]] <- dec[[sample_name]]
	}
	
	# compile list of genes
	dec[[sample_name]] <- do.call(combineVar, collected)
	dec[[sample_name]]$gene_symbol <- rowData(sce)$gene_symbol
	ordering <- order(dec[[sample_name]]$bio, decreasing=TRUE)
	dec[[sample_name]] <- dec[[sample_name]][ordering,]
    }

    # determine which are the top-varying genes and intersect the lists
    top <- lapply(dec, function(x) {y <- rownames(x)[seq_len(1000)]})
    chosen <- Reduce(intersect, top)

    # create a list with the log-counts of the chosen genes
    original <- lapply(normalized, logcounts)

    # and invoke the MNN correction algorithm to remove batch effects
    corrected <- do.call(mnnCorrect, c(original, list(k=20, sigma=0.1)))

    # define a new SingleCellExperiment object, with the values of the corrected
    # expressions
    omat <- do.call(cbind, original)
    mat <- do.call(cbind, corrected$corrected)
    colnames(mat) <- NULL
    sc_hiv <- SingleCellExperiment(list(counts=omat, corrected=mat))
    colData(sc_hiv)$Batch <- rep(sample_names, lapply(corrected$corrected, ncol))
    colData(sc_hiv)$label <- unlist(lapply(normalized,
				     function(x) {as.character(x$label)}),
				     use.names = FALSE)
    # load the list of pairs of genes
    hs.pairs <- readRDS(system.file("exdata", "human_cycle_markers.rds", package="scran"))

    # prepare the data for cyclone
    gene_short_names <- gsub("\\..*", "", rownames(gene_data))
    non_duplicated <- !duplicated(gene_short_names)
    mygenes <- rownames(gene_data)[non_duplicated]
    sce_temp = sc_hiv[mygenes, ]
    rownames(sce_temp) <- gene_short_names[non_duplicated]

    # do the assignment of the cell cycle phases
    assignments <- cyclone(sce_temp , hs.pairs)
    sc_hiv$phases <- assignments$phases
    sc_hiv$scores <- assignments$scores
    sc_hiv
}

# load the WGCNA library and allow multithreading
suppressMessages(library(WGCNA))
suppressMessages(allowWGCNAThreads())

# This function prepares the data structure to be fed to the next function
PrepareDataForClustering <- function (exprMatrix, sampleSheet,
                          ngenes = 3600,
                          cut = 1000) {
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
    
    # plot size
    options(repr.plot.width = 10, repr.plot.height = 6)

    # detect outliers
    sampleTree <- hclust(dist(datExpr0), method = "average");
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
    
    # Determine cluster under the line
    clust <- cutreeStatic(sampleTree, cutHeight = cut, minSize = 10)
    table(clust)
    
    # clust 1 contains the samples we want to keep.
    keepSamples <- (clust == 1)
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
