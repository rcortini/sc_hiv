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
