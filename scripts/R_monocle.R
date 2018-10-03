library(monocle)
library(reshape2)

# file names
matrices_dir <- "/home/rcortini/work/CRG/projects/sc_hiv/data/matrices"
sample_name <- "P2449"
sample_sheet_fname <- sprintf("%s/monocle/%s.pd.tsv", matrices_dir, sample_name)
gene_annotation_fname <- sprintf("%s/monocle/%s.fd.tsv", matrices_dir, sample_name)
expr_matrix_fname <- sprintf("%s/%s.tsv.gz", matrices_dir, sample_name)

# load data
sample_sheet <- read.delim(sample_sheet_fname, header = TRUE, row.names = 1)
gene_annotation <- read.delim(gene_annotation_fname, header = TRUE, row.names = 1)
expr_matrix <- read.table(expr_matrix_fname, header = TRUE, row.names = 1,
			  sep = "\t", check.names = FALSE)

# apply constructs from monocle package
pd <- new("AnnotatedDataFrame", data = sample_sheet)
fd <- new("AnnotatedDataFrame", data = gene_annotation)
HSMM <- newCellDataSet(as.matrix(expr_matrix),
    phenoData = pd, featureData = fd, expressionFamily = negbinomial.size())

# estimate size factors and dispersions
HSMM <- estimateSizeFactors(HSMM)
HSMM <- estimateDispersions(HSMM)

# add total expression to experiment phenoData
pData(HSMM)$Total_mRNAs <- Matrix::colSums(exprs(HSMM))

# this step now is to eliminate the cells from the data set that have too few or
# too many reads
upper_bound <- 10^(mean(log10(pData(HSMM)$Total_mRNAs)) +
            2*sd(log10(pData(HSMM)$Total_mRNAs)))
lower_bound <- 10^(mean(log10(pData(HSMM)$Total_mRNAs)) -
            2*sd(log10(pData(HSMM)$Total_mRNAs)))
HSMM <- HSMM[,pData(HSMM)$Total_mRNAs > lower_bound &
      pData(HSMM)$Total_mRNAs < upper_bound]
HSMM <- detectGenes(HSMM, min_expr = 0.1)

# this generates the list of genes that are expressed in at least 10 cells
expressed_genes <- row.names(subset(fData(HSMM),
    num_cells_expressed >= 10))

