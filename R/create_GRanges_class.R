#' @export
create_GRanges_class <- function(methy, epi_res, sam, chr){
	if(class(methy) == "GenomicRatioSet") {
		betas <- minfi::getBeta(methy)
		pd <- as.data.frame(SummarizedExperiment::colData(methy))
		fd <- as.data.frame(SummarizedExperiment::rowRanges(methy))
		rownames(fd) <- rownames(betas)
	} else if(class(methy) == "ExpressionSet") {
		betas <- Biobase::exprs(methy)
		pd <- Biobase::pData(methy)
		fd <- Biobase::fData(methy)
	} else {
		stop("Input data 'methy' must be a 'GenomicRatioSet' or an 'ExpressionSet'")
	}
	cpg_ids <- epi_res[,"cpg_ids"]
	cpg_ids <- unlist(strsplit(cpg_ids, ","))
	fd <- fd[cpg_ids,]
	betas <- jitter(betas[cpg_ids,])
	ranges <- .Call("C_solve_user_SEW0", as.integer(fd$start), end = as.integer(fd$end), width = as.integer(fd$width), 
					PACKAGE = "IRanges")

	gr <- GenomicRanges::GRanges(seqnames= fd$seqnames, ranges = ranges, strand = fd$strand, metadata = betas)
	names(GenomicRanges::mcols(gr)) <- colnames(betas)
	return(gr)
}