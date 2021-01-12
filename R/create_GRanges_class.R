#' @export
create_GRanges_class <- function(methy, cpg_ids){
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
	cpg_ids <- unlist(strsplit(cpg_ids, ","))
	fd <- fd[cpg_ids,]
	betas <- jitter(betas[cpg_ids,])
	fd  <- fd_cols(fd)
	rownames(fd) <- rownames(betas)
	gr <- GenomicRanges::makeGRangesFromDataFrame(fd)
	S4Vectors::values(gr) <- betas
	return(gr)
}
