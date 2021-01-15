#' @export
create_GRanges_class <- function(methy, cpg_ids){
	if(class(methy) == "GenomicRatioSet") {
		betas <- minfi::getBeta(methy)
		pd <- as.data.frame(SummarizedExperiment::colData(methy))
		fd <- as.data.frame(SummarizedExperiment::rowRanges(methy))
		rownames(fd) <- rownames(betas)
	} else {
		stop("Input data 'methy' must be a 'GenomicRatioSet'")
	}
	
	cpg_ids <- unlist(strsplit(cpg_ids, ","))
	fd <- fd[cpg_ids,]
	fd  <- .fd_cols(fd)
	betas <- jitter(betas[cpg_ids,])
	rownames(fd) <- rownames(betas)
	gr <- GenomicRanges::makeGRangesFromDataFrame(fd)
	S4Vectors::values(gr) <- betas
	return(gr)
}
