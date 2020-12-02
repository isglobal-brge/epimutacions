#' @export
create_GRanges_class <- function(methy, epi_res, sam, chr){
	
	cpg_ids <- epi_res[,"cpg_ids"]
	cpg_ids <- unlist(strsplit(cpg_ids, ","))
	fd <- as.data.frame(SummarizedExperiment::rowRanges(methy))
	fd <- fd[cpg_ids,]
	betas <- minfi::getBeta(methy)
	betas <- jitter(betas[cpg_ids,])
	ranges <- .Call("C_solve_user_SEW0", as.integer(fd$start), end = as.integer(fd$end), width = as.integer(fd$width), 
					PACKAGE = "IRanges")

	gr <- GenomicRanges::GRanges(seqnames= fd$seqnames, ranges = ranges, strand = fd$strand, metadata = betas)
	names(GenomicRanges::mcols(gr)) <- colnames(betas)
	return(gr)
}