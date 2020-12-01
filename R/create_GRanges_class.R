#' @export
create_GRanges_class <- function(methy, epi_res, sam, chr){
	cpg_ids <- epi_res[epi_res$chromosome == chr,"cpg_ids"]
	cpg_ids <- unlist(strsplit(cpg_ids, ","))
	fd <- as.data.frame(SummarizedExperiment::rowRanges(methy))
	fd <- fd[cpg_ids,]
	betas <- minfi::getBeta(methy)
	betas <- jitter(betas[cpg_ids,])
	gr <-GenomicRanges::GRanges(seqnames= fd$seqnames, ranges = IRanges::IRanges(start = fd$start, end = fd$end), strand = fd$strand, metadata = betas)
	names(GenomicRanges::mcols(gr)) <- colnames(betas)
	return(gr)
}
