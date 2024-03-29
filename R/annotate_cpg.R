#' @title Annotate the DMR resulting from epimutacions package
#' 
#' @description This function annotates a differentially 
#' methylated region
#' @param data DataFrame-like object.
#' @param db a character string specifying the 
#' Database to use for annotation. 
#' E.g: \code{'IlluminaHumanMethylationEPICanno.ilm10b2.hg19'}.
#' @param split a character string containing the separator for CpG ids. 
#' Default \code{','}.
#' @param epi_col CpG ids, should be row names in the data base. 
#' @param gene_col column name from where to extract gene names. 
#' Default: \code{'GencodeBasicV12_NAME'}.
#' @param feat_col column name from where to extract CpG feature groups. 
#' Default: \code{'Regulatory_Feature_Group'}.
#' @param relat_col column name from where 
#' to extract relation to island info. 
#' Default: \code{'Relation_to_Island'}.
#' @param build The build for bioMart. Default \code{'37'}.
#' @param omim a boolean, if TRUE will annotate OMIMs as well. 
#' Takes a bit longer. Default TRUE.
#' @return The function returns a 
#' DataFrame-like object annotated.
#' 
#' @importFrom minfi getAnnotation
#' @importFrom biomaRt useEnsembl getBM

annotate_cpg <- function(data, db, split = ',', 
        # illumina annotation parameters
        epi_col = 'cpg_ids', gene_col = 'GencodeBasicV12_NAME', 
        feat_col= 'Regulatory_Feature_Group', relat_col = 'Relation_to_Island',
        # biomart parameters
        build = "37", omim = TRUE)
{
    
    # get illumina data base
    message('Loading DB')
    anno <- minfi::getAnnotation(db)
    epids_list <- data[[epi_col]]
    epids_list <- strsplit(x = epids_list, split = split)
    
    # Basic annotation
    message('Annotating:')
    
    message(gene_col)
    annotated_genes <- lapply(epids_list, function(x) {
            x_in_anno <- which(x %in% rownames(anno))
            x <- x[x_in_anno]
            not_empty <- which(anno[x, ][[gene_col]] != "")
            x <- x[not_empty]
            if (length(x) != 0) {
                paste0(anno[x, ][[gene_col]], collapse = ' /// ')
            }
        })
    
    message(feat_col)
    annotated_regions <- lapply(epids_list, function(x) {
            paste(anno[x, ][[feat_col]], collapse = '///')
        })
    
    message(relat_col)
    annotated_relation <- lapply(epids_list, function(x) {
            paste(anno[x, ][[relat_col]], collapse = '///')
        })
    
    data[[gene_col]] <- annotated_genes
    data[[feat_col]] <- remove_slash(annotated_regions)
    data[[relat_col]] <- remove_slash(annotated_relation)
    data <- data[-which(vapply(data[[gene_col]], is.null, logical(1))), ]

    # OMIM annotation
    if (omim == TRUE) {
        # Using biomart  to extract OMIMs
        message('- Extracting and annotating OMIMs')

        gene_mart <- biomaRt::useEnsembl( biomart = "ENSEMBL_MART_ENSEMBL",
                                        GRCh = build,
                                        host = "https://www.ensembl.org",
                                        #..# host = "www.ensembl.org",
                                        dataset = "hsapiens_gene_ensembl" )
        
        annotated_omim <- lapply(data[[gene_col]], function(gene_list) {
            gene_list <- unlist(strsplit(x = gene_list, split = '///'))
            gene_list <- unlist(strsplit(x = gene_list, split = ';'))
            
            # extract OMIM info from db
            omims <- biomaRt::getBM( mart = gene_mart,
                                    attributes = c(
                                        "hgnc_symbol",
                                        "mim_morbid_accession",
                                        "mim_morbid_description"
                                    ),
                                    filters = "hgnc_symbol",
                                    values = gene_list )
            
            # get accession per gene
            accesions <- lapply(gene_list, function(gene) {
                acc <- omims[omims$hgnc_symbol == gene, ]$mim_morbid_accession
                acc[is.na(acc)] <- ''
                paste(acc, collapse = '/')
            })
            # get description per gene
            descriptions <- lapply(gene_list, function(gene) {
                descr <- omims[omims$hgnc_symbol==gene, ]$mim_morbid_description
                descr[is.na(descr)] <- ''
                paste(descr, collapse = '/')
            })
            # per gene set
            omims_acc <- paste(accesions, collapse = '///')
            omims_acc <- remove_slash(omims_acc)
            omims_desc <- paste(descriptions, collapse = '///')
            omims_desc <- remove_slash(omims_desc)
            return(list('OMIM_ACC' = omims_acc, 'OMIM_DESC' = omims_desc))
        })
        
        # Appending OMIM data
        annotated_omim <- do.call('rbind', annotated_omim)
        data <- cbind(data, annotated_omim)
    }
    return(data)
}

remove_slash <- function( sequence )
{
    sequence <- gsub('(/)\\1+', '\\1', sequence)
    sequence <- gsub('(/$)', '', sequence)
    sequence <- gsub('(^/)', '', sequence)
    
    return(sequence)
    
}