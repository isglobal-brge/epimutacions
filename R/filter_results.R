filter_results <- function(rst, method, pvalue_cutoff = NULL, outlier_score_cutoff = NULL){
  
  if ( method == "manova" | method == "mlm"){
    
    rst <- rst[which(rst$outlier_significance < pvalue_cutoff),]
  }
  
  if (method == "isoforest"){
    rst <- rst[which(rst$outlier_score > outlier_score_cutoff),]
  }
  return(rst)
}