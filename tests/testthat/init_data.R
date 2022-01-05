
# Data
## Create cases and controls dataset
data(GRset)
case <- GRset[,11]
control <- GRset[,1:10]

# Run the epimutations function 
result_manova <- epimutations(case, control, "manova")
result_mlm <- epimutations(case, control, "mlm")
result_isoforest <- epimutations(case, control, "isoforest")
result_mahadistmcd <- epimutations(case, control, "mahdistmcd")
result_quantile <- epimutations(case, control, "quantile")
result_beta <- epimutations(case, control, "beta")

# Set variables to help in the testing
## Colnames
col_names <- c("epi_id", "sample", 
               "chromosome", "start", 
               "end", "sz","cpg_n", 
               "cpg_ids", "outlier_score", 
               "outlier_direction",
               "pvalue", "adj_pvalue",
               "delta_beta",
               "epi_region_id", 
               "CRE", "CRE_type") 