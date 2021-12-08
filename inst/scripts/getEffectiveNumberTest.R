#' Estimate effective number of tests. 
#' 
#' Steps:
#' 1. Make a dataset with two samples. This dataset contains all the CpGs from the array
#' and one samples has all 0s and the other all 1s.
#' 2. Run bumphunter using the same configuration than in the package.
#' 3. Get the total number of regions computed. This is total 
#' number of effective tests for methods with p-value.

# Load libraries
library(minfi)

# Regions 450K
library(IlluminaHumanMethylation450kanno.ilmn12.hg19)

data(Locations)

### Select CpGs (names starting by cg) in autosomic chromosomes
locs.450 <- subset(Locations, 
                   grepl("^cg", rownames(Locations)) & 
                   chr %in% paste0("chr", 1:22))
locs.450GR <- makeGRangesFromDataFrame(locs.450, 
                                       start.field = "pos",
                                       end.field = "pos", 
                                       strand = "*")
locs.450GR <- sort(locs.450GR)
mat <- matrix(0, 
              nrow = length(locs.450GR), 
              ncol = 2, 
              dimnames = list(names(locs.450GR), c("A", "B")))

## Set sample B to all 1
mat[, 2] <- 1

## Define model matrix
pheno <- data.frame(var = c(0, 1))
model <- model.matrix(~ var, pheno)

## Run bumphunter
bumps <- bumphunter(mat, design = model, pos = start(locs.450GR), 
                    chr = as.character(seqnames(locs.450GR)),
                    cutoff = 0.05)$table
bumps.fil <- subset(bumps, L >= 3)
nrow(bumps.fil)
# 44244

## Epic
library(IlluminaHumanMethylationEPICanno.ilm10b2.hg19)
data(Locations)

### Select CpGs (names starting by cg) in autosomic chromosomes
locs.EPIC <- subset(Locations, 
                    grepl("^cg", rownames(Locations)) & 
                    chr %in% paste0("chr", 1:22))
locs.EPICGR <- makeGRangesFromDataFrame(locs.EPIC, 
                                        start.field = "pos", 
                                        end.field = "pos",
                                        strand = "*")
locs.EPICGR <- sort(locs.EPICGR)
matEpic <- matrix(0, nrow = length(locs.EPICGR), ncol = 2, 
              dimnames = list(names(locs.EPICGR), c("A", "B")))

## Set sample B to all 1
matEpic[, 2] <- 1

## Define model matrix
pheno <- data.frame(var = c(0, 1))
model <- model.matrix(~ var, pheno)

## Run bumphunter
bumpsEpic <- bumphunter(matEpic, design = model, pos = start(locs.EPICGR), 
                    chr = as.character(seqnames(locs.EPICGR)),
                    cutoff = 0.05)$table
bumpsEpic.fil <- subset(bumpsEpic, L >= 3)
nrow(bumpsEpic.fil)
# 66324
