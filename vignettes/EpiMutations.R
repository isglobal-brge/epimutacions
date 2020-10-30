## ----setup, include = FALSE---------------------------------------------------
knitr::opts_chunk$set(comment="", warning=FALSE, message=FALSE, cache=TRUE)


## -----------------------------------------------------------------------------
library(EpiMutations)
#Control samples
##GenomicRatioSet
data(grs.control.panel)
##ExpressionSet
data(es.control.panel)

#Case samples
##GenomicRatioSet
data(grs.diseases)
##ExpressionSet
data(es.diseases)

## -----------------------------------------------------------------------------
disease<-es.diseases[,4]
dim(disease)[2]

## -----------------------------------------------------------------------------
#manova (default method)

EpiMutations(diseases = disease)

## -----------------------------------------------------------------------------
#Multivariate inear model

EpiMutations(diseases = disease, method = "mlm")

## -----------------------------------------------------------------------------
#Isolation forest

EpiMutations(diseases = disease, method = "iso.forest")

## -----------------------------------------------------------------------------
#Robust mahalanobis distance

EpiMutations(diseases = disease, method = "Mahdist.MCD")

