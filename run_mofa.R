#!/usr/bin/env Rscript

library(data.table)
library(dplyr)
library(ggplot2); theme_set(theme_bw())
library(limma)  ## for removeBatchEffect
library(matrixStats)
library(MASS)
library(MOFA2)
library(pheatmap)
library(scales)
library(tidyr)
library(WriteXLS)

demo = readRDS("demo.rds")
brNums = readRDS("brNums.rds")
methylData = readRDS("methylData.rds")


## -----------------------------------------------------------------------------
## convert beta values to M-values:
methylM = lapply(methylData, function(meth) {
    meth = as.matrix(meth)
    meth[meth < (1/256)] = 1 / 256
    return(log2(meth / (1-meth)))
})
methylSdCutoffs = lapply(methylM, function(.) {quantile(rowSds(.), 0.90)})
methylForMofas = mapply(
    function(mat, cutoff) {sdm=rowSds(mat); mat[sdm >= cutoff, ]},
    methylM,
    methylSdCutoffs
)

## -----------------------------------------------------------------------------
## relabel samples by brain number and align for mofa:
methylForMofas = mapply(function(meth, brNum) {
    colnames(meth) = gsub("^X", "", colnames(meth))
    colnames(meth) = brNum[colnames(meth)]
    return(meth)
}, methylForMofas, brNums, SIMPLIFY=FALSE)
commonBr = Reduce(intersect, lapply(methylForMofas, colnames))
methylForMofas = lapply(methylForMofas, function(meth) {meth[ , commonBr]})

## -----------------------------------------------------------------------------
methDebatched = lapply(methylForMofas, function(meth) {removeBatchEffect(
    meth,
    batch = demo[commonBr, "year"],
    design = model.matrix(~ Sex + AgeDeath + Race + PrimaryDx,
                          data = demo[commonBr, ])
)})
## replace methylForMofas with debatched version:
methylForMofas = methDebatched

## -----------------------------------------------------------------------------
set.seed(123)
## reticulate::use_python("/home/dennis/bin/python3", required=TRUE)
mofaObj = create_mofa(methylForMofas)
##
modOpt = get_default_model_options(mofaObj)
modOpt$num_factors = 30
##
mofaObj = prepare_mofa(object = mofaObj,
                       data_options = get_default_data_options(mofaObj),
                       model_options = modOpt,
                       training_options = get_default_training_options(mofaObj))
mofaObj = run_mofa(mofaObj, "epidxmofa.hdf5"## ,
                   ## use_basilisk = FALSE
                   )

## -----------------------------------------------------------------------------
n = sum(mofaObj@dimensions$N)
meta = data.frame(
    row.names = samples_names(mofaObj)[[1]],
    sample = samples_names(mofaObj)[[1]],
    subj_id = gsub("_.*", "", samples_names(mofaObj)[[1]])
)
meta$age = demo[meta$subj_id, "AgeDeath"]
meta$race  = demo[meta$subj_id, "Race"]
meta$sex = demo[meta$subj_id, "Sex"]
meta$dx = demo[meta$subj_id, "PrimaryDx"]
meta$year = demo[meta$subj_id, "year"]
meta$trauma_class = factor(demo[meta$subj_id, "trauma_class"])
enrichCols = grep("_enriched_", colnames(demo), value=TRUE)
meta[ , enrichCols] = demo[meta$subj_id, enrichCols]
##
samples_metadata(mofaObj) = meta
##
facs = get_factors(mofaObj)$group1
