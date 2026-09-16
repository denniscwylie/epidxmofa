#!/usr/bin/env Rscript

library(MOFA2)

demo = readRDS("demo.rds")
brNums = readRDS("brNums.rds")
methylData = readRDS("methylData.rds")

mofaObj = load_model("epidxmofa.hdf5")

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
for (clm in c("ex", "inhib", "Tcell")) {
    for (reg in c("Amy", "Hipp", "mPFC")) {
        clmReg = paste0(clm, "_", reg)
        meta[[clmReg]] = demo[meta$subj_id, clmReg]
    }
}
meta$trauma_class = factor(demo[meta$subj_id, "trauma_class"])

samples_metadata(mofaObj) = meta

facs = get_factors(mofaObj)$group1
