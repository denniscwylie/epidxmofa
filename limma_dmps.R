library("WGCNA")
library(data.table)
library(methylCIPHER)
library('limma')
library('ggplot2')
library('ggrepel')
library('dplyr')
library('qqman')
library('tidyr')

beta<-fread("processed_data_debatched_tsv", sep="\t",header=TRUE)
#making the data table a dataframe
beta<-data.frame(beta)
rownames(beta)<-beta$id
beta$id<-NULL
samples<-gsub("X","",colnames(beta))
colnames(beta)<-samples

#read in metadata and grab the columns of interest
metadata<-read.csv("metadata_selected_withcellestimates_allpreds_allregions.csv", sep=",", header=TRUE)
metadata_small<-metadata[,c("Sample_Name","BrNum","Sex.x","PrimaryDx","AgeDeath","Race","Smoking","Manner_Of_Death")]
rownames(metadata_small)<-metadata_small$Sample_Name
which(!(colnames(beta) %in% rownames(metadata_small)))

samplestokeep<-colnames(beta[,-c(66,116)])
beta<-beta[,samplestokeep]
metadata_small<-metadata_small[samplestokeep,]

#probe annotation
annotation<-read.csv("probe_annotation.csv", header=TRUE, row.names=1)
annotation_small<-annotation[,c("UCSC_RefGene_Name","UCSC_RefGene_Name_Simple")]


#MDD vs control
MDD_control_samples<-rownames(metadata_small[metadata_small$PrimaryDx!='PTSD',])
MDD_control_beta<-beta[,MDD_control_samples]
MDD_control_metadata_small<-metadata_small[MDD_control_samples,]
design <- model.matrix(~ PrimaryDx+Sex.x+Race+Smoking+Manner_Of_Death, data=MDD_control_metadata_small)
fit <- lmFit(MDD_control_beta, design)
fit2 <- eBayes(fit)
dat<-topTable(fit2, coef = "PrimaryDxMDD", number = Inf, adjust.method = "fdr")

#how many dmps meet adjpvalue cutoff of 0.05
table(dat$adj.P.Val<=0.05)

#how many dmps meet pvlue cutoff of 0.00005
table(dat$p<=0.00005)

#merge annotation to the dat dataframe
annotation_small<-annotation[,c('chr','pos','strand','UCSC_RefGene_Name','UCSC_RefGene_Name_Simple')]
annotation_small$cpg<-rownames(annotation_small)
dat<-as.data.table(dat)
setkey(dat,cpg)
annotation_small<-as.data.table(annotation_small)
setkey(annotation_small,cpg)
merged_dat<-annotation_small[dat, on="cpg"]

#write out all results
write.csv(merged_dat, "mdd_dmps.all.csv", quote=FALSE)

png("mdd_vs_control_manhattanplot.png", width = 12, height = 6, units = "in",res = 300)
manhattan(merged_dat, chr="chr",bp="pos",p="p",snp="UCSC_RefGene_Name_Simple",annotatePval =1e-5, suggestiveline=-log10(1e-5))
dev.off()

#PTSD vs control
PTSD_control_samples<-rownames(metadata_small[metadata_small$PrimaryDx!='MDD',])
PTSD_control_beta<-beta[,PTSD_control_samples]
PTSD_control_metadata_small<-metadata_small[PTSD_control_samples,]
design <- model.matrix(~ PrimaryDx+Sex.x+Race+Smoking+Manner_Of_Death, data=PTSD_control_metadata_small)
fit <- lmFit(PTSD_control_beta, design)
fit2 <- eBayes(fit)
dat<-topTable(fit2, coef = "PrimaryDxPTSD", number = Inf, adjust.method = "fdr")

#how many dmps meet adjpvalue cutoff of 0.05
table(dat$adj.P.Val<=0.05)

#how many dmps meet pvlue cutoff of 0.00005
table(dat$p<=0.00005)

#merge annotation to the dat dataframe
annotation_small<-annotation[,c('chr','pos','strand','UCSC_RefGene_Name','UCSC_RefGene_Name_Simple')]
annotation_small$cpg<-rownames(annotation_small)
dat<-as.data.table(dat)
setkey(dat,cpg)
annotation_small<-as.data.table(annotation_small)
setkey(annotation_small,cpg)
merged_dat<-annotation_small[dat, on="cpg"]

#write out all results
write.csv(merged_dat, "ptsd_dmps.all.csv", quote=FALSE)

png("ptsd_vs_control_manhattanplot.png", width = 12, height = 6, units = "in",res = 300)
manhattan(merged_dat, chr="chr",bp="pos",p="p",snp="UCSC_RefGene_Name_Simple",annotatePval =1e-5, suggestiveline=-log10(1e-5))
dev.off()
