library('ggplot2')
library(reshape2)
library("ggpubr")

cellCountEstimates<-read.table("metadata_withcellestimates_allpreds_allregions.csv", sep=",", header=TRUE, row.names = 1)
table(cellCountEstimates$brainregion)
amg<-cellCountEstimates[cellCountEstimates$brainregion=='Amygdala',]

amg_mod<-amg[,c("Neuronal_enriched","Oligodendrocyte_enriched","Microglia_enriched","Astrocyte_enriched","Inhibatory_gabaergic_neuronal_enriched","Excitatory_glutamatergic_neuronal_enriched","PrimaryDx")]
colnames(amg_mod)<-c("Neuronal","Oligo","Microglia","Astro","Inhib_gabaergic_neuronal","Excit_glutamatergic_neuronal","PrimaryDx")
melt<-melt(amg_mod, id.vars='PrimaryDx', measure.vars=c("Neuronal","Oligo","Microglia","Astro","Inhib_gabaergic_neuronal","Excit_glutamatergic_neuronal")) 
colnames(melt)<-c("PrimaryDx","CellType", "Cell_Type_Estimate")
melt$PrimaryDx<-factor(melt$PrimaryDx)

#pvalues
compare_means(Oligodendrocyte_enriched ~ PrimaryDx, data =amg)
compare_means(Neuronal_enriched ~ PrimaryDx, data=amg)
compare_means(Microglia_enriched ~ PrimaryDx, data=amg)
compare_means(Astrocyte_enriched ~ PrimaryDx, data=amg)
compare_means(Inhibatory_gabaergic_neuronal_enriched~ PrimaryDx, data=amg)
compare_means(Excitatory_glutamatergic_neuronal_enriched~ PrimaryDx, data=amg)

pdf('boxplots_brainregions.allpred.pdf')
ggboxplot(melt, x = "PrimaryDx", y = "Cell_Type_Estimate",fill="PrimaryDx", facet.by="CellType")+ggtitle("Amygdala")+scale_fill_manual(values=c("slategray","dodgerblue","orangered"))
dev.off()

