library(CETYGO)
library("data.table")
beta<-data.frame(fread("processed_data_debatched.tsv", sep="\t", check.names=FALSE))
rownames(beta)<-beta$id
beta$id<-NULL

predProp<-projectCellTypeWithError(beta,modelBrainCoef[["IDOL"]][[1]])
predProp2<-projectCellTypeWithError(beta,modelBrainCoef[["IDOL"]][[2]])
predProp3<-projectCellTypeWithError(beta,modelBrainCoef[["IDOL"]][[3]])
predProp4<-projectCellTypeWithError(beta,modelBrainCoef[["IDOL"]][[4]])
predProp5<-projectCellTypeWithError(beta,modelBrainCoef[["IDOL"]][[5]])
predProp6<-projectCellTypeWithError(beta,modelBrainCoef[["IDOL"]][[6]])
predProp7<-projectCellTypeWithError(beta,modelBrainCoef[["IDOL"]][[7]])
predProp8<-projectCellTypeWithError(beta,modelBrainCoef[["IDOL"]][[8]])

all_pred<-cbind(predProp, predProp2, predProp4, predProp5, predProp6, predProp7,predProp8)
#using panels 2 and 8
#neuronal enriched, Oligodendrocyte enriched, Microglia enriched, Astrocyte enriched,Inhibatory (gabaergic) neuronal enriched,Excitatory (glutamatergic) neuronal enriched
all_pred_subset<-all_pred[,c(7,8,6,9,10,38,37,41)]

rownames(all_pred_subset)<-gsub("^X","", rownames(all_pred_subset))
write.csv(all_pred_subset, "cellestimates.csv",quote=FALSE)
save.image("cetygo.Rdata")

