library(RRHO2)
library('dplyr')

listA <- read.table("region1_foldchanges.limma.common.csv", sep=',', header = T)
listA <- na.omit(listA)
listA_input <-  data.frame(X=listA$probe, padj=listA$log2ratio)
listB<-read.table("region2_foldchanges.limma.common.csv",sep=',', header = T)
listB<-na.omit(listB)
listB_input<- data.frame(X=listB$probe, padj=listB$log2ratio)
alt='enrichment'
xlab="region1"
ylab="region2"

#Subsetting probes
listA_input_small<-sample_n(listA_input,50000)
probes_to_keep<-listA_input_small$X
listA_input_small<-listA_input_small[listA_input_small$X %in% probes_to_keep,]
listB_input_small<-listB_input[listB_input$X %in% probes_to_keep,]

RRHO_obj <-  RRHO2_initialize(listA_input_small, listB_input_small, labels = c(xlab,ylab), log10.ind=TRUE)
jpeg(filename = "RRHO2_heatmap.limma.jpg")
# RRHO2 Heat map
hypermat <- RRHO_obj$hypermat
labels <- RRHO_obj$labels     
maximum <- max(hypermat,na.rm=TRUE)
minimum <- min(hypermat,na.rm=TRUE)
color.bar <- function(lut, min, max=-min, 
                        nticks=11, 
                        ticks=seq(min, max, len=nticks), 
                        title='RRHO2 heatmap region1 vs region2') {
    scale  <- (length(lut)-1)/(max-min)
    plot(c(0,10), c(min,max), type='n', bty='n', 
         xaxt='n', xlab='', yaxt='n', ylab='')
    mtext(title,2,2.3, cex=0.8)
    axis(2, round(ticks,0), las=1,cex.lab=0.8)
    for (i in 1:(length(lut)-1)) {
      y  <- (i-1)/scale + min
      rect(0,y,10,y+1/scale, col=lut[i], border=NA)}
  }
jet.colors  <- colorRampPalette(
            c("#00007F", "blue", "#007FFF", "cyan", 
              "#7FFF7F", "yellow", "#FF7F00", "red", "#7F0000"))
        
  colorGradient <- jet.colors(101)
  layout(matrix(c(rep(1, 6), 2), 1, 7, byrow = TRUE))
  breaks <- seq(minimum,maximum,length.out = length(colorGradient) + 1)
  image(hypermat, col = colorGradient,breaks=breaks,
        axes = FALSE, useRaster = TRUE)
  mtext(labels[2],2,0.5)
  mtext(labels[1],1,0.5)
atitle <- ifelse(RRHO_obj$log10.ind, "-log10(P-value)", "-log(P-value)")
  color.bar(colorGradient, min = minimum, max = maximum, nticks = 6, title = atitle)
dev.off()

jpeg(filename = "RRHO2_venn_diagram_uu.limma.jpg")
RRHO2_vennDiagram(RRHO_obj, type="uu")
dev.off()
