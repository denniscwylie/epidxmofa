#!/usr/bin/env Rscript

library(dplyr)
library(ggplot2); theme_set(theme_bw())
library(MOFA2)
library(pheatmap)
library(scales)
library(tidyr)

source("load_mofa_from_hdf5.R")

unsquash_trans = function() {trans_new(
    "unsquash",
    function(x) {ifelse(abs(x) < 1, x, sign(x) * sqrt(abs(x)))},
    function(u) {ifelse(abs(u) < 1, u, sign(u) * u * u)}
)}

cellTyped = setdiff(rownames(facs), paste0("Br270", 0:2))
cellTypeColumns = grep("enriched", colnames(demo), value=TRUE)

pdf("plots/mofa_factor_vs_celltype_spearman_heatmap.pdf", h=9, w=5)
pheatmap(cor(facs[cellTyped, ],
             demo[cellTyped, sort(cellTypeColumns)],
             method = "spearman"),
         colorRampPalette(c("darkblue", "skyblue", "white", "goldenrod", "red"))(101),
         cluster_rows = FALSE,
         cluster_cols = FALSE)
garbage = dev.off()

for (theCellType in cellTypeColumns) {
    ggd = data.frame(sample = cellTyped,
                     scale(facs)[cellTyped, ],
                     demo[cellTyped, c("PrimaryDx", "AgeDeath")],
                     check.names = FALSE) %>%
        pivot_longer(c(-sample, -PrimaryDx, -AgeDeath),
                     names_to = "Factor",
                     values_to = "Score_Scaled") %>%
        as.data.frame
    ggd[ , theCellType] = demo[ggd$sample, theCellType]
    thePs = apply(facs[cellTyped, ], 2, function(f) {
        cor.test(f, demo[cellTyped, theCellType],
                 method="spearman")$p.value
    })
    ggd$Factor = factor(as.character(ggd$Factor),
                        levels=paste0("Factor", 1:30))
    levels(ggd$Factor) = paste0("Factor ", 1:30,
                                "\n(p = ", signif(thePs, 2), ")")
    levels(ggd$Factor) = gsub("p = 0\\)", "p ~ 0\\)", levels(ggd$Factor))
    gg = ggplot(ggd, aes(x = .data[[theCellType]],
                         y = Score_Scaled,
                         color = PrimaryDx))
    gg = gg + facet_wrap(~ Factor, scales="free_y")
    gg = gg + geom_point(size=0.3, alpha=0.5)
    gg = gg + stat_smooth(method.args=list(degree=1, span=1.5))
    gg = gg + scale_color_manual(values=c("slategray", "dodgerblue", "orangered"))
    gg = gg + ylab("Score (Scaled)")
    gg = gg + scale_y_continuous(trans="unsquash")
    gg = gg + theme(panel.grid.minor = element_blank())
    pdf(paste0("plots/mofa_factors_vs_",
               tolower(theCellType), ".pdf"),
        h=9, w=13)
    print(gg)
    garbage = dev.off()
}
