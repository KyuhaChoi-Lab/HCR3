.#' loading libraries ,input, goi
# loading libraries{{{
library(tidyverse)
library(colorspace)
library(extrafont)

loadfonts(device="pdf")
#}}}
###### NOTE #####
# RNAseq data info
#
# Caution! RNA-seq data obtained from different date was merged for DESeq analysis. Batch effect may concern.
#################

# inputs{{{
dirdeseq <- "/home/namilhand/01_Projects/hcr3/HCR3_public_release/RNA-seq/meiotic_gene_expression/deseq2/deseq_results"
wd <- "/home/namilhand/01_Projects/hcr3/HCR3_public_release/RNA-seq/meiotic_gene_expression"

# sel.samples <- c("col_s_1", "col_s_2", "col_s_3", "col_s_4", "col0_b_1", "col0_b_2", "col0_b_3", "hcr3_b_1", "hcr3_b_2", "hcr3_b_3")
## select samples from normalized data
# deseq <- dplyr::select(deseq, c("AGI", sel.samples))
# rlog <- dplyr::select(rlog, c("AGI", sel.samples))

## goi
goi <- read_csv("geneset/meiotic_genes_and_hsf.csv") %>%
		mutate(symbol=factor(symbol, levels=rev(symbol)))
colnames(goi)[1] <- "func"

## log2fold changes, DEG= padj < 0.01
# l2fc.col_b_s <- read_csv(file.path(dirdeseq, "deseq_results", "deseqResShrunken_col-b_vs_col-s.csv"), col_names=T)
l2fc.col_b_s <- read_csv(file.path(dirdeseq, "deseqResShrunken_col-b_vs_col-s.csv"), col_names=T)
l2fc.hcr3_col <- read_csv(file.path(dirdeseq, "deseqResShrunken_hcr3_vs_col_bud.csv"), col_names=T)
l2fc.j2_col <- read_csv(file.path(dirdeseq, "deseqResShrunken_j2_vs_col_bud.csv"), col_names=T)
l2fc.j3_col <- read_csv(file.path(dirdeseq, "deseqResShrunken_j3_vs_col_bud.csv"), col_names=T)

colnames(l2fc.col_b_s)[1] <- "AGI"
colnames(l2fc.hcr3_col)[1] <- "AGI"
colnames(l2fc.j2_col)[1] <- "AGI"
colnames(l2fc.j3_col)[1] <- "AGI"

#}}}

#' processing data for analysis
## join log2foldchange data to goi

join_l2fc <- function(query, l2fc){
    # query <- goi
    # l2fc <- l2fc.hcr3_col_b
    result <- left_join(goi, l2fc, by="AGI") %>%
        mutate(log2FoldChange = case_when(padj >= 0.01 ~ 0,
                                        TRUE ~ log2FoldChange)) %>%
        dplyr::select(c("func", "symbol", "AGI", "log2FoldChange"))
}

goi.l2fc.col_b_s <- join_l2fc(goi, l2fc.col_b_s)
goi.l2fc.hcr3_col <- join_l2fc(goi, l2fc.hcr3_col)
goi.l2fc.j2_col <- join_l2fc(goi, l2fc.j2_col)
goi.l2fc.j3_col <- join_l2fc(goi, l2fc.j3_col)

mei.l2fc <- left_join(goi.l2fc.col_b_s, goi.l2fc.hcr3_col, by=c("func", "symbol", "AGI")) %>%
    left_join(goi.l2fc.j2_col, by=c("func", "symbol", "AGI")) %>%
    left_join(goi.l2fc.j3_col, by=c("func", "symbol", "AGI"))

colnames(mei.l2fc)[4:7] <- c("bud/sdl", "hcr3/col", "j2/col", "j3/col")

mei.l2fc <- pivot_longer(mei.l2fc, 4:7, names_to="group", values_to="log2fc") %>%
		filter(func  != "HSF") %>%
        mutate(log2fc = case_when(abs(log2fc) < 0.5 ~ 0,
                                   TRUE ~ log2fc ))
##}}}


#' Draw heatmap
## function: fcHeat{{{
fcHeat <- function(dat){
		p <- ggplot(dat, aes(x=group, y=symbol, fill=log2fc)) +
				geom_tile(colour=NA, size=0) +
				theme_classic() +
				theme(legend.position="top",
					  axis.text.y=element_text(face="italic"),
					  axis.line=element_blank(),
					  axis.text.x=element_text(hjust=1, vjust=1),
					  legend.margin=margin(2, 2, 2, 2),
					  legend.key.size=unit(0.1, "inches"),
					  legend.title=element_blank()
					  ) +
				labs(x=NULL, y=NULL) +
				scale_fill_continuous_diverging("Blue-Red 3", limits=c(-2, 2), oob=scales::squish)
		return(p)
}
#}}}

p.mei.l2fc <- fcHeat(mei.l2fc)

svg(file=file.path(wd, "meiotic-genes_log2fc.svg"), width=1.5)
print(p.mei.l2fc)
dev.off()

