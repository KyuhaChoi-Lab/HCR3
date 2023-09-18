# ---------------------------------------------------------------------------------------------------
# template : featurecount2deseq2_pipeline.R
# date: 2021. 08. 13
# Description :
# Use this script as template for the DESeq analysis. Copy it to the place where DESeq2 output will be stored, then edit the file according to your specific use.
# Edit in the vim, and excute it in that buffer with Nvim-R plugin. 

#==== A. Setting Environment ========
suppressPackageStartupMessages({
  library("DESeq2")
  library("tidyverse")
  library("extrafont")
  library("pheatmap")
  library(org.At.tair.db)
})

font_import()
loadfonts(device="pdf")

#---- output variables --------
outputpath="/datasets/data_4/nison/hcr3/HCR3_public_release/RNA-seq/meiotic_gene_expression/deseq2"
dir.create(outputpath, recursive=TRUE)

#---- read cts file, a output of STAR --------
ctsfile1 <- "/datasets/data_4/nison/hcr3/HCR3_public_release/RNA-seq/meiotic_gene_expression/featureCount/hcr3-paper_total-RNAseq.featureCounts.txt"
ctsfile2 <- "/datasets/data_4/nison/hcr3/HCR3_public_release/RNA-seq/meiotic_gene_expression/featureCount_hei10/hcr3-paper_total-RNAseq_hei10.featureCounts.txt"
ctsfile3 <- "/datasets/data_4/nison/hcr3/HCR3_public_release/RNA-seq/meiotic_gene_expression/featureCount_hcr2/hcr2.featureCounts.txt"
ctsfile4 <- "/datasets/data_4/nison/hcr3/HCR3_public_release/RNA-seq/meiotic_gene_expression/featureCount_hcr2_hei10/hcr2-paper_hei10.featureCounts.txt"
ctsfile5 <- "/datasets/data_4/nison/hcr3/HCR3_public_release/RNA-seq/meiotic_gene_expression/featureCount_clpt/2023-socts.featureCounts.txt"
ctsfile6 <- "/datasets/data_4/nison/hcr3/HCR3_public_release/RNA-seq/meiotic_gene_expression/featureCount_clpt_hei10/Col_hei10_from_clpt-series.featureCounts.txt"
# ctsfile2: HEI10 count
# ctsfile3: hcr2 series

cts <- read_delim(ctsfile1, skip=1, delim="\t")
cts_hei10 <- read_delim(ctsfile2, skip=1, delim="\t")
cts_hcr2 <- read_delim(ctsfile3, skip=1, delim="\t")
cts_hcr2_hei10 <- read_delim(ctsfile4, skip=1, delim="\t")
cts_clpt_col_s <- read_delim(ctsfile5, skip=1, delim="\t") %>%
    dplyr::select(c(1:6, 19:21))
cts_clpt_hei10 <- read_delim(ctsfile6, skip=1, delim="\t")

sample_names <- c("Col_s_A_1", "Col_s_A_2", "Col_s_A_3", "Col_b_A_1", "Col_b_A_2", "Col_b_A_3", "Col_b_B_1", "Col_b_B_2", "Col_b_B_3", "hcr3_b_B_1", "hcr3_b_B_2", "hcr3_b_B_3", "Col_b_C_1", "Col_b_C_2", "Col_b_C_3", "j3_b_C_1", "j3_b_C_2", "j3_b_C_3", "j2_b_C_1", "j2_b_C_2", "j2_b_C_3")
# Set A: hcr11D series
# Set B: hcr3 series (2021)
# Set C: j2, j3 series (2023)

# hcr2 series
sample_names_hcr2 <- c("Col_s_D_1", "Col_s_D_2", "Col_s_D_3", "Col_s_D_4", "Col_b_D_1", "Col_b_D_2", "Col_b_D_3", "Col_b_D_4", "hcr2_b_1", "hcr2_b_2", "hcr2_b_3", "hcr2_b_4")

# Col seedling from clpt series
sample_names_clpt <- c("Col_s_E_1", "Col_s_E_2", "Col_s_E_3")

colnames(cts)[7:ncol(cts)] <- sample_names 
colnames(cts_hei10)[7:ncol(cts_hei10)] <- sample_names 
colnames(cts_hcr2)[7:ncol(cts_hcr2)] <- sample_names_hcr2
colnames(cts_hcr2_hei10)[7:ncol(cts_hcr2_hei10)] <- sample_names_hcr2
colnames(cts_clpt_col_s)[7:ncol(cts_clpt_col_s)] <- sample_names_clpt
colnames(cts_clpt_hei10)[7:ncol(cts_clpt_hei10)] <- sample_names_clpt

# Add Col_s, Col_b, hcr2_b data from hcr2 dataset (set D)
# Add Col_s data from clpt dataset (set E)
cts_wo_hei10 <- left_join(cts, cts_hcr2, by=c("Geneid", "Chr", "Start", "End", "Strand", "Length")) %>%
    left_join(cts_clpt_col_s, by=c("Geneid", "Chr", "Start", "End", "Strand", "Length"))

cts_hei10 <- left_join(cts_hei10, cts_hcr2_hei10, by=c("Geneid", "Chr", "Start", "End", "Strand", "Length")) %>%
    left_join(cts_clpt_hei10, by=c("Geneid", "Chr", "Start", "End", "Strand", "Length"))


cts_total <- cts_wo_hei10 %>%
        filter(Geneid != "AT1G53490") %>%
        bind_rows(cts_hei10) %>%
        dplyr::select(-Chr, -Start, -End, -Strand, -Length) %>%
        column_to_rownames(var="Geneid")

#---- Create 'coldata' containing sample info, which is required for DESeq2 function
# create 'coldata'  containing sample info.
genotype <- c(rep(c("Col_s_A", "Col_b_A", "Col_b_B", "hcr3_b_B", "Col_b_C", "j3_b_C", "j2_b_C"), each=3), rep(c("Col_s_D", "Col_b_D", "hcr2_b"), each=4), rep("Col_s_E", times=3)) # modify 3 
genotype <- factor(genotype)
genotype <- relevel(genotype, ref="Col_s_A")
coldata <- data.frame(genotype=genotype,row.names=colnames(cts_total))
coldata$replicate <- c(rep(c("r1", "r2", "r3"), times=7), rep(c("r1", "r2", "r3", "r4"), times=3), "r1", "r2", "r3") # modify 5

# ...Now inputs for DESeq2 are ready...

#==== B. DESeq2 running ========
deseq <- function(cts, coldata){
    dds <- DESeqDataSetFromMatrix(countData = cts,
                                  colData = coldata,
                                  design = ~ genotype)
    dds
    
    #----pre-filtering--------
    keep <- rowSums(counts(dds)) >= 10
    dds <- dds[keep,]
    
    #----DE analysis--------
    dds <- DESeq(dds)
    return(dds)
}

dds <- deseq(cts_total, coldata) #deseq dataset
ncts <- counts(dds, normalized=TRUE) #normalized counts file (DESeq2's median or ratios method)
write.csv(ncts, file=paste0(outputpath, "/hcr3-paper_full-RNA-seq_normalized_counts.csv"))

#---- QC analysis --------
rld <- rlog(dds, blind=T) # Use rlog transformed data for sample clustering, but not in the other analysis

write.csv(assay(rld), file=file.path(outputpath, "/hcr3-paper_full-RNA-seq_rlog.csv"))
## PCA analysis
dir.create(paste0(outputpath, "/QC"), recursive=TRUE)
samplepca <- plotPCA(rld, intgroup=c("genotype"), returnData=TRUE)
ggpca <- ggplot(samplepca, aes(PC1, PC2, color=genotype)) +
        geom_point(size=1) +
        scale_color_brewer(type="qual", palette="Dark2") +
        theme_bw() +
        theme(text=element_text(size=7, family="Arial", colour="black"), axis.text=element_text(colour="black"), axis.title=element_text(size=7))

pdf(file=paste0(outputpath, "/QC/sample_PCA_anal.pdf"), width=2.25, height=2.25)
print(ggpca)
dev.off()

## Hierarchical clustering of the samples
rld_mat <- assay(rld)
rld_cor <- cor(rld_mat)
samplehc <- pheatmap(rld_cor, silent=TRUE)
pdf(file=paste0(outputpath, "/QC/sample_hc.pdf"), width=4.75, height=4.75)
print(samplehc)
dev.off()

#---- DESeq2 result --------
deseqRes <- function(dat, cgroup, egroup){
		result <- results(dat, contrast=c("genotype", egroup, cgroup), alpha=0.05)
		shres <- lfcShrink(dat, contrast=c("genotype", egroup, cgroup), type="ashr", res=result)
		return(shres)
}
		
## alpha: the significance cutoff used for optimizing the independent filtering (default = 0.1).
## DESeq2 perform independent filtering of low number-reads for efficient multiple testing

## Col_b/Col_s from hcr2 dataset
### HEI10 readcount of Col_s from hcr2 dataset is too high
res_colb_cols_1 <- deseqRes(dds, "Col_s_D", "Col_b_D")

## Col_b/Col_s from hcr11D dataset
## J2 readcount of Col_s from hcr11D data is too high
res_colb_cols_2 <- deseqRes(dds, "Col_s_A", "Col_b_A")

## Col_b from hcr2 dataset, Col_s from clpt series
res_colb_cols_3 <- deseqRes(dds, "Col_s_E", "Col_b_D")
## hcr2/Col bud. Col bud from hcr2 dataset
res_hcr2_colb <- deseqRes(dds, "Col_b_D", "hcr2_b")
## hcr3/Col bud. Col bud from hcr3 dataset (2021)
res_hcr3_colb <- deseqRes(dds, "Col_b_B", "hcr3_b_B")
## j3-3/Col bud. Col bud from j3, j2 dataset (2023)
res_j3_colb <- deseqRes(dds, "Col_b_C", "j3_b_C")
## j2-2/Col bud. Col bud from j3, j2 dataset (2023)
res_j2_colb <- deseqRes(dds, "Col_b_C", "j2_b_C")


dir.create(paste0(outputpath, "/deseq_results"))
write.csv(res_colb_cols_1, file=file.path(outputpath, "deseq_results", "deseqResShrunken_col-b_vs_col-s_from_hcr2.csv"))
write.csv(res_colb_cols_2, file=file.path(outputpath, "deseq_results", "deseqResShrunken_col-b_vs_col-s_from_hcr11-D.csv"))
write.csv(res_colb_cols_3, file=file.path(outputpath, "deseq_results", "deseqResShrunken_col-b-from-_col-s_from_clpt_col-b_from_hcr2.csv"))
write.csv(res_hcr2_colb, file=file.path(outputpath, "deseq_results", "deseqResShrunken_hcr2_vs_col_bud.csv"))
write.csv(res_hcr3_colb, file=file.path(outputpath, "deseq_results", "deseqResShrunken_hcr3_vs_col_bud.csv"))
write.csv(res_j3_colb, file=file.path(outputpath, "deseq_results", "deseqResShrunken_j3_vs_col_bud.csv"))
write.csv(res_j2_colb, file=file.path(outputpath, "deseq_results", "deseqResShrunken_j2_vs_col_bud.csv"))
