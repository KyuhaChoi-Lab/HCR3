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
  library("pheatmap")
  library(org.At.tair.db)
})

#---- output variables --------
outputpath="/datasets/data_4/nison/hcr3/RNA-seq/03_deseq2/"
dir.create(outputpath, recursive=TRUE)

#---- read cts file, a output of STAR --------
ctsfile <- "/datasets/data_4/nison/hcr3/RNA-seq/02_featureCounts/hcr3.featureCounts.txt"
cts <- read_delim(ctsfile, skip=1, delim="\t")

sample_names <- c("Col0_1", "Col0_2", "Col0_3", "hcr3_1", "hcr3_2", "hcr3_3")
colnames(cts)[7:ncol(cts)] <- sample_names 

cts <- cts %>%
        dplyr::select(-Chr, -Start, -End, -Strand, -Length) %>%
        column_to_rownames(var="Geneid")

#---- Create 'coldata' containing sample info, which is required for DESeq2 function
# create 'coldata'  containing sample info.
genotype <- rep(c("Col0", "hcr3"), times=c(3,3)) # modify 3 
genotype <- factor(genotype)
genotype <- relevel(genotype, ref="Col0")
coldata <- data.frame(genotype=genotype,row.names=colnames(cts))
coldata$replicate <- rep(c("r1", "r2", "r3"), times=2) # modify 5

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

dds <- deseq(cts, coldata) #deseq dataset
ncts <- counts(dds, normalized=TRUE) #normalized counts file (DESeq2's median or ratios method)
write.csv(ncts, file=paste0(outputpath, "/hcr3_normalized_counts.csv"))

#---- QC analysis --------
rld <- rlog(dds, blind=T) # Use rlog transformed data for sample clustering, but not in the other analysis

write.csv(assay(rld), file=file.path(outputpath, "/hcr3-bud_rlog.csv"))
## PCA analysis
dir.create(paste0(outputpath, "/QC"), recursive=TRUE)
samplepca <- plotPCA(rld, intgroup=c("genotype"), returnData=TRUE)
ggpca <- ggplot(samplepca, aes(PC1, PC2, color=genotype)) +
        geom_point(size=1) +
        scale_color_brewer(type="qual", palette="Dark2") +
        theme_bw() +
        theme(text=element_text(size=7, family="helvetica", colour="black"), axis.text=element_text(colour="black"), axis.title=element_text(size=7))

svg(file=paste0(outputpath, "/QC/sample_PCA_anal.svg"), width=2.25, height=2.25)
print(ggpca)
dev.off()

## Hierarchical clustering of the samples
rld_mat <- assay(rld)
rld_cor <- cor(rld_mat)
samplehc <- pheatmap(rld_cor, silent=TRUE)
svg(file=paste0(outputpath, "/QC/sample_hc.svg"), width=4.75, height=4.75)
print(samplehc)
dev.off()

#---- DESeq2 result --------
deseqRes <- function(dat, cgroup, egroup){
		result <- results(dat, contrast=c("genotype", egroup, cgroup), alpha=0.05)
		shres <- lfcShrink(dat, contrast=c("genotype", egroup, cgroup), type="ashr", res=result)
		return(shres)
}
		
# tt_c_res <- results(dds, contrast=c("genotype", "tt", "Col0"), alpha=0.05)
## alpha: the significance cutoff used for optimizing the independent filtering (default = 0.1).
## DESeq2 perform independent filtering of low number-reads for efficient multiple testing

res_c_hcr3 <- deseqRes(dds, "Col0", "hcr3")


dir.create(paste0(outputpath, "/deseq_results"))
write.csv(res_c_hcr3, file=file.path(outputpath, "deseq_results", "deseqResShrunken_hcr3_vs_col.csv"))

## MA plot of DESeq2 result
pdf(file=paste0(outputpath, "/MA_plot.pdf"))
hcr3_c_ma <- plotMA(res_c_hcr3, ylim=c(-2, 2), main="hcr3 vs Col0")
dev.off()

#---- DEG call --------

## Set thresholds
padj=0.05

callDEG <- function(res, padj.cutoff, lfc.cutoff){
        deg <- res %>%
                data.frame() %>%
                rownames_to_column(var="gene") %>%
                as_tibble() %>%
                filter(padj < padj.cutoff & abs(log2FoldChange) > lfc.cutoff)
        return(deg)
}

deg_c_hcr3 <- callDEG(res_c_hcr3, padj,0.5)
deg_c_hcr3_lfc1 <- callDEG(res_c_hcr3, padj, 1)
deg_c_hcr3_lfc0 <- callDEG(res_c_hcr3, padj, 0)

deg_c_hcr3$symbol <- mapIds(org.At.tair.db, keys=deg_c_hcr3$gene, keytype="TAIR", column="SYMBOL")
deg_c_hcr3_lfc1$symbol <- mapIds(org.At.tair.db, keys=deg_c_hcr3_lfc1$gene, keytype="TAIR", column="SYMBOL")
deg_c_hcr3_lfc0$symbol <- mapIds(org.At.tair.db, keys=deg_c_hcr3_lfc0$gene, keytype="TAIR", column="SYMBOL")

dir.create(paste0(outputpath, "/deg_call"))
write.csv(deg_c_hcr3, file=file.path(outputpath, "deg_call", "fdr005_lfc05_hcr3_vs_col.csv"))
write.csv(deg_c_hcr3_lfc0, file=file.path(outputpath, "deg_call", "fdr005_lfc00_hcr3_vs_col.csv"))
write.csv(deg_c_hcr3_lfc1, file=file.path(outputpath, "deg_call", "fdr005_lfc10_hcr3_vs_col.csv"))

## DEG call report
degreport <- function(deg){
        tot <- count(deg)[[1]]
        inc <- count(deg[deg$log2FoldChange>0,])[[1]]
        dec <- count(deg[deg$log2FoldChange<0,])[[1]]
        rep <- c(total=tot, up=inc, down=dec)
        return(rep)
}

deg_summary <- bind_rows(degreport(deg_c_hcr3)) %>%
			add_column(group=c("hcr3_vs_c"), .before="total")
write.csv(deg_summary, file=file.path(outputpath, "DESeq2_lfc05_DEG_summary.csv"))



# --- fpkm ---

library("GenomicFeatures")
library("TxDb.Athaliana.BioMart.plantsmart28")
txdb <- TxDb.Athaliana.BioMart.plantsmart28
exonsByGene <- exonsBy(txdb, by="gene")
cts_fpkm <- cts[rownames(cts) %in% names(exonsByGene),]
exonsByGene_filter <- exonsByGene[rownames(cts_fpkm)]

dds_fpkm <- DESeqDataSetFromMatrix(countData=cts_fpkm,
								   colData=coldata,
								   rowData=exonsByGene_filter,
								   design= ~ genotype)
dds_fpkm <- estimateSizeFactors(dds_fpkm)

all_fpkm <- fpkm(dds_fpkm, robust=TRUE)
									
write.csv(all_fpkm, file=file.path(outputpath, "hcr3-bud_fpkm.csv"))
