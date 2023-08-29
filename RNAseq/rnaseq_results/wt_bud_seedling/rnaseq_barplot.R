library(tidyverse)

readcount <- read_tsv("/datasets/data_4/nison/clpt/RNA-seq/2023_clpt-socts/CLPT1_readCount/CLPT1_bamSummary.tab") %>% dplyr::select(-(1:3))
sizefactor <- read_csv("/datasets/data_4/nison/clpt/RNA-seq/2023_clpt-socts/03_deseq/clpt-socts-2023_sizeFactor.csv")

setwd("result")

colnames(readcount) <- paste(rep(c("Col", "clpt", "clptidm3", "clptros1", "clptnrpd1", "clptrdr2", "idm3", "ros1", "rdr2"), each=3), 1:3, sep="_")
colnames(sizefactor) <- c("Sample", "sizeFactor")
sizefactor$Sample <- paste(rep(c("Col", "clpt", "idm3", "ros1", "clptidm3", "clptros1", "rdr2", "clptnrpd1", "clptrdr2"), each=3), 1:3, sep="_")


readtab <- readcount[1,] %>%
		pivot_longer(cols=1:ncol(.), names_to="Sample",
					 values_to="readCount") %>%
		left_join(sizefactor, by="Sample") %>%
		mutate(normCount = readCount/sizeFactor) %>%
		mutate(genotype = str_replace(Sample, "_.", "")) %>%
		mutate(genotype = factor(genotype, levels=c("Col", "idm3", "ros1", "rdr2", "clpt", "clptidm3", "clptros1", "clptnrpd1", "clptrdr2")))
readtab.summary <- group_by(readtab, genotype) %>%
		summarise(meanCount=mean(normCount), sdCount=sd(normCount))

p <- ggplot() +
		geom_point(data=readtab, aes(x=genotype, y=normCount), size=0.2) +
		geom_errorbar(data=readtab.summary, aes(x=genotype, ymin=(meanCount - sdCount), ymax=(meanCount + sdCount)), width=0.2, colour="red") +
		geom_point(data=readtab.summary, aes(x=genotype, meanCount), colour="red", size=0.4) +
		scale_y_continuous(limits=c(0, 1800)) +
		theme_bw() +
		theme(text=element_text(size=8, family="helvetica", colour="black")) +
		theme(axis.text.x=element_text(angle=30, hjust=1, vjust=1)) +
		labs(title="read count CLPT1 exon 1-3\n(Chr4:12972747-12973305)")

svg(file="CLPT1_read_count.svg", width=3.5, height=2.5)
print(p)
dev.off()

ttest <- function(dat, ctrl, test){
#		ctrl <- "clpt"
#		test <- "clptidm3"
#		dat <- readtab

		dat.ctrl <- filter(dat, genotype == ctrl)$normCount
		dat.test <- filter(dat, genotype == test)$normCount

		print(t.test(dat.ctrl, dat.test))
}

sink(file="CLPT1_readcount_t-test.txt")
print("clptidm3 vs clpt")
ttest(readtab, "clpt", "clptidm3")
print("clptros1 vs clpt")
ttest(readtab, "clpt", "clptros1")
print("clptnrpd1 vs clpt")
ttest(readtab, "clpt", "clptnrpd1")
print("clptrdr2 vs clpt")
ttest(readtab, "clpt", "clptrdr2")
sink()
