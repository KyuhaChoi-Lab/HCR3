library(tidyverse)
library(khroma)

dir.create("results", recursive=T)

pal <- color("bright")(7)[c("blue", "green")]
names(pal) <- c("seedling", "bud")
goi <- c("J2"="AT5G22060", "J3"="AT3G44110")
readcount <- read_csv("hcr2-bud-seedling_deseq-normalized.csv") %>%
    dplyr::select(1:9) %>%
    filter(...1 %in% goi) %>%
    pivot_longer(cols=2:ncol(.), names_to="sample", values_to="normCount") %>%
    separate(col=sample, into=c("col", "sample", "rep"), sep="_") %>%
    dplyr::select(-col)

colnames(readcount)[1] <- "agi"
readcount <- mutate(readcount, gene = case_when(agi == goi[1] ~ names(goi[1]),
                                                agi == goi[2] ~ names(goi[2]))) %>%
            mutate(gene = factor(gene, levels=names(goi))) %>%
            mutate(sample = case_when(sample == "b" ~ "bud",
                                    sample == "s" ~ "seedling")) %>%
            mutate(sample = factor(sample, levels=c("seedling", "bud")))

readcount.summary <- group_by(readcount, gene, sample) %>%
    summarise(meanCount = mean(normCount),
            sdCount = sd(normCount))

write_csv(readcount, file="results/j2_j3_normalized_readcount.csv")

# read count bar plots

pBar <- ggplot() +
        geom_col(data=readcount.summary, aes(x=gene, y=meanCount, group=sample), position=position_dodge(width=0.9), colour="black", fill=NA, width=0.8) +
		geom_jitter(data=readcount, aes(x=gene, y=normCount, group=sample, fill=sample), shape=21, colour="black", size=1, position=position_jitterdodge(dodge.width=0.9, jitter.width=0.35, jitter.height=0)) +
		geom_errorbar(data=readcount.summary, aes(x=gene, ymin=(meanCount - sdCount), ymax=(meanCount + sdCount), group=sample), width=0.2, colour="black", position=position_dodge(width=0.9)) +
        scale_fill_manual(values=pal) +
        scale_colour_manual(values=pal) +
		# geom_point(data=readtab.summary, aes(x=genotype, meanCount), colour="red", size=0.4) +
		# scale_y_continuous(limits=c(0, 1800)) +
		theme_classic() +
		theme(text=element_text(size=8, family="Helvetica", colour="black")) +
        labs(y="Normalized Count") +
        theme(axis.title.x = element_blank())
        # theme(legend.position=c(0.3,0.9))
        # theme(legend.key.size = unit(0.1, "in")) +
        # guides(
            # colour = guide_legend(override.aes = list(size=2))
        # )
		# theme(axis.text.x=element_text(angle=30, hjust=1, vjust=1))
		# labs(title="read count CLPT1 exon 1-3\n(Chr4:12972747-12973305)")

addSmallLegend <- function(myplot, pointSize=0.5, textSize=3, spaceLegend = 0.1) {
    myplot +
        guides(shape = guide_legend(override.aes = list(size = pointSize)),
                color = guide_legend(override.aes = list(size = pointSize))) +
        theme(legend.title = element_text(size = textSize),
                legend.text = element_text(size = textSize),
                legend.key.size = unit(spaceLegend, "lines"))
}

pBar.smalllegend <- addSmallLegend(pBar, pointSize=1, textSize=7, spaceLegend = 0.2)
pdf(file="results/j2_j3_expression_normcount.pdf", width=2.5, height=1.5)
print(pBar.smalllegend)
dev.off()
png(file="results/j2_j3_expression_normcount.png", width=2.5, height=1.5, unit="in", res=300)
print(pBar.smalllegend)
dev.off()

# readcount Crossbar 
pCrossbar <- ggplot() +
		geom_jitter(data=readcount, aes(x=gene, y=normCount, group=sample, fill=sample), shape=21, colour="black", size=1, position=position_jitterdodge(dodge.width=0.9, jitter.width=0.35, jitter.height=0)) +
        geom_point(data=readcount.summary, aes(x=gene, y=meanCount, group=sample), shape="\U2014", size=3, colour="black", position=position_dodge(width=0.9)) +
        # stat_summary(fun="mean", geom="segment", mapping=aes(xend=..x.. - 0.25, yend=..y..)) +
        # stat_summary(fun="mean", geom="segment", mapping=aes(xend=..x.. - 0.25, yend=..y..)) +
        # geom_crossbar(data=readcount.summary, aes(x=gene, ymin=meanCount, ymax=meanCount, y=meanCount, group=sample), position=position_dodge(width=0.9), width = 0.5, linewidth=0.5) +
		geom_errorbar(data=readcount.summary, aes(x=gene, ymin=(meanCount - sdCount), ymax=(meanCount + sdCount), group=sample), width=0.2, colour="black", position=position_dodge(width=0.9)) +
        scale_fill_manual(values=pal) +
        scale_colour_manual(values=pal) +
		# geom_point(data=readtab.summary, aes(x=genotype, meanCount), colour="red", size=0.4) +
		# scale_y_continuous(limits=c(0, 1800)) +
		theme_classic() +
		theme(text=element_text(size=8, family="Helvetica", colour="black")) +
        labs(y="Normalized Count") +
        theme(axis.title.x = element_blank())
        # theme(legend.position=c(0.3,0.9))
        # theme(legend.key.size = unit(0.1, "in")) +
        # guides(
            # colour = guide_legend(override.aes = list(size=2))
        # )
		# theme(axis.text.x=element_text(angle=30, hjust=1, vjust=1))
		# labs(title="read count CLPT1 exon 1-3\n(Chr4:12972747-12973305)")

pCrossbar.smalllegend <- addSmallLegend(pCrossbar, pointSize=1, textSize=7, spaceLegend = 0.2)
pdf(file="results/j2_j3_expression_normcount_crossbar.pdf", width=2.5, height=1.5)
print(pCrossbar.smalllegend)
dev.off()
png(file="results/j2_j3_expression_normcount_crossbar.png", width=2.5, height=1.5, unit="in", res=300)
print(pCrossbar.smalllegend)
dev.off()

# read count boxplots

pBox <- ggplot() +
        geom_boxplot(data=readcount, aes(x=gene, y=normCount, fill=sample), position=position_dodge(width=0.9), colour="black", width=0.8) +
		# geom_jitter(data=readcount, aes(x=gene, y=normCount, group=sample, fill=sample), shape=21, colour="black", size=1, position=position_jitterdodge(dodge.width=0.9, jitter.width=0.35, jitter.height=0)) +
		# geom_errorbar(data=readcount.summary, aes(x=gene, ymin=(meanCount - sdCount), ymax=(meanCount + sdCount), group=sample), width=0.2, colour="black", position=position_dodge(width=0.9)) +
        scale_fill_manual(values=pal) +
        scale_colour_manual(values=pal) +
		# geom_point(data=readtab.summary, aes(x=genotype, meanCount), colour="red", size=0.4) +
		# scale_y_continuous(limits=c(0, 1800)) +
		theme_classic() +
		theme(text=element_text(size=8, family="Helvetica", colour="black")) +
        labs(y="Normalized Count") +
        theme(axis.title.x = element_blank())
        # theme(legend.position=c(0.3,0.9))
        # theme(legend.key.size = unit(0.1, "in")) +
        # guides(
            # colour = guide_legend(override.aes = list(size=2))
        # )
		# theme(axis.text.x=element_text(angle=30, hjust=1, vjust=1))
		# labs(title="read count CLPT1 exon 1-3\n(Chr4:12972747-12973305)")

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
