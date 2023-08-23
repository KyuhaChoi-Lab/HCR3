#########################################################################
# THIS IS THE TEMPLATE FOR GBS ANALYSIS.
# SPECIFY INPUT MATERIALS, WHICH ARE *_cotable.csv, AT THE BEGINNING.
# RUN THROUGH THE SCRIPT, THEN ADJUST PLOTS ACCORDINGLY.
# (*.cotable are product of GBS_cotable.R script.)
#
#@ maintainer: Namil Son
#@ date: 2021-04-21
#

##### important note! #####
# 2022.01.14
# In TEL-CEN analysis, north arm of Chr2 and Chr4 are too short for fair comparison with other chromosomes. So I excluded left arms of Chr2 and Chr4 from the analysis
# Also, Chr3 is largely distorted in pericentromeric region probably because of genome rearrangement induced by introducint meiMIGS construct. So I also excluded this region from the TEL-CEN analysis
###########################

#########################################################################

library(tidyverse)
library(colorspace)
library(scales)
library(khroma)
library(lemon)

#' SETTINGS =============================================================
#@ 1. Read samples
dir_cotable="cotables"
file1=file.path(dir_cotable, "WT_2019_cotable.txt")
file2=file.path(dir_cotable, "WT_2020_cotable.txt")
file3=file.path(dir_cotable, "WT_2021_cotable.txt")
file4=file.path(dir_cotable, "20230811_GBS_J3_mJ3_cotable.txt")
file5=file.path(dir_cotable, "20230727_GBS_pSPO11-1_mJ3_cotable.txt")

readCotable <- function(file) {
        dat <- read_csv(file) %>%
                mutate(lib = str_replace(lib, "results/06_tiger/lib", "")) %>%
                mutate(lib = str_replace(lib, "_MappedOn_tair10", "")) %>%
                mutate(lib = as.numeric(lib))
        return(dat)
}

sample1=readCotable(file1)
sample2=readCotable(file2)
sample3=readCotable(file3)
sample4=readCotable(file4)
sample5=readCotable(file5)

#@ 2. Merge cotables from the same genotype
## In case cotable from more than one sample needs to be merged for a genotype,
## merge those files and edit the lib number accordingly to avoid duplicated lib numbers.
## Then put that merged cotable file into genotype[n] variables 
sample2$lib <- sample2$lib + max(sample1$lib)
sample3$lib <- sample3$lib + max(sample2$lib)
wt=bind_rows(sample1, sample2, sample3)
mj3_2=filter(sample4, lib <=48)
mj3_3=filter(sample4, lib <=96 & lib >48)
pspo11=sample5

## genotypes of the samples
genotypes=c("WT", "mJ3-2", "mJ3-3", "pSPO11_mJ3")

#@ 3. colors for each genotype
# below is for four genotypes
light <- colour("light")(9)
pal <- light[c("light blue", "pear", "olive", "orange")] 
names(pal) <- genotypes

#@ 4. Merge the cotables from different genotypes
cos.all.list <- list(wt, mj3_2, mj3_3, pspo11)
names(cos.all.list) <- genotypes
cos.all <- bind_rows(cos.all.list, .id="genotype") %>%
        mutate(genotype = factor(genotype, levels=genotypes))

## count the number of crossover for each individual. 
cos.all.count <- bind_rows(cos.all.list, .id="genotype") %>%
        group_by(genotype) %>%
        count(lib) %>%
        mutate(genotype = factor(genotype, levels=genotypes))
#@ 5. prefix and dirout
prefix="independent_mj3"
dirout=paste0("results")
# setwd(dirout)

#@ 7. Col-0 DNA methylation binnded by 200k/100k base
col0lang15.bin200 <- read_csv("data/col0Lang15_window2e+05_C.csv", col_names = TRUE)
colnames(col0lang15.bin200) <- c("chrs", "bin.start", "bin.end", "width", "strand", "C", "mC", "bin")
col0lang15.bin100 <- read_tsv("data/SRR11304866_C_TAIR10_window100kb_step100kb.tsv", col_names = TRUE)
colnames(col0lang15.bin100) <- c("chrs", "bin.start", "bin.end", "bin")
#=======================================================================================

#' 0. chromosome info 
centromeres <- c(15086045,3607929,13587786,3956021,11725024)
chrs <- c("Chr1","Chr2","Chr3","Chr4","Chr5")
chr.ends <- c(30427671,19698289,23459830,18585056,26975502)

tha.cum <- cumsum(chr.ends)
# cumulative coordinate of chromosomes 
tha.cum <- c(0,tha.cum)
# total length of arabidopsis nuclear genomes
tha.tot <- tha.cum[length(tha.cum)]
# centromere coordinate on the cumulative coordinates
cent <- centromeres
centromeres[2] <- centromeres[2]+tha.cum[2]
centromeres[3] <- centromeres[3]+tha.cum[3]
centromeres[4] <- centromeres[4]+tha.cum[4]
centromeres[5] <- centromeres[5]+tha.cum[5]

# chromosome length. same as chr.ends
chr.lens <- c(30427671,19698289,23459830,18585056,26975502)
chrs <- c("Chr1","Chr2","Chr3","Chr4","Chr5")

# pericentromeric + centromeric region
# (Ref.: Underwood, 2018, Genome Research)
north.start <- c(11420001, 910001, 10390001, 1070001, 8890001)
south.end <- c(18270000, 7320000, 16730000, 6630000, 15550000)

#' 1. Statistical tests of COs per F2 ========
wt.count <- filter(cos.all.count, genotype==genotypes[1])$n
mj3_2.count <- filter(cos.all.count, genotype==genotypes[2])$n
mj3_3.count <- filter(cos.all.count, genotype==genotypes[3])$n
pspo11.count <- filter(cos.all.count, genotype==genotypes[4])$n
count.list <- list(wt.count, mj3_2.count, mj3_3.count, pspo11.count)

con <- file(file.path(dirout, paste0(prefix, "_CO_statistical-test.txt")))

sink(con, append=TRUE)
# students t-test
ttest_mat <- matrix(, nrow=length(count.list)-1, ncol=length(count.list))
rownames(ttest_mat) <- genotypes[1:3]
colnames(ttest_mat) <- genotypes[1:4]
for (i in 1:(length(count.list)-1)){
        ctrl_group <- count.list[[i]]
        for (j in i:length(count.list)) {
                test_group <- count.list[[j]]
                ttest_mat[i,j] <- t.test(ctrl_group, test_group)$p.value
        }
}
print(ttest_mat)
sink()


#' 2. histogram of COs per F2 ============
# calculate mean co
mean_co <- cos.all.count %>%
        group_by(genotype) %>%
        summarise(co=mean(n))

hist_bin <- max(cos.all.count$n)-1
co_hist <- ggplot(cos.all.count, aes(x=n, after_stat(density))) +
        geom_histogram(bins=hist_bin, position="dodge", fill="white", colour="black") +
        geom_vline(data=mean_co, aes(xintercept=co), linetype="dashed", colour="red") +
        facet_rep_grid(genotype ~ ., repeat.tick.labels="all") +
        scale_y_continuous(breaks=c(0, 0.1, 0.2)) +
        # facet_grid(rows="genotype") +
        # scale_fill_manual(values=pal) +
        # scale_colour_manual(values=pal) +
        theme_classic() +
        theme(text=element_text(size=9, family="Helvetica", colour="black"), axis.text=element_text(colour="black")) +
        labs(x="Crossovers", y="Ratio") +
        theme(legend.key.size=unit(0.1, "inches"),
              legend.title=element_text(size=7),
              legend.text=element_text(size=7),
              legend.position="top") +
        theme(strip.background = element_rect(colour=NA, fill = "white")) +
        theme(axis.line = element_line(colour="black", linewidth=0.5))
        # theme(panel.border = element_rect(colour="black", fill=NA, linewidth=1)) +
        # theme(axis.line = element_line(colour="black", linewidth=0))

cos.all.count.stick <- cos.all.count %>%
    group_by(genotype, n) %>%
    summarise(freq=n()) %>%
    mutate(ratio = freq/sum(freq))

co_hist_stick <- ggplot(cos.all.count.stick, aes(x=n, y=ratio)) +
        # geom_histogram(bins=hist_bin, position="dodge", fill="white", colour="black") +
        geom_col(width=0.2, fill="black") +
        geom_vline(data=mean_co, aes(xintercept=co), linetype="dashed", colour="red") +
        facet_rep_grid(genotype ~ ., repeat.tick.labels="all") +
        scale_y_continuous(breaks=c(0, 0.1, 0.2)) +
        # facet_wrap(vars(genotype), ncol=1) +
        # facet_grid(rows="genotype") +
        # scale_fill_manual(values=pal) +
        # scale_colour_manual(values=pal) +
        theme_classic() +
        theme(text=element_text(size=9, family="Helvetica", colour="black"), axis.text=element_text(colour="black")) +
        labs(x="Crossovers", y="Ratio") +
        theme(legend.key.size=unit(0.1, "inches"),
              legend.title=element_text(size=7),
              legend.text=element_text(size=7),
              legend.position="top") +
        theme(strip.background = element_rect(colour=NA, fill = "white")) +
        theme(axis.line = element_line(colour="black", linewidth=0.5))
        # theme(panel.border = element_rect(colour="black", fill=NA, linewidth=1)) +


pdf(file=file.path(dirout, paste0(prefix, "_co_hist.pdf")), width=2, height=3.5)
co_hist
dev.off()
png(file=file.path(dirout, paste0(prefix, "_co_hist.png")), width=2, height=3.5, unit="in", res=300)
co_hist
dev.off()

pdf(file=file.path(dirout, paste0(prefix, "_co_hist_stick.pdf")), width=2, height=3.5)
co_hist_stick
dev.off()
png(file=file.path(dirout, paste0(prefix, "_co_hist_stick.png")), width=2, height=3.5, unit="in", res=300)
co_hist_stick
dev.off()


                

#' 3. crossover summary table =============
cos.width <- cos.all %>%
        group_by(genotype) %>%
        summarise(mean.width=mean(width))
cos.width <- cos.width$mean.width

nlib <- group_by(cos.all.count, by=genotype) %>%
        summarise(nlib=length(unique(lib))) %>%
        dplyr::select(nlib) %>%
        as_vector()
        
cos_summary <- cos.all %>%
        group_by(genotype, chrs) %>%
        summarise(ncos=dplyr::n()) %>%
        spread(key=chrs, value=ncos) %>%
        add_column(n.lib=nlib) %>%
        mutate(n.cos=sum(Chr1, Chr2, Chr3, Chr4, Chr5)) %>%
        mutate('n.cos/n.lib'=n.cos/n.lib) %>%
        add_column(cos.width=cos.width)

write_csv(cos_summary, file="cos.all.counts_summary.csv")



#' 4. Draw crossovers per chromosome plot
mean_cos_chr <- cos_summary[,1:7] %>%
        mutate(across(Chr1:Chr5, ~.x/n.lib)) %>%
        # mutate_at(vars(-c(genotype, n.lib)), funs(./n.lib)) %>%
        pivot_longer(cols=2:6, names_to="chrs", values_to="mean_co") %>%
        add_column(chr.size=rep(chr.ends, times=length(genotypes))) %>%
        mutate(cMMb = mean_co/(chr.size*2/1e6)*100)
ylim.co_by_chr <- c(1, max(mean_cos_chr$cMMb) * 1.2)

p <- ggplot(mean_cos_chr, aes(x=chrs, y=cMMb, colour=genotype)) + 
        # geom_point(aes(group=genotype), size=2) +
        geom_point() +
        theme_classic() +
        scale_colour_manual(values=pal) +
        # scale_x_continuous(name="Chromosome length (Mb)", labels=scales::label_number(scale=1/1000000, accuracy=1), breaks=seq(20000000, 30000000, by=5000000)) +
        labs(x="Chromosome", y="Crossovers (cM/Mb)")  +
        theme(text=element_text(size=9, family="Helvetica", colour="black"), axis.text=element_text(colour="black")) +
        theme(legend.key.size=unit(0.15, "inches"),
                axis.title.x = element_blank(),
              legend.position="top",
              legend.title=element_text(size=7),
              legend.text=element_text(size=7)) + 
		ylim(ylim.co_by_chr)

pdf(file.path(dirout, paste0(prefix, "_mean_co_per_chrlen.pdf")), width=2.75, height=2.25)
print(p)
dev.off()
png(file.path(dirout, paste0(prefix, "_mean_co_per_chrlen.png")), width=2.75, height=2.25, unit="in", res=300)
print(p)
dev.off()

## statistical test for each chromosome
# genotype1.chr <- group_by(genotype1, lib, chrs) %>%
# 		summarise(nco=n())
# genotype2.chr <- group_by(genotype2, lib, chrs) %>%
# 		summarise(nco=n())
# genotype3.chr <- group_by(genotype3, lib, chrs) %>%
# 		summarise(nco=n())

# t_test.chr <- function(chr){
# 		geno2 <- filter(genotype2.chr, chrs==chr)$nco
# 		geno3 <- filter(genotype3.chr, chrs==chr)$nco
# 		result <- t.test(geno2, geno3)
# 		print(result)
# 		}

# t_test_by_chr <- file(paste0(prefix, "_t-test_by_chr.txt"))
# sink(t_test_by_chr, append=TRUE)
# print("Chr1")
# t_test.chr1 <- t_test.chr("Chr1")
# print("Chr2")
# t_test.chr1 <- t_test.chr("Chr2")
# print("Chr3")
# t_test.chr1 <- t_test.chr("Chr3")
# print("Chr4")
# t_test.chr1 <- t_test.chr("Chr4")
# print("Chr5")
# t_test.chr1 <- t_test.chr("Chr5")
# sink()


#' 5. arm versus pericentromere locations =============
# centromere, pericentromere coordinates
cen.pericen <- tibble(Chr=1:5, chr.start=1, north.start=north.start, south.end=south.end, chr.ends=chr.ends)

# armperi(): count the number of crossovers within the centromeric/pericentromeric region or chromosomal arm
armperi <- function(cotable){
        #cotable <- hcr1a.cos
        narm <- NULL
        nperi <- NULL
       for(i in 1:5){
               chr <- cotable %>%
                       filter(chrs==paste0("Chr",i))
               arm <- chr %>%
                       filter((cos>cen.pericen$chr.start[i] & cos<cen.pericen$north.start[i]) | (cos>cen.pericen$south.end[i] & cos<cen.pericen$chr.ends[i])) %>%
                       nrow()
               narm <- c(narm, arm)

               peri <- chr %>%
                       filter(cos>cen.pericen$north.start[i] & cos<cen.pericen$south.end[i]) %>%
                       nrow()
               nperi <- c(nperi, peri)
       }
narmperi <- tibble(chr=1:5, arm=narm, peri=nperi)
narmperi <- narmperi %>%
        add_row(summarise_all(narmperi, sum))
return(narmperi)
}

# apply armperi() to each crossover tables

cos.armperi <- lapply(cos.all.list, armperi)

# chisq.tests
cos.armperi.chisq <- bind_rows(cos.armperi) %>%
        filter(chr==15) %>%
        dplyr::select(-chr)
con <- file(paste0(prefix, "_armperi_chisq-test.txt"))
sink(con, append=TRUE)
chisq.test(cos.armperi.chisq)
sink()

#' 6 plotting crossover/SNP/DNA meth chromatin landscale

# bin_chr(): slice chromosome by window size and record coordinate and cumulative coordinate of each slice
bin_chr <- function(winsize){
        wg.bin <- tibble(chrs=character(), bin.start=numeric(), bin.end=numeric(), cum.start=numeric(), cum.end=numeric())
        for(i in 1:5){
                coord.start <- seq(1, chr.ends[i], by=winsize)
                coord.end <- c(coord.start[-1]-1, chr.ends[i])
                cum.start <- coord.start+tha.cum[i]
                cum.end <- coord.end + tha.cum[i] 
                nchr <- rep(i, length(coord.start))
                bins <- tibble(chrs=paste0("Chr",nchr), bin.start=coord.start, bin.end=coord.end, cum.start=cum.start, cum.end=cum.end)
                wg.bin <- add_row(wg.bin, bins)
        }
        return(wg.bin)
        #wg.bin = whole genome bins
}

bin200k <- bin_chr(200000)

# bin_costable: slice cotable by the product of bin_chr(), then count the number of COS within that slice
## cotable = a product of cotable() function
# binning: slice the data by the bins produced by bin_chr(), then count the number of data within that bin. Can be used for binning CO, SNP, DNAme.
## dat_col: (string); the name of column where the data coordinate is
## wg.bin = whole genome bin; a product of bin_chr() function
binning <- function(dat, dat_col, wg.bin){
        binned <- NULL
        collect.bin <- NULL
        for(i in 1:5){
                chr.dat <- dat %>%
                        filter(chrs==paste0("Chr",i))
                chr.bin <- wg.bin %>%
                        filter(chrs==paste0("Chr",i))
                for(j in 1:nrow(chr.bin)){
                        dat.in.bin <- sum(chr.dat[[dat_col]] >= chr.bin$bin.start[j] & chr.dat[[dat_col]] <= chr.bin$bin.end[j])
                        collect.bin <- c(collect.bin, dat.in.bin)
                }
        }
        #fin <- tibble(chrs=wg.bin$chrs, bin=collect.bin)
        fin <- tibble(wg.bin, bin=collect.bin)

        if(dat_col == "cos"){
                fin <- fin %>%
                        mutate(across(bin, ~./length(unique(dat$lib))))
        }

        return(fin)
}

#col0.co.bin <- binning(cos.all.list[[1]], "cos", bin200k)
#hcr2.co.bin <- binning(hcr2.co, "cos", bin200k, "hcr2")

cos.all.list.bin <- lapply(cos.all.list, binning, "cos", bin200k)

#col0.co.bin$bin <- col0.co.bin$bin/length(unique(col0.co$lib))
#hcr2.co.bin$bin <- hcr2.co.bin$bin/length(unique(hcr2.co$lib))

# Read Col0 mC tiled by 200kb.
col0lang15.bin200 <- col0lang15.bin200 %>%
        filter(!(chrs %in% c("ChrC", "ChrM"))) %>%
        left_join(bin200k)  %>%
        add_column(group="mC") %>%
        select(-c("strand", "C", "mC", "width"))

# mafilter(): linear filtering of binned crossover table. Apply moving average.
## set ma filter
## k=one-side width of moving average filter

mafilter <- function(dat, k){
        filt <- 1/(2*k+1)
        filt <- rep(filt, 2*k+1)
        filt.collect <- NULL
        for(i in 1:5){
                chr <- dat %>%
                        filter(chrs==paste0("Chr",i))
                filt.chr <- stats::filter(chr$bin, filt)
                filt.chr[1:k] <- filt.chr[k+1]
                len <- length(filt.chr)
                filt.chr[(len-k+1):len] <- filt.chr[(len-k)]

                filt.collect <- c(filt.collect, filt.chr)
        }
        fin <- tibble(dat, filt.bin=filt.collect)
        return(fin)
}


cos.all.list.bin.filt <- lapply(cos.all.list.bin, mafilter, 5)
cos.all.list.bin.filt.bind <- bind_rows(cos.all.list.bin.filt, .id="group")

col0lang15.bin200.filt <- mafilter(col0lang15.bin200, 5)

mafilt <- bind_rows(cos.all.list.bin.filt.bind, col0lang15.bin200.filt) 


# Draw chromosome plots
## color setting
## centromere : "Black"
## chromosome end : "Black"
## snp density : "Dark2"
## p.co_chr : drawing CO landscape on the background of SNP density
## p.co_meth_chr : drawing CO landscape on the background of mC density
pri.ymax <- max(filter(mafilt, group != "mC")$filt.bin) * 1.1

interval.len <- mafilt[1,]$bin.end - mafilt[1,]$bin.start +1
toCm <- 1e6 / (2*interval.len) * 100
pri.ymax.cm <- pri.ymax * toCm

p.co_meth_chr <- function(dat){
        p <- ggplot() +
                geom_line(data=filter(dat, !(group %in% c("snp", "mC"))), aes(x=cum.start, y=filt.bin * toCm, colour=group), size=0.4) +
                geom_area(data=filter(dat, group=="mC"), aes(x=cum.start, y=filt.bin*pri.ymax.cm/0.25), fill="grey60", alpha=0.5) +
                scale_y_continuous(name="Crossovers (cM/Mb)", sec.axis=sec_axis(~.*0.25/pri.ymax.cm, name="DNA methylation (mC/C)", breaks=c(0, 0.1, 0.2)), breaks=seq(0, pri.ymax.cm, by=5), limits=c(0, pri.ymax.cm)) +
                scale_x_continuous(name="Coordinates (Mb)", labels=scales::label_number(scale=1/1000000), breaks=c(seq(1, max(dat$cum.start), 20*10^6), 120000000)) +
                scale_colour_manual(values=pal) +
                geom_vline(xintercept=tha.cum, colour="Black", size=0.2) +
                geom_vline(xintercept=centromeres, colour="Black", linetype="dashed", size=0.2) +
                theme_classic() +
                theme(legend.key.size=unit(0.2, "inches"),
                legend.position="top",
                legend.title=element_text(size=7),
                legend.text=element_text(size=7)) +
                theme(text=element_text(size=9, family="Helvetica", colour="black"),
                axis.text=element_text(colour="black", size=7))

        return(p)
}

p.co_meth_land <- p.co_meth_chr(mafilt)

pdf(file=file.path(dirout, paste0(prefix, "_co-meth_landscape_ma5.pdf")), width=8, height=2)
p.co_meth_land
dev.off()
png(file=file.path(dirout, paste0(prefix, "_co-meth_landscape_ma5.png")), width=8, height=2, unit="in", res=300)
p.co_meth_land
dev.off()

#' 7 TEL-CEN plotting
# distFromTel() : calculate the distance from telomere in two directions (north or south).
# It takes *.cotable as input

#-----Note-----
# This code is kind of trikcy at a first glance.
# The first idea of TEL-CEN scaling that pop in my mind is below:
#   1. scale the CO coordinate into proportion
#   2. tile the proportion into bin
#   However, this does not work because when you scale the coordinate before tiling,
# you will eventually mixup the CO coordinates completely. It hides the pattern.
#   Therefore, before converting into proportion, just tile first so that the CO sites fromthe same chromosome can be tied together a little. Maybe it's the point where the method needs to be modified. 
#---------------
bin100k <- bin_chr(100000)

cos.all.list.bin100 <- lapply(cos.all.list, binning, "cos", bin100k)

# read col0lang15 methylation
col0lang15.bin100 <- col0lang15.bin100 %>%
        mutate(bin.end = bin.start + 100000 -1)
        filter(!(chrs %in% c("ChrC", "ChrM"))) %>%
        left_join(bin100k)  %>%
        add_column(group="mC")

# distFromTel() : convert coordinate into proportion of distance from cent to tel
distFromTel <- function(dat, binsize, mafilt.size){
		# dat <- cos.all.list.bin100$hcr2
		# binsize <- 0.01
		# mafilt.size <- 9

        # dat <- col0lang15.bin100
#        dat <- arrange(dat, chrs)
        dat$cent <- rep(cent, times=table(dat$chrs))
        dat$end <- rep(chr.ends, times=table(dat$chrs))
        
        left <- filter(dat, bin.start<cent) %>%
                add_column(arm="north") %>%
                mutate(coord.prop=bin.start/cent) %>%
				filter(chrs %in% c("Chr1", "Chr3", "Chr5"))
		# Excluded north arm of Chr2 and Chr4 as these regions is too short for fair comparison with other chromosomes
        right <- filter(dat, bin.start>=cent) %>%
                add_column(arm="south") %>%
                mutate(coord.prop=(end-bin.start)/(end-cent))
        
        prop <- bind_rows(left, right)
		# prop <- filter(prop, chrs != "Chr3")
		# exclude chr3 from the TEL-CEN analysis
        
        collect_bin <- NULL
        wins <- seq(0, 1, by=binsize)
        for(j in 1:(length(wins)-1)){
                bin <- mean(prop[which(prop$coord.prop >= wins[j] & prop$coord.prop < wins[j+1]),]$bin)
                collect_bin <- c(collect_bin, bin)
        }

        k=mafilt.size

        filt <- 1/(2*k+1)
        filt <- rep(filt, 2*k+1)
        filt.dat <- stats::filter(collect_bin, filt)

        res <- tibble(prop=wins[1:100],
                      bin=filt.dat,
                      )
        
        return(res)
}


cos.all.list.bin100 <- lapply(cos.all.list, binning, "cos", bin100k)
cos.all.list.bin100.prop <- lapply(cos.all.list.bin100, distFromTel, 0.01, 9)
cos.all.list.bin100.prop.bind <- bind_rows(cos.all.list.bin100.prop, .id="group")

col0lang15.meth.prop <- distFromTel(col0lang15.bin100, 0.01, 9)
col0lang15.meth.prop$group="mC"



prop_co <- bind_rows(cos.all.list.bin100.prop.bind, col0lang15.meth.prop)
mean_cos <- prop_co %>%
        group_by(group) %>%
        summarise(meanco=mean(bin, na.rm=TRUE)) %>%
        filter(!(group %in% c("snp", "mC")))

telcen.ymax <- max(filter(prop_co, !(group %in% c("snp", "mC") | is.na(bin)))$bin)*1.1
telcen.ymax.cm <- telcen.ymax * toCm * 2

pTelCen.meth.ma9 <- ggplot() + 
        geom_line(data=filter(prop_co, !(group %in% c("snp", "mC"))), aes(x=prop, y=bin * toCm *2, colour=group), size=0.4) +
        geom_area(data=filter(prop_co, group=="mC"), aes(x=prop, y=bin*telcen.ymax.cm/0.2), fill="grey60", alpha=0.5) +
        geom_hline(data=mean_cos, aes(yintercept=meanco * toCm * 2, colour=group), linetype="dashed", size=0.3) +
        scale_y_continuous(name="Crossovers (cM/Mb)", sec.axis=sec_axis(~.*0.2/telcen.ymax.cm, name="DNA methylation (mC/C)")) +
        scale_x_continuous(name="Distance from the telomere (ratio)",
                           breaks=seq(0, 1, by=0.25),
                           label=c("TEL", 0.25, 0.5, 0.75, "CEN"))+
        scale_colour_manual(values=pal) +
        theme_classic() +
        theme(legend.key.size=unit(0.1, "inches"),
              legend.title=element_text(size=7),
              legend.text=element_text(size=7),
              legend.position="top") +
        theme(text=element_text(size=9, family="Helvetica", colour="black"),
        axis.text=element_text(size=7, colour="black"))

        
#svg(file=paste0(prefix, "_snp_telcen-ma9.svg"), width=3.25, height=2.25)
#pTelCen.ma9
#dev.off()
pdf(file=file.path(dirout, paste0(prefix, "_meth_telcen-ma9.pdf")), width=3.25, height=2.25)
pTelCen.meth.ma9
dev.off()
png(file=file.path(dirout, paste0(prefix, "_meth_telcen-ma9.png")), width=3.25, height=2.25, res=300, unit="in")
pTelCen.meth.ma9
dev.off()

