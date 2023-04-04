library(tidyverse)

#dirin <- "/datasets/data_4/nison/GBS/0_analysis/hcr3-male-female/1_cotables"
#wt_male <- read_csv(file.path(dirin, "wt_male_cotable.csv"), col_names=T)
#input <- wt_male
#wt_female <- read_csv(file.path(dirin, "wt_female_cotable.csv"), col_names=T)
#
#pJ3_mJ3_male <- read_csv(file.path(dirin, "pJ3-mJ3-male_cotable.csv"), col_names=T)
#pJ3_mJ3_female <- read_csv(file.path(dirin, "pJ3-mJ3-female_cotable.csv"), col_names=T)
#setwd("/datasets/data_4/nison/GBS/0_analysis/hcr3-male-female/3_DCO_analysis")


args <- commandArgs(trailingOnly=T)
input <- args[1]
prefix <- args[2]
dirout <- args[3]
chr_size=c("Chr1"=30427671, "Chr2"=19698289, "Chr3"=23459830, "Chr4"=18585056, "Chr5"=26975502)

###
#dirin <- "/datasets/data_4/nison/GBS/0_analysis/hcr3-male-female/1_cotables"
#input <- file.path(dirin, "wt_male_cotable.csv")
#input <- file.path(dirin, "pJ3-mJ3-male_cotable.csv")
#input <- pJ3_mJ3_male
#prefix <- "wt_male"
#dirout <- "/datasets/data_4/nison/hcr3/DCO_analysis"
###

input <- read_csv(input, col_names=T)

# ggplot theme setting
mytheme <- theme_bw() +
		theme(text=element_text(family="helvetica", size=8)) +
		theme(panel.grid = element_blank()) +
		theme(legend.key.size=unit(0.1, "in")) +
		theme(legend.margin=margin(2, 2, 2, 2))
theme_set(mytheme)

#cotable <- wt_male

# function simulateCOS: simulate CO site based on the number of CO in the chr of a library {{{
simulateCOS <- function(dat){
		lib <- dat[1] %>%
				as.numeric()
		chrs <- dat[2]
		nco <- dat[3]
		cos <- round(runif(nco, min=1, max=chr_size[chrs]))
		sim.cotable <- tibble(lib=lib, chrs=chrs, cos=cos)
		return(sim.cotable)
}
#}}}

# function simulateF2: simulate cotable with same CO distribution {{{
simulateF2 <- function(dat){
		dat <- dat %>%
				arrange(lib, chrs)

		cotable.summary <- dat %>%
				group_by(lib, chrs) %>%
				summarise(nco=n())
		f2 <- apply(cotable.summary, 1, simulateCOS) %>%
				bind_rows() %>%
				arrange(lib, chrs, cos)
		return(f2)
}
#}}}

# function filterDCO: filterDCO from cotable {{{
filterDCO <- function(dat){
		#dat <- cotable
		ncos <- dat %>%
				group_by(lib, chrs) %>%
				summarise(nco=n())
		dat.ncos <- left_join(dat, ncos, by=c("lib", "chrs"))
		dat.dco <- filter(dat.ncos, nco >=2)
}
# }}}

# function calcDCOdistance: calculate DCO distance of cotable {{{
calcDCOdistance <- function(dat){
		#dat <- cotable.dco
		libUniq <- unique(dat$lib)
		distance_vector <- c()
		for (i in libUniq){
				#i=1
				byLib <- filter(dat, lib==i)
				chrUniq <- unique(byLib$chrs)
				for(chr in chrUniq){
						byLibChr <- filter(byLib, chrs==chr)
						dist1 <- byLibChr$cos
						dist2 <- c(0, dist1[1:(length(dist1)-1)])
						distance <- dist1-dist2
						distance <- distance[2:length(distance)]
						distance_vector <- c(distance_vector, distance)
				}
		}
		return(distance_vector)
}
#}}}

main <- function(cotable){
		###
#		cotable <- input
		###
		simulation <- simulateF2(cotable)

		obs.dco <- filterDCO(cotable)
		sim.dco <- filterDCO(simulation)

		obs.distance <- calcDCOdistance(obs.dco)
		sim.distance <- calcDCOdistance(sim.dco)

		distance <- list(random=sim.distance, observ=obs.distance)
		return(distance)
}

distance <- main(input)

#sink
sink(file.path(dirout, paste0(prefix, "_dco-distance_summary.txt")))
cat(paste("no CO=", nrow(input), "\n"))
cat(paste("no DCO=", length(distance$observ), "\n"))
lapply(distance, summary)
wilcox.test(distance$observ, distance$random)
sink()

toMb <- function(bp){
		mb <- bp/1e+6
		return(mb)
}

distanceDistPlot <- function(distance, prefix){
		#distance <- distance
		mean.obs <- mean(distance$observ)
		mean.rand <- mean(distance$random)

		distance.tb <- as_tibble(distance) %>%
				pivot_longer(1:2, names_to="type", values_to="distance") %>%
				arrange(type, distance)


		p <- ggplot(distance.tb, aes(x=distance, colour=type)) +
				geom_density() +
				scale_colour_manual(values=c("observ"="red", "random"="blue")) +
				scale_x_continuous(labels=~./1e+6, limits=c(0, 30e+6)) +
				geom_vline(xintercept=mean.obs, colour="red", linetype=2) +
				geom_vline(xintercept=mean.rand, colour="blue", linetype=2) +
				labs(x="Distance (Mb)") +
				ggtitle(prefix)

		svg(file=file.path(dirout, paste0(prefix, "_DCO_distance_distribution.svg")), width=3, height=2.5)
		print(p)
		dev.off()
}

distanceDistPlot(distance, prefix)

# The number of COs per chromosome
coPerChr <- input %>%
		group_by(lib, chrs) %>%
		summarise(nco=n())
pCoPerChr <- ggplot(coPerChr, aes(x=nco)) +
		geom_bar(fill="#575757", colour="black", aes(y=(..prop..)*100)) +
		scale_x_continuous(limits=c(0, 7.5), breaks=0:7, labels=0:7) +
		scale_y_continuous(limits=c(0, 100)) +
		labs(x="Crossovers per Chromatid", y="Percentage") +
		ggtitle(prefix)
svg(file=file.path(dirout, paste0(prefix, "_CO-per-chromatid.svg")), width=2, height=2)
print(pCoPerChr)
dev.off()
