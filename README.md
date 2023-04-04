# HCR3
Codes for double crossover (DCO) analysis used in HCR3 project

# File description
## input: tables of crossover sites analyzed by TIGER pipeline
WT: col48_cotable.csv, col96_cotable.csv, wt2021_cotable.csv
mJ3: pJ3_mJ3_cotable.csv
recq4: recq4_cotable.csv
mj3-recq4_cotable.csv
recq4(Serra17): serra17-recq4_cotable

* Below are sex-specific GBS data *
sex-specific/pJ3-mJ3-female_cotable.csv: female-specific mJ3 CO sites
sex-specific/pJ3-mJ3-male_cotable.csv: male-specific mJ3 CO sites
sex-specific/WT-female_cotable.csv: female-specific mJ3 CO sites
sex-specific/WT-male_cotable.csv: male-specific mJ3 CO sites

## sources
DCO_analysis.R: R script for double crossover analysis. Take two crossover tables as input
run_DCO_analysis.sh: Automating DCO_analysis of various genotypes
