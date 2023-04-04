# HCR3
Codes for double crossover (DCO) analysis used in HCR3 project

# File description
## input: tables of crossover sites analyzed by TIGER pipeline
WT: col48_cotable.csv, col96_cotable.csv, wt2021_cotable.csv<br />
mJ3: pJ3_mJ3_cotable.csv<br />
recq4: recq4_cotable.csv<br />
mj3-recq4_cotable.csv<br />
recq4(Serra17): serra17-recq4_cotable<br />

* Below are sex-specific GBS data <br />
sex-specific/pJ3-mJ3-female_cotable.csv: female-specific mJ3 CO sites<br />
sex-specific/pJ3-mJ3-male_cotable.csv: male-specific mJ3 CO sites<br />
sex-specific/WT-female_cotable.csv: female-specific mJ3 CO sites<br />
sex-specific/WT-male_cotable.csv: male-specific mJ3 CO sites<br />

## sources
DCO_analysis.R: R script for double crossover analysis. Take two crossover tables as input<br />
run_DCO_analysis.sh: Automating DCO_analysis of various genotypes<br />
