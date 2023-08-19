#!/bin/bash

wd="/location/where/the/codes/are/"
src="$wd/DCO_analysis.R"
dircotable="$wd/input"
dirsex="$wd/input/sex-specific"
dirout="$wd"
#mkdir -p $dirout

wt="$dircotable/wt2021_cotable.csv"
recq4="$dircotable/serra17-recq4_cotable.csv"
j3="$dircotable/pJ3_mJ3_cotable.csv"
j3_recq4="$dircotable/mj3-recq4_cotable.csv"

wt_male="$dirsex/wt_male_cotable.csv"
wt_female="$dirsex/wt_female_cotable.csv"
j3_male="$dirsex/pJ3-mJ3-male_cotable.csv"
j3_female="$dirsex/pJ3-mJ3-female_cotable.csv"

#Rscript $src $wt "wt" $dirout &
#Rscript $src $recq4 "recq4" $dirout &
#Rscript $src $j3 "pJ3-mJ3" $dirout &
#Rscript $src $j3_recq4 "pJ3-mJ3_recq4" $dirout &

Rscript $src $wt_male "wt_male" $dirout &
Rscript $src $j3_male "pJ3-mJ3_male" $dirout &
Rscript $src $wt_female "wt_female" $dirout &
Rscript $src $j3_female "pJ3-mJ3_female" $dirout &
wait
echo finished

