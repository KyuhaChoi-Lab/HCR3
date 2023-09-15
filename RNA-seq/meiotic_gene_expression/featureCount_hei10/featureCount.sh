#!/bin/bash
# Maintainer : nison
# template: STAR2featureCounts_v2.sh
# last modified : 2021. 08. 12. 
# Description :
# This script takes bam files as a input, then count the number of reads on each of genomic features.# Genomic features are either "gene" or "all feature" in GFF3 file specified.

# You can specify the type of feature to use by giving the 4th argument (default=0). If 1 is given, all features in GFF3 file will be used including TE.

# 1. Use Araport11.gtf for annotation. gtf format is better for meta-feature counting than gff format.
# 2. -t : exon --> count reads hit on exon
# 3. -g : gene --> group features into gene (meta-feature)

###### INPUT FILES ######
# Col-s, Col-b: hcr11D series (2023) --> set A
# Col-b, hcr3: hcr3 series (2021) --> Set B
# Col-b, j3-3, j2-2: hcr3-j3-j2 series (2023) --> set C
# 

dirIn_A="/datasets/data_4/nison/NGS_libraries/RNAseq_hcr11D-series_bud-seedling/01_STAR"
dirIn_B="/datasets/data_4/nison/hcr3/RNA-seq/hcr3-series_2021/01_STAR"
dirIn_C="/datasets/data_4/nison/NGS_libraries/RNAseq_j2-j3_meiotic-bud/01_STAR"
dirIn_D="/datasets/data_4/nison/hcr2/RNA-seq/02_STAR"

Col_s_A_1=$dirIn_A/Col_s_1_Aligned.sortedByCoord.out.bam
Col_s_A_2=$dirIn_A/Col_s_2_Aligned.sortedByCoord.out.bam
Col_s_A_3=$dirIn_A/Col_s_3_Aligned.sortedByCoord.out.bam

Col_b_A_1=$dirIn_A/Col_b_1_Aligned.sortedByCoord.out.bam
Col_b_A_2=$dirIn_A/Col_b_2_Aligned.sortedByCoord.out.bam
Col_b_A_3=$dirIn_A/Col_b_3_Aligned.sortedByCoord.out.bam

Col_b_B_1=$dirIn_B/col-bud-1_Aligned.sortedByCoord.out.bam
Col_b_B_2=$dirIn_B/col-bud-2_Aligned.sortedByCoord.out.bam
Col_b_B_3=$dirIn_B/col-bud-3_Aligned.sortedByCoord.out.bam

hcr3_b_B_1=$dirIn_B/hcr3-bud-1_Aligned.sortedByCoord.out.bam
hcr3_b_B_2=$dirIn_B/hcr3-bud-2_Aligned.sortedByCoord.out.bam
hcr3_b_B_3=$dirIn_B/hcr3-bud-3_Aligned.sortedByCoord.out.bam

Col_b_C_1=$dirIn_C/Col_bud_1_Aligned.sortedByCoord.out.bam
Col_b_C_2=$dirIn_C/Col_bud_2_Aligned.sortedByCoord.out.bam
Col_b_C_3=$dirIn_C/Col_bud_3_Aligned.sortedByCoord.out.bam

j3_b_C_1=$dirIn_C/j3-3_bud_1_Aligned.sortedByCoord.out.bam
j3_b_C_2=$dirIn_C/j3-3_bud_2_Aligned.sortedByCoord.out.bam
j3_b_C_3=$dirIn_C/j3-3_bud_3_Aligned.sortedByCoord.out.bam

j2_b_C_1=$dirIn_C/j2-2_bud_1_Aligned.sortedByCoord.out.bam
j2_b_C_2=$dirIn_C/j2-2_bud_2_Aligned.sortedByCoord.out.bam
j2_b_C_3=$dirIn_C/j2-2_bud_3_Aligned.sortedByCoord.out.bam

#col_s_D_1=$dirIn_D/Col-0-seedling1_aligned.sortedByCoord.out.bam
#col_s_D_2=$dirIn_D/Col-0-seedling2_aligned.sortedByCoord.out.bam
#col_s_D_3=$dirIn_D/Col-0-seedling3_aligned.sortedByCoord.out.bam
#col_s_D_4=$dirIn_D/Col-0-seedling4_aligned.sortedByCoord.out.bam
#
#col_b_D_1=$dirIn_D/Col-0-bud-19_aligned.sortedByCoord.out.bam
#col_b_D_2=$dirIn_D/Col-0-bud-20_aligned.sortedByCoord.out.bam
#col_b_D_3=$dirIn_D/Col-0-bud-21_aligned.sortedByCoord.out.bam
#col_b_D_4=$dirIn_D/Col-0-bud-22_aligned.sortedByCoord.out.bam
#
#hcr2_b_1=$dirIn_D/hcr2-bud-23_Aligned.sortedByCoord.out.bam
#hcr2_b_2=$dirIn_D/hcr2-bud-24_Aligned.sortedByCoord.out.bam
#hcr2_b_3=$dirIn_D/hcr2-bud-25_Aligned.sortedByCoord.out.bam
#hcr2_b_4=$dirIn_D/hcr2-bud-26_Aligned.sortedByCoord.out.bam

dirOut="/datasets/data_4/nison/hcr3/HCR3_public_release/RNA-seq/meiotic_gene_expression/featureCount_hei10"
mkdir -p $dirOut
outfile="hcr3-paper_total-RNAseq_hei10" # outfile prefix
#featureType is either 0 or 1 (default:0). If 1 is given to the 3rd argument, 'TAIR10_GFF3_all_features.gff' will be used as gff 
gtf='/datasets/data_4/nison/hcr3/HCR3_public_release/RNA-seq/meiotic_gene_expression/featureCount_hei10/HEI10_exons_selected_to_be_counted.gtf'


featureCounts -T 20 -s 2 -g gene_id -p -t exon \
	-a $gtf \
    -o $dirOut/${outfile}.featureCounts.txt \
    $Col_s_A_1 $Col_s_A_2 $Col_s_A_3 \
    $Col_b_A_1 $Col_b_A_2 $Col_b_A_3 \
    $Col_b_B_1 $Col_b_B_2 $Col_b_B_3 \
    $hcr3_b_B_1 $hcr3_b_B_2 $hcr3_b_B_3 \
    $Col_b_C_1 $Col_b_C_2 $Col_b_C_3 \
    $j3_b_C_1 $j3_b_C_2 $j3_b_C_3 \
    $j2_b_C_1 $j2_b_C_2 $j2_b_C_3 \
	1> $dirOut/${outfile}.featureCounts.log \
	2> $dirOut/${outfile}.featureCounts.log2

#	$dirIn/*_Aligned.sortedByCoord.out.bam \

# -T : number of the threads. The value should be between 1 and 32
# -s : 0, unstranded; 1, stranded; 2, reversely stranded
# -g : attribute type used to group features into meta-features. 'gene_id" by default. This attribute type is usually the gene identifier.a This argument is useful for the meta-feature level summarization.
# -t : Specify the feature type. Only rows which have the matched feature type in the provided GTF annotation file will be included for read counting. 'exon' by default.
