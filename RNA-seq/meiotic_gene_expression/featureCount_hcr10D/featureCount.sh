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

dirIn="/datasets/data_4/nison/NGS_libraries/RNAseq_hcr10D-series_bud-seedling/01_STAR"

Col_s_1=$dirIn/Col_s_1_Aligned.sortedByCoord.out.bam
Col_s_2=$dirIn/Col_s_2_Aligned.sortedByCoord.out.bam
Col_s_3=$dirIn/Col_s_3_Aligned.sortedByCoord.out.bam

Col_b_1=$dirIn/Col_1_Aligned.sortedByCoord.out.bam
Col_b_2=$dirIn/Col_2_Aligned.sortedByCoord.out.bam
Col_b_3=$dirIn/Col_3_Aligned.sortedByCoord.out.bam

dirOut="/datasets/data_4/nison/hcr3/HCR3_public_release/RNA-seq/meiotic_gene_expression/featureCount_hcr10D"
mkdir -p $dirOut
outfile="Col_s_b_from_hcr10D" # outfile prefix
#featureType is either 0 or 1 (default:0). If 1 is given to the 3rd argument, 'TAIR10_GFF3_all_features.gff' will be used as gff 
gtf='/home/nison/work/refgenome/araport11/Araport11_GTF_genes_transposons.Mar92021.gtf'


featureCounts -T 20 -s 2 -g gene_id -p -t exon \
	-a $gtf \
    -o $dirOut/${outfile}.featureCounts.txt \
    $Col_s_1 $Col_s_2 $Col_s_3 \
    $Col_b_1 $Col_b_2 $Col_b_3 \
	1> $dirOut/${outfile}.featureCounts.log \
	2> $dirOut/${outfile}.featureCounts.log2

#	$dirIn/*_Aligned.sortedByCoord.out.bam \

# -T : number of the threads. The value should be between 1 and 32
# -s : 0, unstranded; 1, stranded; 2, reversely stranded
# -g : attribute type used to group features into meta-features. 'gene_id" by default. This attribute type is usually the gene identifier.a This argument is useful for the meta-feature level summarization.
# -t : Specify the feature type. Only rows which have the matched feature type in the provided GTF annotation file will be included for read counting. 'exon' by default.
