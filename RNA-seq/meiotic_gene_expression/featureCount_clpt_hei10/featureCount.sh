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

dirIn="/datasets/data_4/nison/clpt/RNA-seq/2023_clpt-socts/01_STAR"

col_s_1=$dirIn/Col_1_Aligned.sortedByCoord.out.bam
col_s_2=$dirIn/Col_2_Aligned.sortedByCoord.out.bam
col_s_3=$dirIn/Col_3_Aligned.sortedByCoord.out.bam


dirOut="/datasets/data_4/nison/hcr3/HCR3_public_release/RNA-seq/meiotic_gene_expression/featureCount_clpt_hei10"
mkdir -p $dirOut
outfile="Col_hei10_from_clpt-series" # outfile prefix
#featureType is either 0 or 1 (default:0). If 1 is given to the 3rd argument, 'TAIR10_GFF3_all_features.gff' will be used as gff 
gtf='/datasets/data_4/nison/hcr3/HCR3_public_release/RNA-seq/meiotic_gene_expression/featureCount_hei10/HEI10_exons_selected_to_be_counted.gtf'


featureCounts -T 20 -s 2 -g gene_id -p -t exon \
	-a $gtf \
    -o $dirOut/${outfile}.featureCounts.txt \
    $col_s_1 $col_s_2 $col_s_3 \
	1> $dirOut/${outfile}.featureCounts.log \
	2> $dirOut/${outfile}.featureCounts.log2

#	$dirIn/*_Aligned.sortedByCoord.out.bam \

# -T : number of the threads. The value should be between 1 and 32
# -s : 0, unstranded; 1, stranded; 2, reversely stranded
# -g : attribute type used to group features into meta-features. 'gene_id" by default. This attribute type is usually the gene identifier.a This argument is useful for the meta-feature level summarization.
# -t : Specify the feature type. Only rows which have the matched feature type in the provided GTF annotation file will be included for read counting. 'exon' by default.
