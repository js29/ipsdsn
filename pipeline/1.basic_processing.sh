#!/bin/bash
# This file is only meant to be used as an example. Commands like these can be run one at
# a time on the command line, but need to be adapted for a given compute job submission
# system.

STARDIR=STAR

# Convert cram files into fastq
cut -f1 data/metadata.all.txt | submitJobs.py --MEM 4000 --jobname cramToFastq --command "python utils/cramToFastq.py --inputDir cram/ --outputDir fastq/"

# Align reads to the genome using STAR
cut -f1 data/metadata.all.txt | submitJobs.py --MEM 32000 --jobname star_align --ncores 8 --command "python utils/STAR-align.py --outputDir $STARDIR --fastqDir fastq/ --genomeDir annotations/GRCh38/STAR_index_79/ --runThreadN 8"

#Index bam files
cut -f1 data/metadata.all.txt | submitJobs.py --MEM 1000 --jobname index_bams --command  "python utils/index-bams.py --bamdir $STARDIR --insuffix .Aligned.sortedByCoord.out.bam --execute True"

cut -f1 data/metadata.all.txt | submitJobs.py --MEM 1000 --jobname featureCounts --command "python utils/bam2counts.py --sampleDir $STARDIR --gtf annotations/GRCh38/Homo_sapiens.GRCh38.79.gencode_basic.gtf --strand 2 --countsSuffix .gencode_basic.counts.txt --execute True"

# See how many mapped reads per sample:
find $STARDIR -name '*gencode_basic.counts.txt.summary' -print0 | xargs -r0 grep -H 'Assigned'
find $STARDIR -name '*gencode_basic.counts.txt.summary' -print0 | xargs -r0 grep -H 'Unassigned_MultiMapping'
find $STARDIR -name '*gencode_basic.counts.txt.summary' -print0 | xargs -r0 grep -H 'Unassigned_NoFeatures'
find $STARDIR -name '*gencode_basic.counts.txt.summary' -print0 | xargs -r0 grep -H 'Unassigned_Ambiguity'


# Merge together replicates from the same donor
cat data/samples_to_merge.txt | submitJobs.py --MEM 100 --jobname mergeBams --command "python utils/mergeBams.py --indir $STARDIR --outdir $STARDIR --insuffix .Aligned.sortedByCoord.out.bam --outsuffix .merged.bam"

# Get expression counts for merged samples, to use in QTL calling
cut -f1 data/samples_to_merge.txt | submitJobs.py --MEM 1000 --jobname featureCounts.merged --command "python utils/bam2counts.py --sampleDir $STARDIR --bamSuffix .merged.bam --gtf annotations/GRCh38/Homo_sapiens.GRCh38.79.gencode_basic.gtf --strand 2 --countsSuffix .gencode_basic.counts.txt --execute True"

# Sort the bam files
cut -f1 data/samples_to_merge.txt | submitJobs.py --MEM 6000 --jobname sortbams.merged --command  "python utils/bamSortCoord.py --indir $STARDIR --outdir $STARDIR --insuffix .merged.bam --outsuffix .merged.sorted.bam"

#Index bam files
cut -f1 data/samples_to_merge.txt | submitJobs.py --MEM 1000 --jobname index_bams.merged --command  "python utils/index-bams.py --bamdir $STARDIR --insuffix .merged.sorted.bam --execute True"

# Aggregate counts from all samples

# Aggregate counts from samples for QTL calling
submitJobs.py --MEM 1000 -j "combineCounts" -c "Rscript utils/combineCounts.R data/metadata.all.txt data/all_basic_counts.txt"

submitJobs.py --MEM 1000 -j "combineCounts.qtl" -c "Rscript utils/combineCounts.R data/metadata.qtl_samples.txt data/all_basic_counts.qtl.txt"


# This step takes a while... over an hour
submitJobs.py --MEM 2000 -j "processExpressionData" "Rscript utils/processExpressionData.R eqtl/all_basic_counts.txt data/metadata.all.txt data/combined_expression_data"

submitJobs.py --MEM 2000 -j "processExpressionData.qtl" "Rscript utils/processExpressionData.R eqtl/all_basic_counts.qtl.txt data/metadata.all.txt data/combined_expression_data.qtl"

