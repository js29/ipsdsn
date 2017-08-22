#!/bin/bash
# This file is only meant to be used as an example. Similar commands can be run one at
# a time on the command line.
STARDIR=STAR
VCFBASE=imputed.97_samples
#mkdir rasqual
RASQUALDIR=/full/path/here/ipsdsn/rasqual

# Tried using Rasqual's ASVCF to create the VCF file with allele-specific counts, but
# ended up using a different method.
#export RASQUALDIR=/nfs/users/nfs_j/js29/software/rasqual/
#cd genotypes/GRCh38/imputed_20151005/
#submitJobs.py --MEM 5000 -j "createASVCF" -c "bash ~/software/rasqual/src/ASVCF/createASVCF.sh rasqual.bam_list.txt $VCFBASE.snps_only.INFO_08.named.vcf.gz"

# Add read group to the BAM files to make them work with GATK
# Pass the sample ID and HIPSCI ID
cut -f 1,4 data/metadata.qtl_samples.txt | sed '1d' | submitJobs.py --MEM 1000 --jobname bamAddRG --command "python utils/bamAddRG.py --indir $STARDIR --outdir $STARDIR --insuffix .Aligned.sortedByCoord.out.bam --outsuffix .Aligned.sortedByCoord.RG.bam --execute True"

cut -f 1 data/metadata.qtl_samples.txt | sed '1d' | submitJobs.py --MEM 1000 --jobname index_bams --command "python utils/index-bams.py --bamdir $STARDIR --insuffix .Aligned.sortedByCoord.RG.bam --execute True"

# We use the snps_only VCF to get allele-specific counts for Rasqual, because it can't
# use indels as feature snps. But we used both SNPs and indels for QTL calling
# i.e. in the VCF that we pass to Rasqual as input.
VCF_ASE=$SN/genotypes/GRCh38/imputed_20151005/imputed.98_samples.snps_only.INFO_08.named.vcf.gz
VCF_QTL_IN=$SN/genotypes/GRCh38/imputed_20151005/imputed.98_samples.snps_indels.INFO_08.named.vcf.gz
GRCh38_FASTA=annotation/Homo_sapiens.GRCh38.dna.primary_assembly.fa

# Use GATK's ASEReadCounter to count allele-specific expression
cut -f1 data/metadata.qtl_samples.txt | sed '1d' | submitJobs.py --MEM 4000 -j bamCountASE -c "python utils/bamCountASE.py --indir $STARDIR --outdir $STARDIR --insuffix .Aligned.sortedByCoord.RG.bam --reference $GRCh38_FASTA --sites $VCF_ASE --Xmx 2500m --execute True"

# Merge ASE read counts per condition into one file:
paste <(cut -f 1 data/metadata.qtl_samples.txt) <(cut -f 1 data/metadata.qtl_samples.txt) > rasqual/sensoryneurons.sample_map.txt
submitJobs.py --MEM 6000 -j mergeASEcounts -q yesterday -c "python utils/mergeASECounts.py --sample_list rasqual/sensoryneurons.sample_map.txt --indir $STARDIR --suffix .ASEcounts | gzip > rasqual/$VCFBASE.RNAseq_ASE_counts.txt.gz"


cd $RASQUALDIR
mkdir input
mkdir output
#Construct genotype list
cut -f 4 ../data/metadata.qtl_samples.txt > sensoryneurons.genotype_list.txt 

#Extract only relevant genotypes from the vcf file (in the correct order)
bcftools view -O z -S sensoryneurons.genotype_list.txt $VCF_QTL_IN > $VCFBASE.snps_indels.INFO_08.selected_samples.vcf.gz

#Add ASE counts to the VCF file
VCF_QTL=$VCFBASE.snps_indels.INFO_08.MAF_0.05.RASQUAL.vcf.gz
cut -f 1,4 ../data/metadata.qtl_samples.txt > sensoryneurons.sample_genotype_map.txt
submitJobs.py --MEM 12000 -j vcfAddASE -q yesterday -c "python ../utils/vcfAddASE.py --ASEcounts $VCFBASE.RNAseq_ASE_counts.txt.gz --VCFfile $VCFBASE.snps_indels.INFO_08.selected_samples.vcf.gz --ASESampleGenotypeMap sensoryneurons.sample_genotype_map.txt | bgzip > input/$VCF_QTL"

#Index the VCF files using tabix
jsub.sh "tabix" "tabix -p vcf input/$VCF_QTL"

# Extract SNP coords from a vcf file
zgrep -v "#" input/$VCF_QTL | cut -f 1,2,3 > $VCFBASE.snp_coords.txt

# Check that there are no duplicate SNP coordinates, so we don't overcount SNPs
# (these should give the same number)
wc -l rasqual/$VCFBASE.snp_coords.txt
cat rasqual/$VCFBASE.snp_coords.txt | uniq | wc -l

# File ../annotations/Homo_sapiens.GRCh38.79.gene_exon_start_end.filtered_genes.txt has
# the format:
# gene_id	chromosome_name	strand	exon_starts	exon_ends
# ENSG00000044012	1	1	42153421,42154680,42155535	42153540,42154866,42155824
# ...

#Count SnpPerGene
submitJobs.py --MEM 5000 -j countsSnpsPerGene.100k -q yesterday -c "Rscript ../utils/countsSnpsPerGene.R -s $VCFBASE.snp_coords.txt -w 100000 -e $SN/annotations/Homo_sapiens.GRCh38.79.gene_exon_start_end.filtered_genes.txt -c input/$VCFBASE.snps_per_gene.100k.txt"
submitJobs.py --MEM 5000 -j countsSnpsPerGene.500k -q yesterday -c "Rscript ../utils/countsSnpsPerGene.R -s $VCFBASE.snp_coords.txt -w 500000 -e $SN/annotations/Homo_sapiens.GRCh38.79.gene_exon_start_end.filtered_genes.txt -c input/$VCFBASE.snps_per_gene.500k.txt"

# Get counts for QTL samples, as was done for FastQTL
cut -f 2 --complement data/all_basic_counts.qtl.txt > $BASE.basic_counts.txt

# Run Rasqual's scripts to get the "offset" for each sample
cd ~/software/rasqual
R --vanilla --quiet --args $RASQUALDIR/$VCFBASE.basic_counts.txt < R/makeOffset.R > $RASQUALDIR/makeOffset.log.txt
mv data/your.K.txt $RASQUALDIR/$VCFBASE.sample_offsets.txt

# Use Rasqual's script to get covariates (PCs)
R --vanilla --quiet --args $RASQUALDIR/$VCFBASE.basic_counts.txt $RASQUALDIR/$VCFBASE.sample_offsets.txt < R/makeCovariates.R > $RASQUALDIR/makeCovariates.log.txt
mv data/your.X.txt $RASQUALDIR/$VCFBASE.covariates.txt
 
R --vanilla --quiet --args $RASQUALDIR/$VCFBASE.basic_counts.txt $RASQUALDIR/$VCFBASE.sample_offsets.txt $RASQUALDIR/$VCFBASE.covariates.txt < R/txt2bin.R > $RASQUALDIR/txt2bin.log.txt
# Copy .bin files generated in rasqual/data directory to $SN/rasqual/input
cd $RASQUALDIR
mv $VCFBASE.basic_counts.bin input/
mv $VCFBASE.sample_offsets.bin input/
mv $VCFBASE.covariates.bin input/
cut -f 1 $VCFBASE.basic_counts.txt > input/gene_id_list.txt


###############################################################################
######### Run RASQUAL
VCFBASE=imputed.97_samples
cd rasqual
VCF_QTL=$VCFBASE.snps_indels.INFO_08.MAF_0.05.RASQUAL.vcf.gz

# Prepare job batches

# Split genes into those where we have feature SNPs (Rasqual can use AS model) and those
# with no feature SNPs (Rasqual can only use population model) 
DIST=500k
cat input/$VCFBASE.snps_per_gene.$DIST.txt | sed '1d' | awk '$2 !~ /[0-9]+/' | perl -ane 'chomp; print join("\t", @F, 0)."\n"' > input/non_autosome_genes.$DIST.txt
cat input/$VCFBASE.snps_per_gene.$DIST.txt | sed '1d' | awk '($2 ~ /[0-9]+/ && $8 == 0)' | perl -ane 'chomp; print join("\t", @F, $F[7]*$F[8])."\n"' > input/autosome_genes.no_fsnps.$DIST.txt
cat input/$VCFBASE.snps_per_gene.$DIST.txt | sed '1d' | awk '($2 ~ /[0-9]+/ && $8 > 0)' | perl -ane 'chomp; print join("\t", @F, $F[7]*$F[8])."\n"' | sort -nk10,10 > input/autosome_genes.with_fsnps.$DIST.txt
cat input/non_autosome_genes.$DIST.txt | perl ../utils/rasqual.makeBatches.pl --size 50000 > input/non_autosome_genes.$DIST.batches.txt
cat input/autosome_genes.with_fsnps.$DIST.txt | perl ../utils/rasqual.makeBatches.pl --size 50000 --excludelong 50000 > input/autosome_genes.with_fsnps.$DIST.normalbatches.txt
cat input/autosome_genes.with_fsnps.$DIST.txt | perl ../utils/rasqual.makeBatches.pl --size 50000 --onlylong 50000 > input/autosome_genes.with_fsnps.$DIST.longbatches.txt
cat input/autosome_genes.no_fsnps.$DIST.txt | perl ../utils/rasqual.makeBatches.pl --size 50000 > input/autosome_genes.no_fsnps.$DIST.batches.txt


# Run RASQUAL with covariates
# We also use the option --no-posterior-update. This prevents RASQUAL from adjusting the
# likelihood based on its estimated corrections to sample genotypes. I found that even
# with permuted data, quite a few genes would have extreme p values (e.g. 1e-100) when
# this option was not set.
D=$RASQUALDIR
NUMSAMPLES=97
DIST=500k
mkdir output.$DIST
cat input/non_autosome_genes.$DIST.batches.txt | python ../utils/submitJobArray.py --MEM 200 --jobname runRasqual.$DIST.non_autosome_genes --arraymaxjobs 1000 --farmout $D/FarmOut.$DIST --queue normal --command "python ../utils/runRasqual.py --rasqualargs '--population-only --no-posterior-update' --outprefix $D/output.$DIST/rasqual.$DIST.non_autosome_genes --readCounts $D/input/$VCFBASE.basic_counts.bin --offsets $D/input/$VCFBASE.sample_offsets.bin --covariates $D/input/$VCFBASE.covariates.bin --n $NUMSAMPLES --geneids $D/input/gene_id_list.txt --vcf $D/input/$VCF_QTL --geneMetadata $D/input/$VCFBASE.snps_per_gene.$DIST.txt --execute True"
cat input/autosome_genes.no_fsnps.$DIST.batches.txt | python ../utils/submitJobArray.py --MEM 200 --jobname runRasqual.$DIST.no_fsnps --arraymaxjobs 1000 --farmout $D/FarmOut.$DIST --queue normal --command "python ../utils/runRasqual.py --rasqualargs '--n-threads 2 --no-posterior-update' --outprefix $D/output.$DIST/rasqual.$DIST.no_fsnps --readCounts $D/input/$VCFBASE.basic_counts.bin --offsets $D/input/$VCFBASE.sample_offsets.bin --covariates $D/input/$VCFBASE.covariates.bin --n $NUMSAMPLES --geneids $D/input/gene_id_list.txt --vcf $D/input/$VCF_QTL --geneMetadata $D/input/$VCFBASE.snps_per_gene.$DIST.txt --execute True"
cat input/autosome_genes.with_fsnps.$DIST.normalbatches.txt | python ../utils/submitJobArray.py --MEM 200 --jobname runRasqual.$DIST.with_fsnps --arraymaxjobs 1000 --farmout $D/FarmOut.$DIST --queue normal --command "python ../utils/runRasqual.py --rasqualargs '--n-threads 2 --no-posterior-update' --outprefix $D/output.$DIST/rasqual.$DIST.with_fsnps --readCounts $D/input/$VCFBASE.basic_counts.bin --offsets $D/input/$VCFBASE.sample_offsets.bin --covariates $D/input/$VCFBASE.covariates.bin --n $NUMSAMPLES --geneids $D/input/gene_id_list.txt --vcf $D/input/$VCF_QTL --geneMetadata $D/input/$VCFBASE.snps_per_gene.$DIST.txt --execute True"
cat input/autosome_genes.with_fsnps.$DIST.longbatches.txt | python ../utils/submitJobArray.py --MEM 200 --jobname runRasqual.$DIST.with_fsnps.long --arraymaxjobs 1000 --farmout $D/FarmOut.$DIST --queue long --command "python ../utils/runRasqual.py --rasqualargs '--n-threads 8 --no-posterior-update' --outprefix $D/output.$DIST/rasqual.$DIST.with_fsnps.long --readCounts $D/input/$VCFBASE.basic_counts.bin --offsets $D/input/$VCFBASE.sample_offsets.bin --covariates $D/input/$VCFBASE.covariates.bin --n $NUMSAMPLES --geneids $D/input/gene_id_list.txt --vcf $D/input/$VCF_QTL --geneMetadata $D/input/$VCFBASE.snps_per_gene.$DIST.txt --execute True"

grep TERM FarmOut.$DIST/*.txt
grep -i "segmentation fault" FarmOut.500k/*.txt > FarmOut.500k.segfaults.txt

submitJobs.py --MEM 200 -j rasqual.combineGenes -c "cat output.$DIST/rasqual.$DIST.*.txt | grep -v SKIPPED | bgzip > output/rasqual.$DIST.all.txt.gz"

# Get the lead SNP per gene
submitJobs.py --MEM 200 -j getLeadSNPs.$DIST -q yesterday -c "gzip -cd output/rasqual.$DIST.all.txt.gz | grep -v \"SKIPPED\" | perl ../utils/getLeadSnps.pl -f stdin --genecol 1 --statcol 11 > output/rasqual.$DIST.all.lead.txt"

# Filter Rasqual results based on raw p value (chisq value >= 6.63)
submitJobs.py --MEM 100 -j filterRasqualResults -q yesterday -c "gzip -cd output/rasqual.$DIST.all.txt.gz | awk '\\\$11 >= 6.634896' | gzip > output/rasqual.$DIST.pthreshold.0.01.txt.gz"
# This gave us a file with all SNPs having p < 0.01. Some genes will not be included in the
# file at all, but this is a problem for downstream analysis, since all genes were tested.
# So we want to add in at least the one best SNP per gene, even for those genes with lead
# SNP p value > 0.01.
# Write out the set of genes present in the filtered results file.
zcat output/rasqual.$DIST.pthreshold.0.01.txt.gz | cut -f 1 | sort | uniq >> output/rasqual.$DIST.pthreshold.0.01.uniquegenes.txt
cat output/rasqual.$DIST.all.lead.txt | awk '$11 < 6.634896' > output/rasqual.$DIST.lead.abovethreshold.txt
wc -l output/rasqual.$DIST.pthreshold.0.01.uniquegenes.txt
wc -l output/rasqual.$DIST.lead.abovethreshold.txt
# The sum of number of lines in the above two files should equal the total number of genes
# for which we have results (32682).
(zcat output/rasqual.$DIST.pthreshold.0.01.txt.gz; cat output/rasqual.$DIST.lead.abovethreshold.txt) | gzip > output/rasqual.$DIST.pthreshold.0.01.allgenes.txt.gz

echo -e "gene\tsnpid\tchr\tpos\tref\talt\taf\thwechisq\timputation_quality\tchisq\teffect_size\tmap_error_rate\tref_mapping_bias\toverdispersion\tfsnps\ttestedsnps\tconvergence\tr2_genotype_fsnps\tr2_genotype_rsnps" | gzip > output/rasqual.$DIST.pthreshold.0.01.allgenes.cut.txt.gz
zcat output/rasqual.$DIST.pthreshold.0.01.allgenes.txt.gz | cut -f 1-9,11-15,17,18,23-25 | gzip >> output/rasqual.$DIST.pthreshold.0.01.allgenes.cut.txt.gz

cat input/autosome_genes.with_fsnps.500k.longbatches.txt input/autosome_genes.with_fsnps.500k.normalbatches.txt input/autosome_genes.no_fsnps.500k.batches.txt input/non_autosome_genes.500k.batches.txt | tr ',' '\n' | wc -l
# Count how many genes were "skipped"
cat output.500k/rasqual.500k.*.txt | grep SKIPPED | wc -l

tar -zcvf FarmOut.500k.tar.gz FarmOut.500k
tar -zcvf output.500k.tar.gz output.500k


###############################################################################
######### Run RASQUAL with permutations to get null distribution

SN=/lustre/scratch109/sanger/dg13/share/jeremy/sensoryneurons
VCFBASE=imputed.97_samples
VCF_QTL=$VCFBASE.snps_indels.INFO_08.MAF_0.05.RASQUAL.vcf.gz
D=$RASQUALDIR
NUMSAMPLES=97
cd $RASQUALDIR

DIST=500k
cat input/non_autosome_genes.$DIST.batches.txt | python ../utils/submitJobArray.py --MEM 200 -j runRasqual.$DIST.non_autosome_genes --arraymaxjobs 1000 -o $D/FarmOut.$DIST.perm.nopu -q normal -c "python ../utils/runRasqual.py --rasqualargs '--population-only --random-permutation --no-posterior-update' --outprefix $D/output.$DIST.perm.nopu/rasqual.$DIST.non_autosome_genes --readCounts $D/input/$VCFBASE.basic_counts.bin --offsets $D/input/$VCFBASE.sample_offsets.bin --covariates $D/input/$VCFBASE.covariates.bin --n $NUMSAMPLES --geneids $D/input/gene_id_list.txt --vcf $D/input/$VCF_QTL --geneMetadata $D/input/$VCFBASE.snps_per_gene.$DIST.txt --execute True"
cat input/autosome_genes.no_fsnps.$DIST.batches.txt | python ../utils/submitJobArray.py --MEM 200 -j runRasqual.$DIST.no_fsnps --arraymaxjobs 1000 -o $D/FarmOut.$DIST.perm.nopu -q normal -n 2 -c "python ../utils/runRasqual.py --rasqualargs '--n-threads 2 --random-permutation --no-posterior-update' --outprefix $D/output.$DIST.perm.nopu/rasqual.$DIST.no_fsnps --readCounts $D/input/$VCFBASE.basic_counts.bin --offsets $D/input/$VCFBASE.sample_offsets.bin --covariates $D/input/$VCFBASE.covariates.bin --n $NUMSAMPLES --geneids $D/input/gene_id_list.txt --vcf $D/input/$VCF_QTL --geneMetadata $D/input/$VCFBASE.snps_per_gene.$DIST.txt --execute True"
cat input/autosome_genes.with_fsnps.$DIST.normalbatches.txt | python ../utils/submitJobArray.py --MEM 200 -j runRasqual.$DIST.with_fsnps --arraymaxjobs 1000 -o $D/FarmOut.$DIST.perm.nopu -q normal -n 2 -c "python ../utils/runRasqual.py --rasqualargs '--n-threads 2 --random-permutation --no-posterior-update' --outprefix $D/output.$DIST.perm.nopu/rasqual.$DIST.with_fsnps --readCounts $D/input/$VCFBASE.basic_counts.bin --offsets $D/input/$VCFBASE.sample_offsets.bin --covariates $D/input/$VCFBASE.covariates.bin --n $NUMSAMPLES --geneids $D/input/gene_id_list.txt --vcf $D/input/$VCF_QTL --geneMetadata $D/input/$VCFBASE.snps_per_gene.$DIST.txt --execute True"
cat input/autosome_genes.with_fsnps.$DIST.longbatches.txt | python ../utils/submitJobArray.py --MEM 200 -j runRasqual.$DIST.with_fsnps.long --arraymaxjobs 1000 -o $D/FarmOut.$DIST.perm.nopu -q long -n 8 -c "python ../utils/runRasqual.py --rasqualargs '--n-threads 8 --random-permutation --no-posterior-update' --outprefix $D/output.$DIST.perm.nopu/rasqual.$DIST.with_fsnps.long --readCounts $D/input/$VCFBASE.basic_counts.bin --offsets $D/input/$VCFBASE.sample_offsets.bin --covariates $D/input/$VCFBASE.covariates.bin --n $NUMSAMPLES --geneids $D/input/gene_id_list.txt --vcf $D/input/$VCF_QTL --geneMetadata $D/input/$VCFBASE.snps_per_gene.$DIST.txt --execute True"

grep TERM FarmOut.$DIST.perm.nopu/*.txt
grep -i "segmentation fault" FarmOut.$DIST.perm.nopu/*.txt > FarmOut.$DIST.perm.nopu.segfaults.txt

submitJobs.py --MEM 200 -j combineGenes.$DIST.perm.nopu -q yesterday -c "cat output.$DIST.perm.nopu/rasqual.$DIST.*.txt | grep -v SKIPPED | bgzip > output/rasqual.$DIST.perm.nopu.all.txt.gz"
submitJobs.py --MEM 200 -j getLeadSNPs.$DIST.perm.nopu -q yesterday -c "gzip -cd output/rasqual.$DIST.perm.nopu.all.txt.gz | grep -v SKIPPED | perl ../utils/getLeadSnps.pl -f stdin --genecol 1 --statcol 11 > output/rasqual.$DIST.perm.nopu.all.lead.txt"


###############################################################################
######### Run eigenMT on Rasqual output

##### Prepare eigenMT input
mkdir input/eigenMT
submitJobs.py --MEM 1000 -j rasqualToEigenMT -q yesterday -c "python ../utils/rasqualToEigenMT.py --rasqualOut output/rasqual.500k.all.txt.gz > input/eigenMT/rasqual.500k.all.eigenMT_input.txt"

# Get genotypes split up by chromosome
submitJobs.py --MEM 60000 -j getGenotypesByChr -q hugemem -c "Rscript ../utils/getGenotypesByChr.R ../genotypes/imputed.97_samples.snps_indels.INFO_08.selected_samples.rds input/eigenMT/$VCFBASE.genotypes"
cat input/chromosome_list.txt | submitJobs.py --MEM 500 -j splitVCFByChr -q normal -c "python ../utils/vcfSplitByChromosome.py --vcf $D/input/$VCF_QTL --outdir input/vcf --outPrefix imputed.97_samples.snps_indels. --execute True"

# Convert chromosome VCF files to the more compact GDS format
Rscript ~/src/utils/vcf/vcfToGds.R --vcf-directory input/vcf --chr-list input/chromosome_list.txt
# Then I had to run the function eigenMTExportGenotypesByChr from qtl_eigenMT.R
# to get text file versions of the genotype matrices (GDS files)


echo -e "gene_id\tchrom_probe\ts1\ts2" > input/eigenMT/gene.positions.txt
zcat ../fastqtl/input/imputed.97_samples.fastqtl.phenotypes.allgenes.bed.gz | cut -f 1-4 | sed '1d' | awk 'BEGIN{OFS="\t"}{print $4,$1,$2,$3}' >> input/eigenMT/gene.positions.txt

CHR=19
submitJobs.py --MEM 2000 -j eigenMT.chr$CHR -n 7 -q yesterday -c "python ~/software/eigenMT/eigenMT_fixed.py --CHROM $CHR --QTL input/eigenMT/rasqual.500k.all.eigenMT_input.txt --GEN input/eigenMT/chr_$CHR.genotypes.txt --GENPOS input/eigenMT/chr_$CHR.snp_positions.txt --PHEPOS input/eigenMT/gene.positions.txt --OUT output/eigenMT.output.$CHR.txt --cis_dist 1e5"

CHRS=( 6 )
WINDOW=200
for CHR in "${CHRS[@]}"; do
    submitJobs.py --MEM 2000 -j eigenMT.500k.chr$CHR -n 15 -q normal -c "python ~/software/eigenMT/eigenMT_fixed.py --CHROM $CHR --QTL input/eigenMT/rasqual.500k.all.eigenMT_input.txt --GEN input/eigenMT/chr_$CHR.genotypes.txt --GENPOS input/eigenMT/chr_$CHR.snp_positions.txt --PHEPOS input/eigenMT/gene.positions.txt --OUT output/eigenMT.output.500k.win_$WINDOW.$CHR.txt --cis_dist 1e7 --window $WINDOW"
done
CHRS=( 1 2 3 4 5 6 7 8 9 10 11 12 13 14 15 16 17 18 19 20 21 22 X )
WINDOW=200
for CHR in "${CHRS[@]}"; do
    submitJobs.py --MEM 2000 -j eigenMT.500k.chr$CHR -n 7 -q normal -c "python ~/software/eigenMT/eigenMT_fixed.py --CHROM $CHR --QTL input/eigenMT/rasqual.500k.all.eigenMT_input.txt --GEN input/eigenMT/chr_$CHR.genotypes.txt --GENPOS input/eigenMT/chr_$CHR.snp_positions.txt --PHEPOS input/eigenMT/gene.positions.txt --OUT output/eigenMT.output.500k.win_$WINDOW.$CHR.txt --cis_dist 1e7 --window $WINDOW"
done
head -n 1 output/eigenMT.output.500k.win_$WINDOW.1.txt > output/eigenMT.output.500k.win_$WINDOW.allchrs.txt
CHRS=( 1 2 3 4 5 6 7 8 9 10 11 12 13 14 15 16 17 18 19 20 21 22 X )
for CHR in "${CHRS[@]}"; do
    sed '1d' output/eigenMT.output.500k.win_$WINDOW.$CHR.txt >> output/eigenMT.output.500k.win_$WINDOW.allchrs.txt
done

# Compare Rasqual p values corrected by bonferroni vs. by EigenMT
Rscript ~/src/utils/rasqual/compareRasqualEigenMT.R output/rasqual.500k.all.lead.txt output/eigenMT.output.500k.win_200.allchrs.txt


########################## Get genotypes (VCF) of lead SNPs ######################

cat rasqual.500k.97samples.leadSNPs.fdr0.1.chrpos.txt | cut -f 1-2 > rasqual.500k.97samples.leadSNPs.regionsFile.txt
head rasqual.500k.97samples.leadSNPs.regionsFile.txt > rasqual.500k.97samples.leadSNPs.regionsFile.txt.head

VCFBASE=imputed.97_samples
VCF_QTL=$VCFBASE.snps_indels.INFO_08.MAF_0.05.RASQUAL.vcf.gz
bcftools view -R rasqual.500k.97samples.leadSNPs.regionsFile.txt -O v -o imputed.97_samples.leadSNPs.tmp.vcf input/$VCF_QTL

# Note that some SNPs are eQTLs for more than one gene, i.e. there are duplicates
# Also, bcftools view returns indels that overlap the lead SNP positions we specified,
# which is not what we want and needs to be corrected.
zcat input/$VCFBASE.snps_indels.INFO_08.MAF_0.05.RASQUAL.vcf.gz | head -n 200 | grep "^##" > imputed.97_samples.leadSNPs.vcf
Rscript ../utils/getCorrectLeadSNPs.R imputed.97_samples.leadSNPs.tmp.vcf rasqual.500k.eigenMT.97samples.merged.leadSNPs.fdr0.1.txt >> imputed.97_samples.leadSNPs.vcf
grep -v "^#" imputed.97_samples.leadSNPs.vcf | wc -l

