#!/bin/bash
# This file is only meant to be used as an example. Similar commands can be run one at
# a time on the command line.
STARDIR=STAR/
cd leafcutter

# Runs leafcutter's bam2junc.sh script for each sample and outputs the paths in 
# all.juncfiles.txt.
cut -f1 ../data/metadata.qtl_samples.txt | sed '1d' | submitJobs.py --MEM 200 --jobname bam2junc --command  "python ../utils/leafcutter.bam2junc.py --indir $STARDIR --insuffix .qtl.bam --juncfilelist all.juncfiles.txt"
grep "Successfully completed" FarmOut/bam2junc.*.txt | wc -l
grep -i "Failed|TERM|error" FarmOut/bam2junc.*.txt | wc -l

submitJobs.py --MEM 1000 -q yesterday -j leafcutter_cluster -c "python ~/software/leafcutter/clustering/leafcutter_cluster.py -j all.juncfiles.txt -o snqtl.v5.clusters"

# Fix header of clusters output so that it is readable by data.table fread
(echo -n "chrom "; zcat snqtl.v5.clusters_perind_numers.counts.gz | head -n 1) | gzip > snqtl.v5.clusters_perind_numers.counts.fixed.gz
zcat snqtl.v5.clusters_perind_numers.counts.gz | sed '1d' | gzip >> snqtl.v5.clusters_perind_numers.counts.fixed.gz

submitJobs.py --MEM 3000 -q yesterday -j filterClusterCounts -c "Rscript ../utils/leafcutter.filterClusterCounts.R --file snqtl.v5.clusters_perind_numers.counts.fixed.gz --min_intron_fraction 0.02 --outputroot snqtl.v5.clusters_perind_numers"
gzip snqtl.v5.clusters_perind_numers.filtered_counts.txt snqtl.v5.clusters_perind_numers.ratios.txt

# Count how many clusters there are
zcat snqtl.v5.clusters_perind_numers.counts.gz | sed '1d' | cut -f 1 -d ' ' | cut -f 4 -d ':' | sort | uniq | wc -l
zcat snqtl.v5.clusters_perind_numers.filtered_counts.txt.gz | sed '1d' | cut -f 1 | cut -f 4 -d ':' | sort | uniq | wc -l

# Plot the distribution of intron sizes
Rscript ../utils/leafcutter.plotIntronSizes.R --file snqtl.v5.clusters_perind_numers.counts.fixed.gz -o snqtl.v5.clusters.intronsizes
# This shows that 96% of introns are < 50 kb in size. So if we place our "phenotype" at
# the center of the intron and use a 25 kb window, we keep 96% of intron/exon boundaries.
# And 90% of introns are < 25 kb in size.


# This writes out the normalized ratio matrix, as well as the PCA of this matrix
# This also sets the column IDs to be the HIPSCI genotype IDs, important for
# later use in FastQTL
zcat snqtl.v5.clusters_perind_numers.ratios.txt.gz | grep -v "^Y" | gzip > snqtl.v5.clusters_perind_numers.ratios.nochrY.txt.gz
submitJobs.py --MEM 2000 -q yesterday -j leafcutter_normalize -c "Rscript ../utils/leafcutter.normalize.R snqtl.v5.clusters_perind_numers.ratios.nochrY.txt.gz ../data/metadata.qtl_samples.txt snqtl.v5.ratios"


###############################################################################
# Run FastQTL on the normalized ratios for each intron
BASE=imputed.97_samples
VCF_QTL=../data/$BASE.snps_indels.INFO_08.MAF_0.05.vcf.gz
NAME=snqtl.v5.ratios

mkdir fastqtl
cd fastqtl
mkdir input
mkdir output

F=../$NAME.normalized.txt
head -n 1 $F > input/$NAME.normalized.sorted.txt
sed '1d' $F | sort -k1,1 -k2,2n >> input/$NAME.normalized.sorted.txt
bgzip input/$NAME.normalized.sorted.txt && tabix -p bed input/$NAME.normalized.sorted.txt.gz

for i in {1..30}; do
   head -n $(($i+1)) ../$NAME.pcs.txt > input/clusters.allnonNA.PCs.$i.txt
done


WINDOW=15000
i=0
OUTNAME=fastqtl.nominals.PCs
cat input/fastQTL.50chunks.txt | submitJobs.py --MEM 300 -j $OUTNAME.$i -q normal -c "python ../utils/runFastQTL.py --out output/$OUTNAME.$i --vcf $VCF_QTL --bed input/$NAME.normalized.sorted.txt.gz --window $WINDOW --threshold 0.01 --execute True"
for i in {1..30}; do
    COV=input/clusters.allnonNA.PCs.$i.txt
    cat input/fastQTL.50chunks.txt | submitJobs.py --MEM 300 -j $OUTNAME.$i -q normal -c "python ../utils/runFastQTL.py --out output/$OUTNAME.$i --cov $COV --vcf $VCF_QTL --bed input/$NAME.normalized.sorted.txt.gz --window $WINDOW --threshold 0.01 --execute True"
done

# Check output to see if any jobs failed
grep -i TERM FarmOut/$OUTNAME.*.txt > FarmOut.$OUTNAME.TERM.txt
grep -iP "error|fault|failed" FarmOut/$OUTNAME.*.txt > FarmOut.$OUTNAME.error.txt
grep "Successfully completed" FarmOut/$OUTNAME.*.txt | wc -l

for i in {0..30}; do
    zcat output/$OUTNAME.$i.chunk_*.txt.gz | bgzip > output/$OUTNAME.$i.pvalues.txt.gz
    cat FarmOut/$OUTNAME.$i.*.txt > FarmOut/all.$OUTNAME.$i.txt
done
rm output/$OUTNAME.*.chunk_*.txt.gz
rm FarmOut/$OUTNAME.*.txt

# Count eQTLs
for f in output/$OUTNAME.*.pvalues.txt.gz
do
    (echo -ne "$f\t"; zcat $f | awk '$4 < 1e-5 { print $1 }' | sort | uniq | wc -l) >> sqtl_counts.1e-5.txt
done
# This shows a clear maximum in number of QTLs with p<1e-5 when 5 PCs are used (~2500).
# We'll just use 5 PCs for our permutations.

i=5
submitJobs.py --MEM 5000 -j fastQTL_add_coords -c "python ../utils/fastqtlAddSnpCoordinates.py --vcf $VCF_QTL --fastqtl output/$OUTNAME.$i.pvalues.txt.gz | bgzip > output/$OUTNAME.$i.pvalues.coords.txt.gz"


###############################################################################
# Run FastQTL WITH PERMUTATIONS
i=5
WINDOW=15000
OUTNAME=permutations.10k.PCs

COV=input/clusters.allnonNA.PCs.$i.txt
cat input/fastQTL.200chunks.txt | submitJobs.py --MEM 500 -j $OUTNAME.$i -q normal -c "python ../utils/runFastQTL.py --permute 10000 --out output/$OUTNAME.$i --cov $COV --vcf $VCF_QTL --bed input/$NAME.normalized.sorted.txt.gz --window $WINDOW --execute True"

# Check output to see if any jobs failed
grep -i TERM FarmOut/$OUTNAME.*.txt > FarmOut.$OUTNAME.TERM.txt
grep -iP "error|fault|failed" FarmOut/$OUTNAME.*.txt > FarmOut.$OUTNAME.error.txt
grep "Successfully completed" FarmOut/$OUTNAME.*.txt | wc -l

zcat output/$OUTNAME.$i.chunk_*.txt.gz | bgzip > output/$OUTNAME.$i.pvalues.txt.gz
cat FarmOut/$OUTNAME.$i.*.txt > FarmOut/all.$OUTNAME.$i.txt
rm output/$OUTNAME.*.chunk_*.txt.gz
rm FarmOut/$OUTNAME.*.txt

submitJobs.py --MEM 5000 -j fastQTL_add_coords -c "python ../utils/fastqtlAddSnpCoordinates.py --vcf $VCF_QTL --fastqtl output/$OUTNAME.$i.pvalues.txt.gz --snpcol 6 | bgzip > output/$OUTNAME.$i.pvalues.coords.txt.gz"


# See how many clusters we had as input to FastQTL
zcat input/snqtl.v5.ratios.normalized.sorted.txt.gz | sed '1d' | cut -f 4 | tr ':' '\t' | cut -f 4 | sort | uniq | wc -l
# 30591

# And how many as output
zcat output/permutations.10k.PCs.$i.pvalues.coords.txt.gz | tr ' ' '\t' | cut -f 1 | tr ':' '\t' | cut -f 4 | sort | uniq | wc -l
# 30591
# Noted in the R script below, bonferroniCorrectClusters.R:
# Strangely, we have some duplicates of cluster introns / rsid combinations.
# This shouldn't happen, but my best guess is that it's due to there being
# SNPs with duplicate IDs in the VCF. To get around this we just remove the
# relatively small number of duplicates.


# Get significant sQTLs at FDR 10%
# First get the min P value of any intron per cluster, and do Bonferroni correction
# on these min P values.
Rscript ../utils/leafcutter.bonferroniCorrectClusters.R output/$OUTNAME.$i.pvalues.coords.txt.gz > output/$OUTNAME.$i.pvalues.coords.bonf.txt

Rscript ../utils/addFDRCol.R output/$OUTNAME.$i.pvalues.coords.bonf.txt 14 TRUE | sort -k14,14g > output/$OUTNAME.$i.fdr.txt
cat output/$OUTNAME.$i.fdr.txt | awk '$16 < 0.1' > output/$OUTNAME.$i.fdr0.1.txt
wc -l output/$OUTNAME.$i.fdr0.1.txt
# 2079 lead sQTLs at FDR 10%


###############################################################################
# Annotate sQTLs (genes, splice site dist)

cd leafcutter/fastqtl
i=5

NAME=sqtls.fdr0.1
# Get genes associated with the intron clusters
cp output/permutations.10k.PCs.$i.fdr0.1.txt $NAME.txt

cat $NAME.txt | perl -ane '@l = split(/:/, $F[0]); print join("\t", @l[0..2])."\n";' > $NAME.introns.bed

bedtools closest -a $NAME.introns.bed -b ../../annotations/Homo_sapiens.GRCh38.84.gene_start_end.bed > $NAME.introns.withgenes.tmp.txt

echo -e "chr\tstart\tend\tgene_chr\tgene_start\tgene_end\tgene\tsymbols" > $NAME.introns.withgenes.tmp.txt
Rscript ../utils/collapseCol.R --file $NAME.introns.withgenes.tmp.txt --dupcols "1,2,3" --collapsecols "7,8" --sep "," >> $NAME.introns.withgenes.tmp2.txt

cat ../../annotations/Homo_sapiens.GRCh38.84.exon_start_end.txt | awk 'BEGIN{OFS="\t"}{print $2,$3,$4,$1}' > ../../annotations/Homo_sapiens.GRCh38.84.exon_start_end.bed

Rscript ../utils/leafcutter.refineOverlappingGenes.R $NAME.introns.withgenes.tmp2.txt ../../annotations/Homo_sapiens.GRCh38.84.exon_start_end.bed > $NAME.introns.withgenes.tmp3.txt

paste <(echo -e "cluster\tnumtested\tmleshape1\tmleshape2\tdummy\trsid\tdist\tnominal_pval\tslope\tperm_pval\tbeta_pval\tchr\tpos\tbonf_pval\tclustersize\tfdr") <(head -n 1 $NAME.introns.withgenes.tmp3.txt | cut -f 7,8) > $NAME.genes.txt
paste $NAME.txt <(cut -f 7,8 $NAME.introns.withgenes.tmp3.txt | sed '1d') >> $NAME.genes.txt

cat $NAME.genes.txt | awk 'BEGIN{OFS="\t"}{print $1,$2,$6,$12,$13,$7,$11,$15,$16,$17,$18}' > $NAME.genes.cut.txt

rm $NAME.introns.withgenes.tmp.txt $NAME.introns.withgenes.tmp2.txt $NAME.introns.withgenes.tmp3.txt

# Now annotate other info - DRG expression, literature associations, etc.


########################## Get genotypes (VCF) of lead SNPs ######################

cat $NAME.txt | cut -f 12-13 > $NAME.leadSNPs.regionsFile.txt

VCF_QTL=../data/imputed.97_samples.snps_indels.INFO_08.MAF_0.05.RASQUAL.vcf.gz
bcftools view -R $NAME.leadSNPs.regionsFile.txt -O v -o imputed.97_samples.leadSNPs.tmp.vcf $VCF_QTL

# Note that some SNPs are eQTLs for more than one gene, i.e. there are duplicates
# Also, bcftools view returns indels that overlap the lead SNP positions we specified,
# which is not what we want and needs to be corrected.
zcat $VCF_QTL | head -n 200 | grep "^##" > imputed.97_samples.sqtl.leadSNPs.vcf
Rscript ../utils/getCorrectLeadSNPs.R imputed.97_samples.leadSNPs.tmp.vcf $NAME.txt >> imputed.97_samples.sqtl.leadSNPs.vcf
grep -v "^#" imputed.97_samples.sqtl.leadSNPs.vcf | wc -l
# 2059 unique lead SNPs for 2079 splicing clusters

