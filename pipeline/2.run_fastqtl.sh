#!/bin/bash
# This file is only meant to be used as an example. Similar commands can be run one at
# a time on the command line.
DATADIR=data
BASE=imputed.97_samples
OUTPUTDIR=fastqtl

# This VCF file not provided publicly
VCF_QTL=data/$BASE.snps_indels.INFO_08.MAF_0.05.vcf.gz

cut -f 2 --complement data/all_basic_counts.qtl.txt > $BASE.basic_counts.txt

# Get phenotype data (gene counts) for FastQTL. Ignore chrY and MT for now.
Rscript utils/getFastqtlPhenotypes.R data/combined_expression_data.qtl.rds ExpressedOnly | grep -v -P "^Y|^MT" > fastqtl/input/$BASE.fastqtl.phenotypes.expressedgenes.bed
Rscript utils/getFastqtlPhenotypes.R data/combined_expression_data.qtl.rds | grep -v -P "^Y|^MT" > fastqtl/input/$BASE.fastqtl.phenotypes.allgenes.bed
bgzip fastqtl/input/$BASE.fastqtl.phenotypes.expressedgenes.bed && tabix -p bed fastqtl/input/$BASE.fastqtl.phenotypes.expressedgenes.bed.gz
bgzip fastqtl/input/$BASE.fastqtl.phenotypes.allgenes.bed && tabix -p bed fastqtl/input/$BASE.fastqtl.phenotypes.allgenes.bed.gz

# Get principal components from the expression matrix, using either all genes
# or only expressed genes
Rscript $SN/eqtl/R/getExpressionPCs.R $SN/eqtl/combined_expression_data.v5.qtl.rds ExpressedOnly > $SN/eqtl/combined_expression_data.v5.qtl.expressed_genes.PCs.hipsciGT.txt
Rscript utils/getExpressionPCs.R data/combined_expression_data.qtl.rds > fastqtl/input/combined_expression_data.qtl.allgenes.PCs.hipsciGT.txt
for i in {1..30}; do
   head -n $(($i+1)) $SN/eqtl/combined_expression_data.v5.qtl.expressed_genes.PCs.hipsciGT.txt > input/expressed_genes.PCs.$i.txt
   head -n $(($i+1)) $SN/eqtl/combined_expression_data.qtl.allgenes.PCs.hipsciGT.txt > fastqtl/input/all_genes.PCs.$i.txt
done

Rscript utils/writePCSelectionFiles.R 30 fastqtl/input/includePCs



######### Run FastQTL
cd fastqtl

# Run FastQTL for all genes (~35000) with PCs generated using expressed genes
# File fastQTL.25chunks.txt is a file with 25 lines like this:
# 1 25
# 2 25
# ...
# 25 25
mkdir output
i=0
cat input/fastQTL.25chunks.txt | submitJobs.py --MEM 200 -j fastqtl.allgenes.ePCs.$i -q normal -c "python ../utils/runFastQTL.py --out output/nominals.allgenes.ePCs.$i --vcf $VCF_QTL --bed input/$BASE.fastqtl.phenotypes.allgenes.bed.gz --window 500000 --threshold 0.01 --execute True"
for i in {1..30}; do
    COV=input/expressed_genes.PCs.$i.txt
    #COV=input/all_genes.PCs.$i.txt
    cat input/fastQTL.25chunks.txt | submitJobs.py --MEM 200 -j fastqtl.allgenes.ePCs.$i -q normal -c "python ../utils/runFastQTL.py --out output/nominals.allgenes.ePCs.$i --cov $COV --vcf $VCF_QTL --bed input/$BASE.fastqtl.phenotypes.allgenes.bed.gz --window 500000 --threshold 0.01 --execute True"
done

# Check output to see if any jobs failed
grep TERM FarmOut/fastqtl.allgenes.ePCs.*.txt | wc -l
grep -i error FarmOut/fastqtl.allgenes.ePCs.*.txt | wc -l
grep -i "successfully completed" FarmOut/fastqtl.allgenes.ePCs.*.txt | wc -l

for i in {0..30}; do
    zcat output/nominals.allgenes.ePCs.$i.chunk_*.txt.gz | bgzip > output/nominals.allgenes.ePCs.$i.pvalues.txt.gz
    cat FarmOut/fastqtl.allgenes.ePCs.$i.*.txt > FarmOut/fastqtl.allgenes.ePCs.$i.all.txt
done
rm output/nominals.allgenes.ePCs.*.chunk_*.txt.gz
rm FarmOut/fastqtl.allgenes.ePCs.

# Count eQTLs
for f in output/nominals.allgenes.ePCs.*.pvalues.txt.gz
do
    (echo -ne "$f\t"; zcat $f | awk '$4 < 1e-5 { print $1 }' | sort | uniq | wc -l) >> eqtl_counts.imputed.allgenes.ePCs.1e-5.txt
done

# FastQTL doesn't include the coordinates of SNPs in its output files. So we get these
# from the VCF and add them in to the FastQTL output files.
i=20
submitJobs.py --MEM 5000 -j fastQTL_add_coords -c "python ../utils/fastqtlAddSnpCoordinates.py --vcf $VCF_QTL --fastqtl output/nominals.ePCs.$i.pvalues.txt.gz | bgzip > output/nominals.ePCs.$i.pvalues.coords.txt.gz"


# Run FastQTL for all genes (~35000) with PCs generated using expressed genes only
# WITH PERMUTATIONS
i=0
OUTNAME=permutations.10k.allgenes.PCs
cat input/fastQTL.25chunks.txt | submitJobs.py --MEM 500 -j $OUTNAME.$i -q normal -c "python ../utils/runFastQTL.py --permute 10000 --out output/$OUTNAME.$i --vcf $VCF_QTL --bed input/$BASE.fastqtl.phenotypes.allgenes.bed.gz --window 500000 --execute True"
for i in {1..30}; do
    COV=input/expressed_genes.PCs.$i.txt
    cat input/fastQTL.25chunks.txt | submitJobs.py --MEM 500 -j $OUTNAME.$i -q normal -c "python ../utils/runFastQTL.py --permute 10000 --out output/$OUTNAME.$i --cov $COV --vcf $VCF_QTL --bed input/$BASE.fastqtl.phenotypes.allgenes.bed.gz --window 500000 --execute True"
done

# Check output to see if any jobs failed
grep TERM FarmOut/$OUTNAME.*.txt | wc -l
grep -i error FarmOut/$OUTNAME.*.txt | wc -l
grep -i "successfully completed" FarmOut/$OUTNAME.*.txt | wc -l

for i in {0..30}; do
    zcat output/$OUTNAME.$i.chunk_*.txt.gz | bgzip > output/$OUTNAME.$i.pvalues.txt.gz
    cat FarmOut/$OUTNAME.$i.*.txt > FarmOut/all.$OUTNAME.$i.txt
done
rm output/$OUTNAME.*.chunk_*.txt.gz
rm FarmOut/$OUTNAME.*.txt


Rscript ../utils/addFDRCol.R $OUTNAME.20.pvalues.txt.gz 11 FALSE > $OUTNAME.20.pvalues.fdr.txt
# Count the number of QTLs at FDR 10% for different numbers of PCs
for f in output/$OUTNAME.*.pvalues.txt.gz
do
    (echo -ne "$f\t"; Rscript ../utils/addFDRCol.R $f 11 FALSE | sort -nk12,12 | awk '$12 < 0.1' | wc -l) >> eqtl_counts.$OUTNAME.fdr0.1.txt
done


# Add SNP coords to fastq output file
i=20
OUTNAME=nominals.allgenes.cis500k.ePCs
submitJobs.py --MEM 5000 -j fastQTL_add_coords -q yesterday -c "python ../utils/fastqtlAddSnpCoordinates.py --vcf $VCF_QTL --fastqtl output/$OUTNAME.$i.pvalues.txt.gz |  awk '$2 > 0' | bgzip > output/$OUTNAME.$i.pvalues.notskipped.coords.txt.gz"

