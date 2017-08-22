#!/usr/bin/perl -w
use strict;
use warnings;
use Getopt::Long;

my ($size, $excludeLong, $onlyLong);
GetOptions('size=i' => \$size,
           'excludelong=i' => \$excludeLong,
           'onlylong=i' => \$onlyLong) or die("Failed to get command line arguments");

my $MAX_GENES_IN_BATCH = 100;

# Bug: this won't ever treat the first gene if it is long running
my $curGenes = "";
my $curSum = 0;
my $numGenesInBatch = 1;
my $batchNum = 1;
# Get the first gene to include
while (<>) {
	chomp();
	my @l = split();	
	next if ($excludeLong and $l[9] > $excludeLong);
	next if ($onlyLong and $l[9] <= $onlyLong);
	$curGenes = $l[0];
	$curSum = $l[9];
	last;
}

while (<>) {
	chomp();
	my @l = split();	
	next if ($excludeLong and $l[9] > $excludeLong);
	next if ($onlyLong and $l[9] <= $onlyLong);
	$curSum += $l[9];
	if ($curSum > $size or $numGenesInBatch >= $MAX_GENES_IN_BATCH) {		
		print "Batch".$batchNum."\t".$curGenes."\n";
		$curSum = $l[9];
		$curGenes = $l[0];
		$batchNum += 1;
		$numGenesInBatch = 1;
	} else {
		$curGenes = $curGenes.",".$l[0];
		$numGenesInBatch += 1;
	}
}

if ($curGenes) {
	print "Batch".$batchNum."\t".$curGenes."\n";
}
