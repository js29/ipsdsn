#!/usr/bin/perl -w

use strict;
use Getopt::Long;
use Common qw(open_for_read2 open_for_write2 verbose checkFileExists checkDirExists);

sub usage($);
my $usage = "$0 - Get the lead SNP and the p value difference to the second-best SNP.
    usage:
       $0 [OPTIONS] <inFile|stdin>
    options:
       -verbose=XXX - set verbosity level
       -XXX\n";

## Read options from CL
my $FvalField = 11;
my ($inFile, $statcol, $pcol, $genecol);
GetOptions(	'statcol=i' => \$statcol,
			'genecol=i' => \$genecol,
			'pcol=i' => \$pcol,
			'file|f=s' => \$inFile);

(!defined $inFile) and usage("Error: required argument --file is missing.\n");
(!defined $pcol and !defined $statcol) and usage("Error: either --pcol or --statcol argument is required.\n");
(!defined $genecol) and usage("Error: --genecol argument is required.\n");

(defined $genecol) and $genecol -= 1;
(defined $pcol) and $pcol -= 1;
(defined $statcol) and $statcol -= 1;

## Check CL args
#&checkFileExists($inFile) unless $inFile eq 'stdin';

## Read input files
my @bestLine;
my @secondBestLine;
my $curGene = 'NA';
$" = "\t";
*IN = &open_for_read2($inFile);
while(<IN>) {
	chomp();
	my @l = split(/\t/, $_, -1);
	if ($curGene ne 'NA') {
		if ($curGene eq $l[$genecol]) {
			if (defined $pcol) {
				if ($l[$pcol] < $bestLine[$pcol]) {
					@secondBestLine = ();
					@secondBestLine = @bestLine;
					@bestLine = @l;
				} elsif (scalar(@secondBestLine) < 2 or $l[$pcol] < $secondBestLine[$pcol]) {
					@secondBestLine = @l;
				}
			} else {
				if ($l[$statcol] > $bestLine[$statcol]) {
					@secondBestLine = ();
					@secondBestLine = @bestLine;
					@bestLine = @l;
				} elsif (scalar(@secondBestLine) < 2 or $l[$statcol] > $secondBestLine[$statcol]) {
					@secondBestLine = @l;
				}
			}
		} else {
			print STDERR ($l[$genecol]."\n");
			my $secondBest = 1;
			if (scalar(@secondBestLine) > 1) {
				if (defined $pcol) {
					$secondBest = $secondBestLine[$pcol];
				} else {
					$secondBest = $secondBestLine[$statcol];
				}
			}
			print "@bestLine\t$secondBest\n";
			@bestLine = ();
			@bestLine = @l;
			@secondBestLine = ();
			$curGene = $l[$genecol];
		}
	} else {
		print STDERR ($l[$genecol]."\n");
		@bestLine = @l;
		@secondBestLine = ();
		$curGene = $l[$genecol];
	}
}
close(IN) unless $inFile eq 'stdin'; ## Avoids perl warning if we try and close STDIN;

my $secondBest = 1;
if (scalar(@secondBestLine) > 1) {
	if (defined $pcol) {
		$secondBest = $secondBestLine[$pcol];
	} else {
		$secondBest = $secondBestLine[$statcol];
	}
}
print "@bestLine\t$secondBest\n";

###############################################################################
## Functions
sub usage($) {
	print $_[0]."\n";
    print STDERR $usage;
    exit;
}

