#!/usr/bin/perl -w
use strict;
use warnings;
use Getopt::Long; 

sub prompt {
	print <<EOF;
	Usage: perl fa2fq.pl --fa <fasta> --qual <>
	- to convert fasta to fastq with fixed Phred33 or Phred64 quality
	[Options]
	--fa: fasta file; required
	--phred: coding method: 33 or 64; default=33
	--qual: a number of quality values from 0-40; default=20;
	--help: help information
EOF
exit 1;
}

my ($fa, $qual, $phred, $help);
&GetOptions("fa|f=s" => \$fa,
            "phred|p=i" => \$phred,
            "qual|q=i" => \$qual,
			"help|h" => \$help);

$phred = 33 if (!defined $phred);
$qual = 20 if (!defined $qual);

&prompt if $help or !defined $fa;

if ($phred != 33 and $phred != 64) {
	print STDERR "only 33 or 64 can be assigned to --phred\n";
	&prompt;
}

if ($qual > 40) {
	print STDERR "only 0-40 can be assigned to --qual\n";
	&prompt;
}

my $qual_char = chr($qual + $phred);

&prompt if $help;

if (!defined $fa) {
	print STDERR "--fa is required\n";
	&promp;
}


my ($name, $seq);

# check fasta file:
if (-z $fa) {
	print "$fa is empty\n";
	exit 1;
}

# read fasta and process data:
open(IN, $fa) || die;
while (<IN>) {
	chomp;
	if (/^>(.+)/) {
		if (defined $name) {
			print "\@$name\n$seq\n";
			my $qual_string = "";
			for (my $i=0; $i<length($seq); $i++) {
				$qual_string .= $qual_char;
			}
			print "\+\n$qual_string\n";
		}
		$name = $1;
		$seq = '';
	} else {
		$seq .= $_;
	}
}
print "\@$name\n$seq\n";
my $qual_string = "";
for (my $i=0; $i<length($seq); $i++) {
	$qual_string .= $qual_char;
}
print "\+\n$qual_string\n";
close IN;

