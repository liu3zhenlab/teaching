#!/usr/bin/perl -w
use strict;
use warnings;
use Getopt::Long;

my ($fasterq_dump, $mem, $time, $intable, $srrcol, $fq_para, $threads, $help);
my $result = &GetOptions(
             "fdpath=s" => \$fasterq_dump,
             "mem=s" => \$mem,
             "time=s" => \$time,
			 "threads=i" => \$threads,
			 "in=s" => \$intable,
			 "srrcol=i" => \$srrcol,
			 "fqpara=s" => \$fq_para,
			 "help|h" => \$help
);

# print help information if errors occur:
if ($help or !defined $intable) {
	&errINF;
	exit;
}

if (! defined $fasterq_dump) {
	$fasterq_dump = "/homes/liu3zhen/software/sra_tools/current/fasterq-dump";
}

$fq_para = "--split-files" if (! defined $fq_para);
if (!defined $mem) {
	$mem = "4G"
} else {
	$mem = $mem."G";
}
$time = "1-00:00:00" if (!defined $time);
$threads = 1 if (!defined $threads);
$srrcol = 16 if (!defined $srrcol); 
&errINF if (! defined $intable);

open (IN, $intable) || die;
while (<IN>) {
	chomp;
	my @line = split(/\t/, $_);
	my $srr = $line[$srrcol - 1];
	if ($srr =~ /[DES]RR/) {
		print "$srr\n";
		my $outfile = $srr.".sbatch";
		open(OUT, ">", $outfile) || die;
		print OUT "#!/bin/bash -l\n";
		print OUT "#SBATCH --cpus-per-task=$threads\n";
		print OUT "#SBATCH --mem-per-cpu=$mem\n";
		print OUT "#SBATCH --time=$time\n";
		print OUT "$fasterq_dump -f $fq_para $srr\n";
		close OUT;
		my $sbatch_cmd = sprintf("sbatch %s", $outfile);
		print "$sbatch_cmd\n";
		system($sbatch_cmd);
	}
}

close IN;

sub errINF {

	print <<EOF;
    Usage: perl fasterq_dump --fdpath <path-to-fasterq_dump> --in <meta> [options]
    [Options]
    --in <file>     tab-delimit flat file containing a column for SRA accessions
                    (e.g., SRRxxxxx ERRxxxxx, Dxxxxx); required
    --fdpath <path> path-to-fasterq_dump ("/homes/liu3zhen/software/sra_tools/current/fasterq-dump")
    --mem <num>     Gb memory per thread requested (4)
    --time <time>   running time requested; (0-23:00:00)
    --threads <num> number of threads (1)
    --srrcol <num>  column number of SRA accession ID in meta file (4)
    --fqpara <para> parameters to pass to fasterq_dump (--split-files)
    --help
EOF
exit;
}


