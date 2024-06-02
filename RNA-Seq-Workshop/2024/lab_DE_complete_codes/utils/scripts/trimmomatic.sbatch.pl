#!/usr/bin/perl -w
use strict;
use warnings;
use Getopt::Long;

my ($trimfilecol, $java_version, $mem, $time, $trim, $adpt, $indir, $outdir, $help);
my ($trimmomatic_jar, $fq1feature, $fq2feature, $threads, $min_len);
my $result = &GetOptions("java=s" => \$java_version,
			 "mem=s" => \$mem,
			 "time=s" => \$time,
			 "trim_shell=s" => \$trim,
			 "adaptor_file=s" => \$adpt,
			 "indir=s" => \$indir,
			 "outdir=s" => \$outdir,
			 "trimmomatic=s" => \$trimmomatic_jar,
			 "fq1feature=s" => \$fq1feature,
			 "fq2feature=s" => \$fq2feature,
			 "threads=i" => \$threads,
			 "min_len=i" => \$min_len,
			 "help|h" => \$help
);

# print help information if errors occur:
if ($help) {
	&errINF;
	exit;
}

$java_version = "Java" if (!defined $java_version);
$mem = "16G" if (!defined $mem);
$time = "12:00:00" if (!defined $time);
$trimmomatic_jar = "/homes/liu3zhen/software/trimmomatic/Trimmomatic-0.38/trimmomatic-0.38.jar" if (!defined $trimmomatic_jar);
$fq1feature = "_1.fastq" if (!defined $fq1feature);
$fq2feature = "_2.fastq" if (!defined $fq2feature);
$outdir = "." if (!defined $outdir);
$threads = 1 if (!defined $threads);
$min_len = 60 if (!defined $min_len);

open (IN,"ls \"$indir\" -1 |");
while (<IN>) {
	chomp;
	my $trimfile = $_;
	if ($trimfile =~ $fq1feature) {
		print "$trimfile\n";
		my $outfile = $trimfile.".sbatch";
		open(OUT, ">", $outfile) || die;
		print OUT "#!/bin/bash -l\n";
		print OUT "#SBATCH --mem-per-cpu=$mem\n";
		print OUT "#SBATCH --time=$time\n";
		print OUT "#SBATCH --cpus-per-task=$threads\n";
		printf OUT "module load $java_version\n";
		print OUT "bash $trim \\\n";
		print OUT "$trimmomatic_jar \\\n";
		print OUT "$adpt \\\n";
		print OUT "$indir \\\n";
		print OUT "$outdir \\\n";
		print OUT "$fq1feature $fq2feature \\\n";
		print OUT "$threads $min_len $trimfile";
		close OUT;
		my $sbatch_cmd = sprintf("sbatch %s", $outfile);
		print "$sbatch_cmd\n";
		system($sbatch_cmd);
	}
}

close IN;

