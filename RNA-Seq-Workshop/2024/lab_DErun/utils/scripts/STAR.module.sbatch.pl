#!/usr/bin/perl -w
use strict;
use warnings;
use Getopt::Long;

my $outdir=".";
my ($fqfilecol, $mem, $time, $indir, $dbdir, $help);
my ($fq1feature, $fq2feature, $threads, $module);
my ($alignIntronMax, $alignMatesGapMax, $outSAMattrIHstart);
my ($outSAMmultNmax, $outSAMstrandField, $outFilterIntronMotifs, $outSAMtype);
my ($quantMode, $outFilterMismatchNmax, $outFilterMismatchNoverLmax, $outFilterMatchNmin);
my ($outSJfilterReads, $outFilterMultimapNmax, $outFilterMultimapScoreRange, $outFilterMatchNminOverLread);
my ($outSAMmapqUnique);

my $result = &GetOptions("module=s" => \$module,
                         "mem=i" => \$mem,
                         "time=s" => \$time,
			 "indir=s" => \$indir,
			 "outdir=s" => \$outdir,
			 "dbdir=s" => \$dbdir,
			 "fq1feature=s" => \$fq1feature,
			 "fq2feature=s" => \$fq2feature,
			 "threads=i" => \$threads,
			 "alignIntronMax=i" => \$alignIntronMax,
			 "alignMatesGapMax=i" => \$alignMatesGapMax,
			 "outSAMattrIHstart=i" => \$outSAMattrIHstart,
			 "outSAMmultNmax=i" => \$outSAMmultNmax,
			 "outSAMmapqUnique=i" => \$outSAMmapqUnique,
			 "outSAMstrandField=s" => \$outSAMstrandField,
			 "outFilterIntronMotifs=s" => \$outFilterIntronMotifs,
			 "outSAMtype=s" => \$outSAMtype,
			 "quantMode=s" => \$quantMode,
			 "outFilterMismatchNmax=i" => \$outFilterMismatchNmax,
			 "outFilterMismatchNoverLmax=f" => \$outFilterMismatchNoverLmax,
			 "outFilterMatchNmin=i" => \$outFilterMatchNmin,
			 "outFilterMatchNminOverLread=f" => \$outFilterMatchNminOverLread,
			 "outSJfilterReads=s" => \$outSJfilterReads,
			 "outFilterMultimapNmax=i" => \$outFilterMultimapNmax,
			 "outFilterMultimapScoreRange=i" => \$outFilterMultimapScoreRange,
			 "help|h" => \$help
);

# print help information if errors occur:
if ($help) {
	&errINF;
	exit;
}

$module = "STAR" if (!defined $module);
$mem = 6 if (!defined $mem);
$time = "24:00:00" if (!defined $time);
$fq1feature = ".R1.pair.fq" if (!defined $fq1feature);
$fq2feature = ".R2.pair.fq" if (!defined $fq2feature);
$threads = 8 if (!defined $threads);
$alignIntronMax = 100000 if (!defined $alignIntronMax);
$alignMatesGapMax = 100000 if (!defined $alignMatesGapMax);
$outSAMattrIHstart = 0 if (!defined $outSAMattrIHstart);
$outSAMmultNmax = 1 if (!defined $outSAMmultNmax);
$outSAMstrandField = "intronMotif" if (!defined $outSAMstrandField);
$outFilterIntronMotifs = "RemoveNoncanonicalUnannotated" if (!defined $outFilterIntronMotifs);
$outSAMtype = "BAM SortedByCoordinate" if (!defined $outSAMtype);
$quantMode = "GeneCounts" if (!defined $quantMode);
$outFilterMismatchNmax = 2 if (!defined $outFilterMismatchNmax);
$outFilterMismatchNoverLmax = 0.02 if (!defined $outFilterMismatchNoverLmax);
$outFilterMatchNmin = 50 if (!defined $outFilterMatchNmin);
$outSJfilterReads = "Unique" if (!defined $outSJfilterReads);
$outFilterMultimapNmax = 1 if (!defined $outFilterMultimapNmax);
$outFilterMultimapScoreRange = 2 if (!defined $outFilterMultimapScoreRange);
$outFilterMatchNminOverLread = 0.95 if (!defined $outFilterMatchNminOverLread);
$outSAMmapqUnique = 60 if (!defined $outSAMmapqUnique);

open (IN,"ls \"$indir\" -1 |");
while (<IN>) {
	chomp;
	my $fqfile = $_;
	if ($fqfile =~ $fq1feature) {
		my $sample = $fqfile;
		$sample =~ s/$fq1feature//g;
		my $fq1 = $fqfile;
		my $fq2 = $fq1;
		$fq2 =~ s/$fq1feature/$fq2feature/g;
		print "$fqfile\n";
		print "$sample\n";
		my $outfile = $sample.".sbatch";
		open(OUT, ">", $outfile) || die;
		print OUT "#!/bin/bash -l\n";
		print OUT "#SBATCH --mem-per-cpu=$mem"."G\n";
		print OUT "#SBATCH --time=$time\n";
		print OUT "#SBATCH --nodes=1\n";
		print OUT "#SBATCH --ntasks-per-node=$threads\n";
		printf OUT "module load $module\n";
		
		print OUT "STAR \-\-runThreadN $threads \\
			\-\-genomeDir $dbdir \\
			\-\-readFilesIn $indir\/$fq1 $indir\/$fq2 \\
			\-\-alignIntronMax $alignIntronMax \\
			\-\-alignMatesGapMax $alignMatesGapMax \\
			\-\-outFileNamePrefix $outdir\/$sample \\
			\-\-outSAMattrIHstart $outSAMattrIHstart \\
			\-\-outSAMmultNmax $outSAMmultNmax \\
			\-\-outSAMstrandField $outSAMstrandField \\
			\-\-outFilterIntronMotifs $outFilterIntronMotifs \\
			\-\-outSAMtype $outSAMtype \\
			\-\-quantMode $quantMode \\
			\-\-outFilterMismatchNmax $outFilterMismatchNmax \\
			\-\-outFilterMismatchNoverLmax $outFilterMismatchNoverLmax \\
			\-\-outFilterMatchNmin $outFilterMatchNmin \\
			\-\-outSJfilterReads $outSJfilterReads \\
			\-\-outFilterMultimapNmax $outFilterMultimapNmax \\
			\-\-outFilterMultimapScoreRange $outFilterMultimapScoreRange \\
			\-\-outSAMmapqUnique $outSAMmapqUnique \\
			\-\-outFilterMatchNminOverLread $outFilterMatchNminOverLread";

		close OUT;
		my $sbatch_cmd = sprintf("sbatch %s", $outfile);
		print "$sbatch_cmd\n";
		system($sbatch_cmd);
	}
}

close IN;

sub errINF {
	# to be added (SL)
}
