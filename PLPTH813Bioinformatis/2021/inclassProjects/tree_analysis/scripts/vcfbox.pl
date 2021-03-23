#!/use/bin/perl -w
# Author: Sanzhen Liu <liu3zhen@ksu.edu>
# 
use strict;
use warnings;
use Getopt::Std;

&master;
exit;

sub master {
	&usage if (@ARGV < 1);
	my $cmd = shift(@ARGV);
	my %modules = (summary    => \&summary,
				   recall     => \&recall,
				   allele     => \&allele,
				   genomatch  => \&genomatch,
				   sitefilt   => \&sitefilt,
				   taxafilt   => \&taxafilt,
				   compare    => \&compare,
				   select     => \&select,
				   recode     => \&recode
				   ); # to add more
	die ("Unknown module \"$cmd\".\n") if (!defined ($modules{$cmd}));
	&{$modules{$cmd}};
}


#########################################################################
# usage
#########################################################################
sub usage {
	die(qq/
Usage: vcfbox.pl <module> [arguments]\n
[Modules]
summary   : summarize variant types
recall    : recall genotypes based on read counts of alleles
            -vcf needs to have depth data for each taxon
allele    : convert genotyping to read counts of alleles for each taxon
            - vcf needs to have depth data for each taxon
genomatch : filter sites that have unmatched geno for select taxa
sitefilt  : filter markers\/sites based on #makers, rates of missing,
            homozygosity, and heterozygosity
taxafilt  : filter taxa based on #makers, missing rate, and heterozygosity
recode    : convert vcf to a standard table* and use defined genotype codes
compare   : compare genotypes of two selected taxa
select    : extract data of specified taxa

* stardard table:
	1st column:   chr
	2nd column:   position
	3rd column:   REF allele
	4th column:   ALT allele
	5th- columns: Genotypes of taxa
\n/);
}

#########################################################################
# summary
#########################################################################
sub summary {
	die(qq/Usage: vcfbox.pl summary <vcf>\n/) if (@ARGV==0 && -t STDIN);
	
	my %vartype_stat;
	my %chr_vartype_stat;
	while (<>) {
		if (!/^#/) {
			my @line = split("\t", $_);
			my $chr = $line[0]; # chr
			my ($ref, $alt) = @line[3..4]; # ref and alt alleles
			
			# common sequence at the beginning of both ref and alt
			my $ref_alt_common_head = &common_head_str($ref, $alt);
			
			# reformat ref and alt
			# remove common head sequence
			$ref =~ s/^$ref_alt_common_head//g;
			$ref = "." if $ref eq "";

			$alt =~ s/^$ref_alt_common_head//g;
			$alt = "." if $alt eq "";

			# type assessment
			my $ref_len = count_nuc($ref);
			my $alt_len = count_nuc($alt);

			if ($ref eq ".") { # insertion
				$vartype_stat{insert}++;
				$chr_vartype_stat{$chr}{insert}++
			} elsif ($alt eq ".") { # deletion
				$vartype_stat{deletion}++;
				$chr_vartype_stat{$chr}{deletion}++
			} elsif ($ref_len == 1 and $alt_len == 1) {
				$vartype_stat{substitution}++;
				$chr_vartype_stat{$chr}{substitution}++
			} elsif ($ref_len == $alt_len) {
				$vartype_stat{replace_equal}++;
				$chr_vartype_stat{$chr}{replace_equal}++
			} elsif ($ref_len < $alt_len) {
				$vartype_stat{replace_plus}++;
				$chr_vartype_stat{$chr}{replace_plus}++
			} elsif ($ref_len > $alt_len) {
				$vartype_stat{replace_minus}++;
				$chr_vartype_stat{$chr}{replace_minus}++
			} else {
				my $other_var = $ref."_to_".$alt;
				push(@{$vartype_stat{others}}, $other_var);
			}
		}
	}
	
	# report
	my @all_vartypes = sort {$vartype_stat{$b} <=> $vartype_stat{$a}} keys %vartype_stat;
	# print header and genomic summary data
	print "chr\t";
	print join("\t", @all_vartypes);
	print "\ttotal\nGenome";
	
	my $total_count = 0;
	foreach my $vartype (@all_vartypes) {
		$total_count += $vartype_stat{$vartype};
		print "\t$vartype_stat{$vartype}";
		print STDERR "$vartype\t$vartype_stat{$vartype}\n";
	}
	print "\t$total_count\n";
	
	# print summary data for each chr
	foreach my $echr (sort {$a cmp $b} keys %chr_vartype_stat) {
		print $echr;
		my $chr_total_count = 0;
		foreach my $vartype (@all_vartypes) {
			if (exists $chr_vartype_stat{$echr}{$vartype}) {
				print "\t$chr_vartype_stat{$echr}{$vartype}";
				$chr_total_count += $chr_vartype_stat{$echr}{$vartype};
			} else {
				print "\t0";
			}
		}
		print "\t$chr_total_count\n";
	}
}

sub count_nuc {
# count number of ATGC or atgc
	my $inseq = shift;
	$inseq =~ s/[^ATGCatgc]//g; # remove non-ATGC
	my $inseq_len = length($inseq);
	return $inseq_len;
}

sub common_head_str {
# determine shared sequences at the beginning of two input sequences
	my $share_head_bases = "";
	my ($inseq1, $inseq2) = @_;
	my @inseq1 = split(//, $inseq1);
	my @inseq2 = split(//, $inseq2);
	my $min_num_bases = ($#inseq1 >= $#inseq2) ? $#inseq2 : $#inseq1;
	for (my $i=0; $i<=$min_num_bases; $i++) {
		if ($inseq1[$i] eq $inseq2[$i]) {
			$share_head_bases .= $inseq1[$i];
		} else {
			last;
		}
	}

	return $share_head_bases;
}


#########################################################################
# genomatch
#########################################################################
sub genomatch {
	my %opts = ();
	&getopts('t:g:o:r', \%opts);

	die(qq/
	Usage: $0 genomatch [options] <vcf> 
	{arguments]
		-t: taxon\/sample names separated by comma; names need to match names in the vcf file; required
		-g: expected genotypes for entered taxa; number of geno must equal to number of taxa; required
		-r: if specified, selected taxa will be removed (unspecified by default)
		-o: output vcf file
	\n/) if (@ARGV==0 && -t STDIN);

	my $vcf = $ARGV[0];
	my $output;
	if (!defined $opts{o}) {
		$output = $vcf;
		$output =~ s/.vcf//g;
		$output .= ".genomatch.vcf";
	} else {
		$output = $opts{o};
	}
	open(OUT, ">", $output) || die;

	my ($taxa, $geno);

	if (!defined $opts{t} | !defined $opts{g}) {
		print STDERR "ERROR: -t and -g are required.\n";
		exit;
	} else {
		$taxa = $opts{t};
		$geno = $opts{g};
	}

	my @taxa = split(/,/, $taxa);
	my @geno = split(/,/, $geno);

	if ($#taxa != $#geno) {
		print STDERR "ERROR: number of taxa != number of geno.\n";
		exit;
	}

	my %inputgeno;
	my @tokeep_geno_cols;
	# open the file:
	open(IN, $vcf) || die;
	while (<IN>) {
		if (!/\#\#/) {
			chomp;
			my %fd = (); # initiate fd, format-depristo
			my @t = split("\t", $_);
   			if (/^\#CHROM/) { # header
				for (my $i=9; $i<= $#t; $i++) {
					my $sample_name = $t[$i];
					my $geno_match = 0;
					for (my $j=0; $j<=$#taxa; $j++) {
						if ($sample_name eq $taxa[$j]) {
							$inputgeno{$i} = $geno[$j];
							$geno_match = 1;
						}
					}
					# if not match, add the col to keep
					if (!$geno_match) {
						push(@tokeep_geno_cols, $i);
					}
				}
				
				if ($opts{r}) { # removed data of specified taxa
					print OUT join("\t", @t[0..8]);
					for (my $m=0; $m<= $#tokeep_geno_cols; $m++) {
						print OUT "\t$t[$tokeep_geno_cols[$m]]";
					}
				} else {
					print OUT $_;
				}
				print OUT "\n";
				next;
			}
			my $discarded = 0;
			my @format = split(/:/, $t[8]); # format column
			my @genocols = keys %inputgeno;
			foreach my $k (@genocols) {
				my @depristo = split(/:/, $t[$k]);
				for (my $i=0; $i<=$#format; $i++) {
					if (!exists $depristo[$i]) {
						$depristo[$i] = "NA";
					}
					$fd{$format[$i]} = $depristo[$i];
				}
			
				if (exists $fd{GT}) {
					if ($fd{GT} ne $inputgeno{$k}) {
						$discarded = 1;
					}
				} else {
					print STDERR "ERROR: No GT data";
					exit;
				}
			}	
		
			# output
			if (!$discarded) {
				if ($opts{r}) { # removed data of specified taxa
					print OUT join("\t", @t[0..8]);
					for (my $m=0; $m<=$#tokeep_geno_cols; $m++) {
						print OUT "\t$t[$tokeep_geno_cols[$m]]";
					}
				} else {
					print OUT $_;
				}
				print OUT "\n";
			}
		} else {
			print OUT $_;
		}
	}
	close IN;
	close OUT;
}


########################################################################
# recall
#########################################################################
sub recall {
	my %opts = ();
	&getopts('A:a:B:b:s:m:o:', \%opts);
	my $homoMinAD = (defined $opts{A}) ? $opts{A} : 4;
	my $homoMinADperc = (defined $opts{a}) ? $opts{a} : 0.9;
	my $heteroMinAD = (defined $opts{B}) ? $opts{B} : 1;
	my $heteroMinADperc = (defined $opts{b}) ? $opts{b} : 0.2;
	my $maxMissing = (defined $opts{m}) ? $opts{m} : 1;
	my $alleleseparator = (defined $opts{s}) ? $opts{s} : "/";
	my $missing_geno = ".".$alleleseparator.".";
	if ($alleleseparator ne "/" and $alleleseparator ne "|") {
		print STDERR "ERROR: only | or / are allowed for --alleleseparator\n";
	}

	die(qq/
	Usage: $0 recall [options] <vcf> 
	{arguments]
	-A: minimum counts of an allele for homozygous genotype calls (integer >0) (4)
	-a: minimum percentage of counts of an allele out of read depth of the site for homozygous genotype calls (float 0-1) (0.9)
	-B: minimum counts of each allele for heterozygous genotype calls (integer >0) (1)
	-b: minimum percentage of counts of each allele out of read depth of the site for heterozygous genotype calls (float 0-1) (0.2)
	-s: alleleseparator, either | or \/; ("\/")
	-m: maximum missing rate after genotype recalling (float 0-1); (1)	
	-o: output vcf file

	examples:
	DH and advanced RIL population:
		-A 2 -a 0.8 -B 1 -b 0.2
	RIL population:
		-A 3 -a 0.9 -B 1 -b 0.2
	F2 population:
		-A 4 -a 0.9 -B 1 -b 0.1
	[note]: this script only works for biallele vcf results.
	\n/) if (@ARGV==0 && -t STDIN);

	my $vcf = $ARGV[0];
	my $output;
	if (!defined $opts{o}) {
		$output = $vcf;
		$output =~ s/.vcf//g;
		$output .= ".recall.vcf";
	} else {
		$output = $opts{o};
	}
	# output
	open(OUT, ">", $output) || die;

	# print out parameters:
	print STDERR "Parameters:\n";
	print STDERR "input vcf file: $vcf\n";
	print STDERR "min homo depth: $homoMinAD\n";
	print STDERR "min homo depth%: $homoMinADperc\n";
	print STDERR "min hetero depth: $heteroMinAD\n";
	print STDERR "min hetero depth%: $heteroMinADperc\n";
	print STDERR "allele separator: $alleleseparator\n";
	print STDERR "max missing rate per site: $maxMissing\n\n";
	# 1	14	.	C	A	5539.93	.	AC=89;AF=0.420;AN=212;BaseQRankSum=-3.185;DP=318;ExcessHet=-0.0000;FS=5.002;InbreedingCoeff=0.4833;MLEAC=137;MLEAF=0.646;MQ=57.87;MQRankSum=-7.644;QD=25.63;ReadPosRankSum=-0.033;SOR=0.445	GT:AD:DP:GQ:PL	1|1:0,1:1:3:45,3,0
	
	# open the file:
	open(IN, "<", $vcf) || die;
	while (<IN>) {
		if (!/\#\#/) {
			chomp;
			my %fd = (); # initiate fd, format-depristo
			my @t = split;
  	 		if (/^\#CHROM/) {
				print OUT "$_\n";
				next; chomp;
			}
			my @format = split(/:/, $t[8]); # format column	
			my $genodiff_count = 0;
			my $missing_count = 0;
			my @newgeno = ();
			print STDERR "===$t[0]-$t[1]===\n";
			for (my $k=9; $k<=$#t; $k++) {
				my @depristo = split(/:/, $t[$k]);
				for (my $i=0; $i<=$#format; $i++) {
					if (!exists $depristo[$i]) {
						$depristo[$i] = "NA";
					}
					$fd{$format[$i]} = $depristo[$i];
				}
			
				# original GT data
				my $originalgeno;
				my $newgeno = $missing_geno;
				if (exists $fd{GT}) {
					$originalgeno = $fd{GT};
				} else {
					print STDERR "ERROR: No GT data";
					exit;
				}
			
				# DP
				my $dp = -1; # initiate DP value
				if (exists $fd{DP}) {
					$dp = $fd{DP};
				}
				# AD
				if (exists $fd{AD}) {
					my $ad = $fd{AD};

					if ($ad ne "NA") {
						my @ad = split(/,/, $ad);
						if ($dp == -1) {
							$dp = $ad[0] + $ad[1]; # sum of counts of two allels
						}
					
						# recall
						$newgeno = &genocall($missing_geno, $alleleseparator, $ad[0], $ad[1], $dp, $homoMinAD, $homoMinADperc, $heteroMinAD, $heteroMinADperc);
					
						# diff?
						my $genodiff = &genocompare($originalgeno, $newgeno);
						$genodiff_count += $genodiff;
						if ($genodiff == 1) {
							print STDERR "   $originalgeno to $newgeno - $t[$k]\n";
						}		
					}
				} else {
					print STDERR "ERROR: No AD data";
					exit;
				}
				
				push(@newgeno, $newgeno);
				$missing_count++ if ($newgeno eq $missing_geno);
			} # end for
			if ($missing_count / ($#t - 8) <= $maxMissing) {
				print OUT join("\t", @t[0..7]);
				print OUT "\tGT\t"; # new format column
				print OUT join("\t", @newgeno);
				print OUT "\n";
			}
			print STDERR "\* $genodiff_count genotypes were modified\n";
		} else {
			print OUT $_;
		}
	}
	close IN;

	sub genocall {
	# recall genotype based on counts of alleles
		#my $geno = $missing_geno;
		my ($geno, $sepchar, $ad1, $ad2, $tdp, $homoAD, $homoADp, $heteroAD, $heteroADp) = @_;
		# homo 1
		if ($ad1 >= $homoAD and ($ad1 / $tdp) >= $homoADp) {
			$geno = "0".$sepchar."0";
		}
		# homo 2
		if ($ad2 >= $homoAD and ($ad2 / $tdp) >= $homoADp) {
			$geno = "1".$sepchar."1";
		}
		# hetero
		if ($ad1 >= $heteroAD and $ad2 >= $heteroAD and
		    ($ad1 / $tdp) >= $heteroADp and ($ad2 / $tdp) >= $heteroADp) {
			$geno = "0".$sepchar."1";
		}
		return $geno;
	}

	sub genocompare {
	# compare two genotypes, return 0 or 1
		my $diff = 0;
		my ($g1, $g2) = @_;
		my @g1 = split(/\||\//, $g1);
		my @g2 = split(/\||\//, $g2);
		my @sortg1 = sort {$a cmp $b} @g1;
		my @sortg2 = sort {$a cmp $b} @g2;
		my $sortg1 = join("", @sortg1);
		my $sortg2 = join("", @sortg2);
		if ($sortg1 ne $sortg2) {
			$diff = 1;
		}
		return $diff;
	}
}


#########################################################################
# recode
#########################################################################
sub recode {
	my %opts = (r=>1, a=>2, h=>3, m=>0);
	getopts('r:a:h:m:', \%opts);
	my $refhomo = $opts{r};
	my $althomo = $opts{a};
	my $hetero = $opts{h};
	my $missing = $opts{m};

	my %genocode = ("0/0" => $refhomo, "0|0" => $refhomo, "0" => $refhomo,
	                "1/1" => $althomo, "1|1" => $althomo, "1" => $althomo,
					"0/1" => $hetero, "0|1" => $hetero,
					"1/0" => $hetero, "1|0" => $hetero,
					"./." => $missing, ".|." => $missing, "." => $missing,
					"2" => $missing, "3" => $missing);

	die(qq/
Usage: vcfbox.pl recode [options] <vcf> 
[arguments]
  -r: code for homozygous genotypes of ref alleles; default=1
  -a: code for homozygous genotypes of alt alleles; default=2
  -h: code for heterozygous genotypes; default=3
  -m: code for missing data; default=0
\n/) if (@ARGV==0 && -t STDIN);

	my $vcf = $ARGV[0];
	open(IN, "<", $vcf) || die;
	while(<IN>) {
		chomp;
		if (! /^\#\#/) {
			my @line = split;
			if (/^\#/) {
				$line[0] =~ s/^\#//;
				print join("\t", @line[0,1,3,4,9..$#line]);
				print "\n";
			} else {
				print join("\t", @line[0,1,3,4]);
				
				for (my $i=9; $i<=$#line; $i++) {
					my $original_geno = $line[$i];
					$original_geno =~ s/\:.*//g;
					if (exists $genocode{$original_geno}) {
						print "\t$genocode{$original_geno}";
					} else {
						print STDERR "ERROR: $line[$i] is NOT a standard genotype format.\n";
						exit;
					}
				}
				
				print "\n";
			}
		}
	}
	close IN;
}


######################################################################### 
# sitefilt
######################################################################### 
sub sitefilt {
	my %opts = ();
	&getopts('v:m:P:p:H:h:o:', \%opts);
	my $minValid = (defined $opts{v}) ? $opts{v} : 1;
	if ($minValid < 1) {
		print STDERR "-v must be >0\n";
		exit;
	}
	my $maxMissingPerc = (defined $opts{m}) ? $opts{m} : 1;
	my $homoPercMin = (defined $opts{p}) ? $opts{p} : 0;
	my $homoPercMax = (defined $opts{P}) ? $opts{P} : 1;
	my $heteroPercMin = (defined $opts{h}) ? $opts{h} : 0;
	my $heteroPercMax = (defined $opts{H}) ? $opts{H} : 1;


	die(qq/
	Usage: $0 sitefilt [options] <vcf> 
	{arguments]
		-v: minimum datapoints\/sites with non-missing genotype (integer, >0); default=1
		-m: maximum missing rate allowed (float 0-1) (1)
		-p: minimum homozygous rate (float 0-1); (0)
		-P: maximum homozygous rate (float 0-1); (1)
		-h: minimum heterozygous rate (float 0-1); (0)
		-H: maximum heterozygous rate (float 0-1); (1)
		-o: output vcf file
	[note]: this script only works for biallele vcf results.
	\n/) if (@ARGV==0 && -t STDIN);

	my $vcf = $ARGV[0];
	my $output;
	if (!defined $opts{o}) {
		$output = $vcf;
		$output =~ s/.vcf//g;
		$output .= ".sitefilt.vcf";
	} else {
		$output = $opts{o};
	}
	
	open(OUT, ">", $output) || die;

	# print out parameters:
	print STDERR "Parameters:\n";
	print STDERR "- input vcf file: $vcf\n";
	print STDERR "- minimum valid data points: $minValid\n";
	print STDERR "- max missing rate per site: $maxMissingPerc\n";
	print STDERR "- minimum homozygous rate: $homoPercMin\n";
	print STDERR "- maximum homozygous rate: $homoPercMax\n";
	print STDERR "- minimum heterozygous rate: $heteroPercMin\n";
	print STDERR "- maximum heterozygous rate: $heteroPercMax\n";
	print STDERR "\n";

	# 1	146	.	C	A	5539.93	.	AC=89;AF=0.420;AN=212;BaseQRankSum=-3.185;DP=318;ExcessHet=-0.0000;FS=5.002;InbreedingCoeff=0.4833;MLEAC=137;MLEAF=0.646;MQ=57.87;MQRankSum=-7.644;QD=25.63;ReadPosRankSum=-0.033;SOR=0.445	GT:AD:DP:GQ:PL	1|1:0,1:1:3:45,3,0
	my @samples; # sample array
	my ($homo1, $homo2, $hetero, $missing);

	open(IN, $vcf) || die;
	while (<IN>) {
		if (!/\#\#/) {
			$homo1 = 0;
			$homo2 = 0;
			$hetero = 0;
			$missing = 0;
			my %fd = (); # initiate fd, format-depristo
			chomp;
			my @t = split;
   			if (/^\#CHROM/) {
				for (my $i=9; $i<=$#t; $i++) {
					push(@samples, $t[$i]);
				}
				next; chomp;
			}
		
			# format in column 9
			my @format = split(/:/, $t[8]); # format column	
		
			# for each genotyped sample
			for (my $k=9; $k<=$#t; $k++) {
				my @depristo = split(/:/, $t[$k]);
				for (my $i=0; $i<=$#format; $i++) {
					if (!exists $depristo[$i]) {
						$depristo[$i] = "NA";
					}
					$fd{$format[$i]} = $depristo[$i];
				}
			
				# GT data
				my $geno;
				if (exists $fd{GT}) {
					$geno = $fd{GT};
				} else {
					$geno = $depristo[0];
					print STDERR "ERROR: No GT indicator in column 9.\n";
				}
			
				my $cursample = $samples[$k - 9];
				$geno =~ s/\||\///g;

				$homo1++ if ($geno eq "00");
				$homo2++ if ($geno eq "11");
				$hetero++ if ($geno eq "01" or $geno eq "10");
				$missing++ if ($geno eq "..");
			} # end for
			
			my $ntaxa = $#t - 8;
			my $nvalid = $homo1 + $homo2 + $hetero;
			
			if ($nvalid >= $minValid) {
				my $missing_perc = $missing / $ntaxa;
				my $homo1_perc = $homo1 / $nvalid;
				my $homo2_perc = $homo2 / $nvalid;
				my $hetero_perc = $hetero / $nvalid;
				if ($missing_perc <= $maxMissingPerc and
				    $homo1_perc >= $homoPercMin and $homo1_perc <= $homoPercMax and
					$homo2_perc >= $homoPercMin and $homo2_perc <= $homoPercMax and
				    $hetero_perc >= $heteroPercMin and $hetero_perc <= $heteroPercMax) {
				
					print OUT "$_\n";
				}

			}
		}
	}
	close IN;
	close OUT;
}


######################################################################### 
# taxafilt
######################################################################### 
sub taxafilt {
	#my %opts = (v=>1000, m=>1, f=>0, c=>1, o=>"test");
	my %opts = ();
	&getopts('v:m:f:c:o:', \%opts);
	my $minValid = (defined $opts{v}) ? $opts{v} : 1000;
	my $maxMissingPerc = (defined $opts{m}) ? $opts{m} : 1;
	my $heteroPercMin = (defined $opts{f}) ? $opts{f} : 0;
	my $heteroPercMax = (defined $opts{c}) ? $opts{c} : 1;

	die(qq/
	Usage: $0 taxafilt [options] <vcf> 
	{arguments]
		-v: minimum datapoints\/sites with non-missing genotype (integer); default=1000
		-m: maximum missing rate allowed (float 0-1) (1)
		-f: minimum heterozygous rate (float 0-1); (0)
		-c: maximum heterozygous rate (float 0-1); (1)
		-o: output vcf file
	[note]: this script only works for biallele vcf results.
	\n/) if (@ARGV==0 && -t STDIN);

	my $vcf = $ARGV[0];
	my $output;
	if (!defined $opts{o}) {
		$output = $vcf;
		$output =~ s/.vcf//g;
		$output .= ".taxafilt.vcf";
	} else {
		$output = $opts{o};
	}

	# print out parameters:
	print STDERR "Parameters:\n";
	print STDERR "- input vcf file: $vcf\n";
	print STDERR "- minimum valid data points: $minValid\n";
	print STDERR "- max missing rate per site: $maxMissingPerc\n";
	print STDERR "- minimum heterozygous rate: $heteroPercMin\n";
	print STDERR "- maximum heterozygous rate: $heteroPercMax\n";
	print STDERR "\n";

	# 1	146	.	C	A	5539.93	.	AC=89;AF=0.420;AN=212;BaseQRankSum=-3.185;DP=318;ExcessHet=-0.0000;FS=5.002;InbreedingCoeff=0.4833;MLEAC=137;MLEAF=0.646;MQ=57.87;MQRankSum=-7.644;QD=25.63;ReadPosRankSum=-0.033;SOR=0.445	GT:AD:DP:GQ:PL	1|1:0,1:1:3:45,3,0
	my @samples; # sample array
	my (%homo1, %homo2, %hetero, %missing);

	print STDERR "taxon\tNum_homo1\tNum_homo2\tNum_hetero\tNum_missing\tPerc_missing\tFilter\n";

	open(IN, $vcf) || die;
	while (<IN>) {
		if (!/\#\#/) {
			my %fd = (); # initiate fd, format-depristo
			chomp;
			my @t = split;
   			if (/^\#CHROM/) {
				for (my $i=9; $i<=$#t; $i++) {
					push(@samples, $t[$i]);
				}
				next; chomp;
			}
		
			# format in column 9
			my @format = split(/:/, $t[8]); # format column	
		
			# for each genotyped sample
		
			for (my $k=9; $k<=$#t; $k++) {
				my @depristo = split(/:/, $t[$k]);
				for (my $i=0; $i<=$#format; $i++) {
					if (!exists $depristo[$i]) {
						$depristo[$i] = "NA";
					}
					$fd{$format[$i]} = $depristo[$i];
				}
			
				# GT data
				my $geno;
				if (exists $fd{GT}) {
					$geno = $fd{GT};
				} else {
					$geno = $depristo[0];
					print STDERR "WARNING: No GT indicator in column 9.\n";
				}
			
				my $cursample = $samples[$k - 9];
				$geno =~ s/\||\///g;

				$homo1{$cursample}++ if ($geno eq "00");
				$homo2{$cursample}++ if ($geno eq "11");
				$hetero{$cursample}++ if ($geno eq "01" or $geno eq "10");
				$missing{$cursample}++ if ($geno eq "..");
			} # end for
		}
	}
	close IN;

	my $npass = 0;
	my $cmd = sprintf("%s%s%s", "cut ", $vcf, " -f 1-9");
	for (my $i=0; $i<=$#samples; $i++) {
		my $taxon_filter = "discarded";
		my $homo1c = 0;
		my $homo2c = 0;
		my $heteroc = 0;
		my $missingc = 0;
		$homo1c = $homo1{$samples[$i]} if (exists $homo1{$samples[$i]});
		$homo2c = $homo2{$samples[$i]} if (exists $homo2{$samples[$i]});
		$heteroc = $hetero{$samples[$i]} if (exists $hetero{$samples[$i]});
		$missingc = $missing{$samples[$i]} if (exists $missing{$samples[$i]});
		my $total = $homo1c + $homo2c + $heteroc + $missingc;
		my $nvalid = $homo1c + $homo2c + $heteroc;
		my $missing_perc = 1;
		my $hetero_perc = "NA";
		if ($nvalid > 0) {
			my $missing_perc_raw = $missingc / $total;
			$hetero_perc = $heteroc / $nvalid;
			$missing_perc = sprintf("%.3f", $missing_perc_raw);	
			# judge using input criteria
			if ($nvalid >= $minValid and $missing_perc <= $maxMissingPerc and
			    $hetero_perc >= $heteroPercMin and $hetero_perc <= $heteroPercMax) {
				$taxon_filter = "pass";
				$npass++;
				my $selected_col = $i + 10;
				$cmd .= ",";
				$cmd .= $selected_col;
			}
		}
		# report
		print STDERR "$samples[$i]\t$homo1c\t$homo2c\t$heteroc\t$missingc\t$missing_perc\t$taxon_filter\n";
	}

	$cmd .= " > ";
	$cmd .= $output;
	# execute
	print STDERR "$npass samples passed criteria\n";
	print STDERR "$cmd\n";
	system($cmd);
}

#########################################################################
# compare
#########################################################################
sub compare {
	my %opts = ();
	&getopts('a:b:', \%opts);

	if (!defined $opts{a} or (!defined $opts{b})) {
		print STDERR "both -a and -b are required\n";
		exit;
	}
	
	my $taxon_1 = $opts{a};
	my $taxon_2 = $opts{b};
	
	die(qq/
	Usage: $0 compare [options] <vcf> 
	{arguments]
	-a: exact name of taxon 1, required
	-b: exact name of taxon 1, required
	
	[note]: this script only works for biallele vcf results.
	\n/) if (@ARGV==0 && -t STDIN);

	my $vcf = $ARGV[0];
	
	# print out parameters:
	print STDERR "Parameters:\n";
	print STDERR "input vcf file: $vcf\n";
	print STDERR "taxon 1: $taxon_1\n";
	print STDERR "taxon 2: $taxon_2\n";
	
	# open the file:
	my $taxoncol_1 = 0;
	my $taxoncol_2 = 0;

	my $total_sites = 0;
	my $both_missing = 0;
	my $taxon_1_missing = 0;
	my $taxon_2_missing = 0;
	my ($equal, $diff, %genocomb);

	open(IN, "<", $vcf) || die;
	while (<IN>) {
		if (!/^\#\#/) {
			chomp;
			my %fd = (); # initiate fd, format-depristo
			my @t = split("\t", $_);
   			if (/^\#CHROM/) { # header
				# determine which cols containing taxon_1 and _2:
				for (my $i=9; $i<= $#t; $i++) {
					my $sample_name = $t[$i];
					if ($sample_name eq $taxon_1) {
						$taxoncol_1 = $i;
					} elsif ($sample_name eq $taxon_2) {
						$taxoncol_2 = $i;
					}
				}
				# if not match, add the col to keep
				if ($taxoncol_1 == 0) {
					print STDERR "ERROR: $taxon_1 did not match any taxa\n";
				}

				if ($taxoncol_2 == 0) {
					print STDERR "ERROR: $taxon_2 did not match any taxa\n";
				}
				
				if ($taxoncol_1 == 0 or $taxoncol_2 == 0) {
					exit;
				}
			} else {
				$total_sites++;
				my @format = split(/:/, $t[8]); # format column	
				my $gt_pos = undef; # position of GT data in depristo
				for (my $i=0; $i<=$#format; $i++) {
					if ($format[$i] eq "GT") {
						$gt_pos = $i;
					}
				}
				
				my @depristo_1 = split(/:/, $t[$taxoncol_1]);
				my @depristo_2 = split(/:/, $t[$taxoncol_2]);
				
				my ($geno_1, $geno_2);
				if (defined $gt_pos) {
					if (exists$depristo_1[$gt_pos]) {
						$geno_1 = $depristo_1[$gt_pos];
					} else {
						print STDERR "GT infomraton missed for $taxon_1\n";
						print STDERR "-- $t[$taxoncol_1]\n";
					}

					if (exists$depristo_1[$gt_pos]) {
						$geno_2 = $depristo_2[$gt_pos];
					} else {
						print STDERR "GT infomraton missed for $taxon_2\n";
						print STDERR "-- $t[$taxoncol_2]\n";
					}
				} else {
					print STDERR "ERROR: GT information in missing at column 9\n";
					exit;
				}
				
				# check if missing
				if ($geno_1 =~ /\./ and ($geno_2 =~ /\./)) {
					$both_missing++;
				} elsif ($geno_1 =~ /\./) {
					$taxon_1_missing++;
				} elsif ($geno_2 =~ /\./) {
					$taxon_2_missing++;
				} else { # if not missing, ...
					$genocomb{$geno_1."\t".$geno_2}++;
					if ($geno_1 eq $geno_2) {
						$equal++;
					} else {
						$diff++;
					}
				}
			}
		}
	}
	close IN;
	
	if ($total_sites == 0) {
		print "no variant entries\n";
		exit;
	}


	# summary
	my $total_valid = $equal + $diff;
	my $equal_perc = "NA";
	my $diff_perc = "NA";
	if ($total_valid > 0) {
		$equal_perc = $equal / $total_valid * 100;
		$diff_perc = $diff / $total_valid * 100;
	}
	print "#Genotype comparison of $taxon_1 and $taxon_2\n";
	print "Category\tCount\tPercentage\n";
	printf ("%s\t%d\t%.2f\n", "Both missing", $both_missing, $both_missing / $total_sites * 100);
	printf ("%s %s\t%d\t%.2f\n", $taxon_1, "missing", $taxon_1_missing, $taxon_1_missing / $total_sites * 100);
	printf ("%s %s\t%d\t%.2f\n", $taxon_2, "missing", $taxon_2_missing, $taxon_2_missing / $total_sites * 100);
	printf ("%s\t%d\t%.2f\n", "Both nonmissing", $total_valid, $total_valid / $total_sites * 100);
	printf("%s\t%d\t%.2f\n", "equal", $equal, $equal_perc);
	printf("%s\t%d\t%.2f\n", "diff", $diff, $diff_perc);

	# details:
	print "\#Details of the comparion of valid data\n";
	print "\#$taxon_1\t$taxon_2\tCount\n";
	foreach my $comb (sort {$a cmp $b} keys %genocomb) {
		print "\#$comb\t$genocomb{$comb}\n";
	}

}

#########################################################################
# select
#########################################################################
sub select {
	my %opts = ();
	&getopts('f:l:ia', \%opts);

	die(qq/
	Usage: $0 select [options] <vcf> 
	{arguments]
		-f: a file containing taxon\/sample names; one per line
		-l: taxon\/sample names separated by comma; names need to match names in the vcf file
		-i: ignore taxon redundancy in the vcf input if specified (off by default)
		-a: allow redundancy of input taxa or taxa not in the vcf file
		only either -f or -l can be input but at least one is required
	\n/) if (@ARGV==0 && -t STDIN);

	my $vcf = $ARGV[0];
	my @taxa = ();

	if (defined $opts{f} and defined $opts{l}) {
		print STDERR "ERROR: only either -f or -l can be input.\n";
		exit;
	}
	
	# if input is a file
	if (defined $opts{f}) {
		open(IN, $opts{f}) || die;
		while(<IN>) {
			chomp;
			push(@taxa, $_);
		}
	} elsif (defined $opts{l}) {
		@taxa = split(/,/, $opts{l});
	}

	# check redundancy
	if ($#taxa==0) {
		print STDERR "ERROR: no taxa was input\n";
	}

	my %taxa;
	foreach (@taxa) {
		$taxa{$_}++;
		if ($taxa{$_} > 1) {
			print STDERR "WARNING: redundancy was identified in input taxa\n";
			if (!$opts{a}) {
				exit;
			}
		}
	}
	
	@taxa = keys %taxa; # ensure no redundancy

	# open the vcf file:
	my (@tokeep_geno_cols, @tokeep_taxa);
	my %tokeep_taxa;

	open(IN, $vcf) || die;
	while (<IN>) {
		if (!/\#\#/) {
			chomp;
			my @t = split("\t", $_);
   			if (/^\#CHROM/) { # header
				for (my $i=9; $i<= $#t; $i++) {
					my $sample_name = $t[$i];
					if (exists $taxa{$sample_name}) {
						push(@tokeep_taxa, $sample_name);
						push(@tokeep_geno_cols, $i);
					}
				}
				
				# check redundancy of identified taxa from vcf
				foreach (@tokeep_taxa) {
					$tokeep_taxa{$_}++;
					if ($tokeep_taxa{$_}>1) {
						print STDERR "WARNING: redundancy was identified in taxa of $vcf\n";
						if (!$opts{i}) {
							exit;
						}
					}
				}
				
				# check if all selected taxa were extracted 
				if ($#taxa != $#tokeep_geno_cols) {
					foreach (@taxa) {
						if (!exists $tokeep_taxa{$_}) {
							print STDERR "WARNING: $_ was not found in $vcf\n";
						}
					}
					exit if (!$opts{a}); # quit if not allow missed taxa
				}
			}	
			# output
			print join("\t", @t[0..8]);
			print "\t";
			print join("\t", @t[@tokeep_geno_cols]);
			print "\n";
		} else {
			print $_;
		}
	}
	close IN;
}

#########################################################################
# allele
#########################################################################
sub allele {
	my %opts = ();
	&getopts('m:o:', \%opts);
	
	my $missing_str = (defined $opts{m}) ? $opts{m} : 0;

	die(qq/
	Usage: $0 allele [options] <vcf> 
	[arguments]
	-m: string to represent missing data (0) 
	-o: output filename
	\n/) if (@ARGV==0 && -t STDIN);

	my $vcf = $ARGV[0];
	my $output;
	if (!defined $opts{o}) {
		$output = $vcf;
		$output =~ s/.vcf//g;
		$output .= ".allele.txt";
	} else {
		$output = $opts{o};
	}

	open(OUT, ">", $output) || die;

	# open the vcf file:
	open(IN, "<", $vcf) || die;
	while (<IN>) {
		if (!/\#\#/) {
			my %fd = (); # initiate fd, format-depristo
			my @t = split;
   			if (/^\#CHROM/) {
				print OUT "CHR\t";
				print OUT join("\t", @t[1,3,4,5]);
				for (my $i=9; $i<= $#t; $i++) {
					printf OUT "\t%s%s\t%s%s", $t[$i], "_REF", $t[$i], "_ALT";
				}
				print OUT "\n";
				next;
			}
			next if ($t[4] eq '.'); # skip non-var sites
 	  		next if ($t[3] eq 'N'); # skip sites with unknown ref ('N')
		
			print OUT join("\t", @t[0, 1, 3, 4, 5]);

			my @format = split(/:/, $t[8]); # format column
			for (my $k = 9; $k <= $#t; $k++) {
				my @depristo = split(/:/, $t[$k]);
				for (my $i=0; $i<=$#format; $i++) {
					if (!exists $depristo[$i]) {
						$depristo[$i] = "NA";
					}
					$fd{$format[$i]} = $depristo[$i];
				}
				if (exists $fd{AD}) {
					my @ad = split(/,/, $fd{AD}); # allele depth
					my $refc = $missing_str;
					my $altc = $missing_str;
					if ($#ad == 1) {
						$refc = $ad[0]; # ref count
						$altc = $ad[1]; # alt count
					}
					print OUT "\t$refc\t$altc";
				} else {
					print STDERR "ERROR:";
					print STDERR "$_\n";
					print STDERR "No AD data, AD=allele depth\n";
					exit;
				}
			}
			print OUT "\n";
		}
	}
	close IN;
}

