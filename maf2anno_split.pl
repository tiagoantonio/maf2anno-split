###############################################################################
#
# April 26 2016 
#
# maf2anno.split
#
# This Perl script splits a TGCA .maf file into separated sample files and 
# annotate each of them using ANNOVAR (hg19) instaled at /share/apps/ANNOVAR.
# Outputs: ANNOVAR output files and a log file containing the total number of 
# samples annotated and total number of SNVs annotated.
#
# 
# USAGE: perl maf2anno.split <TGCAfile.maf>
# 
# Required modules: IPC::System::Simple, Parallel::ForkManager, List::MoreUtils, 
#                   List::Util
#
# If you have any question or suggestions please contact the author:
#
# Author: Tiago A. de Souza
# tiagoantonio@gmail.com
# github.com/tiagoantonio
#
# Licensed under MIT License
#
#################################################################################

#! /usr/bin/perl
use IPC::System::Simple qw(system capture);
use Parallel::ForkManager;
use List::MoreUtils qw(uniq);
use List::Util qw(sum);
use strict;
use warnings;

# ANNOVAR INSTALL DIR
my $annovar_dir;
$annovar_dir="\/share\/apps\/annovar";

# variables
my $inputRawSNV; # loading files and parameters
my @rawFile; # loading files and parameters
my $rawLine; my (@chr, @pos, @snp, @ref, @alt, @sample, @ID); # reading maf file
my (@groupseq, @group); # spliting sample names
my $snp; # creating an ANNOVAR input file for each sample
my @uniq_sample; my $uniq_sample; # saving unique sample names
my $pm; # Using a max of 8 threads to call Annovar 
my $pid; my @ARGS; # Annotation with Annovar installed at /share/apps/annovar
my (@snplineannovar, @countSNV, @i2); my $countfile; # opening ANNOVAR log files
my $sumSNV;# counting the total number of SNVs annotate
my @i; my$i; # index

unless (@ARGV){
	die "\nUSAGE : perl maf2anno.split <TGCAfile.maf>\n";
}

# loading files and parameters
$inputRawSNV = $ARGV[0]; #<tabular_maf_file>
open (inputRawSNV, $inputRawSNV);
@rawFile=<inputRawSNV>;
unless (@rawFile){
	die "\nERROR : Could not load <TGCAfile.maf>\n";
}

# reading maf file
foreach $rawLine (@rawFile){
	@i = split (/\t/, $rawLine);
	chomp (@i);
	push (@chr, "$i[4]");
	push (@pos, "$i[5]");
	push (@snp, "$i[9]");
	push (@ref, "$i[10]");
	push (@alt, "$i[12]");
	push (@sample, "$i[15]");
	push (@ID, "$i[35]");
}

# spliting sample names
@groupseq = split (/\_/, "$ID[1]");
@group = split (/\./, $groupseq[1]);

# creating output directory
mkdir("$group[0]", 0755);

# creating an ANNOVAR input file for each sample
for $i (0..$#snp){
	if ($snp[$i] eq "SNP"){
		open (SAMPLE, ">>$group[0]/$sample[$i].txt");
		print SAMPLE "chr$chr[$i]\t$pos[$i]\t$pos[$i]\t$ref[$i]\t$alt[$i]\n";
		close (">>$group[0]/$sample[$i].txt");
	}
}

# accessing sample names
open (INPUTTABLE, ">>$group[0]/input-table-$group[0].txt");

# saving unique sample names
@uniq_sample = uniq @sample;

# creating input table for woland
for $i (1..$#uniq_sample){
	print INPUTTABLE "$group[0]\t$uniq_sample[$i].txt.variant_function\n";
}
close ("$group[0]/input-table-$group[0].txt");

# Using a max of 8 threads to call Annovar 
$pm = Parallel::ForkManager->new(8);

# Annotation with Annovar installed at /share/apps/annovar
for $i (1..$#uniq_sample){
	@ARGS=("-geneanno", "-buildver", "hg19", "$group[0]/$uniq_sample[$i].txt", "$annovar_dir/humandb");
	$pid=$pm->start and next;
	system ($^X, "$annovar_dir/annotate_variation.pl", @ARGS);
	$pm->finish;
	@ARGS=();
}
# wait all children processes to finish
$pm->wait_all_children;

# opening ANNOVAR log files
for $i (1..$#uniq_sample){
	open (ANNOVARLOG, "<$group[0]/$uniq_sample[$i].txt.log");
	@snplineannovar=<ANNOVARLOG>;
	@i2 = split (/\ /, $snplineannovar[11]);
	push (@countSNV, $i2[5]);
	++$countfile;
}

# counting the total number of SNVs annotated
$sumSNV = sum(@countSNV);

# creating an unique log file cointaing the number of samples annotated and total number of SNVs annotated.
open (LOG, ">>$group[0]/input-table-$group[0].txt.log");
print LOG "Number of samples annotated: $countfile\n";
print LOG "Total number of SNPs annotated: $sumSNV";

exit;