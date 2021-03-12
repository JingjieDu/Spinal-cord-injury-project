#!/usr/bin/perl
use warnings;
use strict;

my $path = "/fs/project/PAS1117/SCI/SCI_ion_all_coassembly/results/03.BBKW_annotations/";
chdir($path) or die "Could not chdir to $path \n";


my %seq_genomes;
my $fasta_file = "01.SCI_ion_all_coassembly.fasta";
open(FASTA, $fasta_file) or die "Could not open $fasta_file \n";
{
   	local $/ = ">";
   	my $number = 0;
   	while(my $record = <FASTA>){
   		chomp $record;
   		my ($defLine, @seqLines) = split /\n/, $record;
   		my @rip0 = split(/ /, $defLine);
   		my $seqid = $rip0[0];			
   		my $sequence = join('',@seqLines);
   		$seq_genomes{$number}=">".$seqid."\n".$sequence;
   		$number++;
    }
}
close FASTA;	

my $total_number_genomes = scalar(keys %seq_genomes);
my $number_fasta_files_split = $total_number_genomes/10000;

for (my $i = 0; $i < $number_fasta_files_split; $i++) {
	my $outfile = "01.SCI_ion_all_coassembly_split_".$i.".fasta";
	open(OUTFILE, '>', $outfile) or die "Could not open $outfile \n";
	my $j = $i*10000;
	for (my $k = $j; (($k < $j+10000) && ($k < $total_number_genomes)); $k++) {
		print OUTFILE $seq_genomes{$k}."\n";
	}
	close OUTFILE;
}
