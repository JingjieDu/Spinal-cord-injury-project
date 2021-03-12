#!/usr/bin/perl
use warnings;
use strict;

my $directory = "<your directory>";
chdir($directory) or die "Could not chdir to $directory \n";

my $infile = "<input>.fna";
my $outfile = "<output>.fna";
open(OUTFILE, '>', $outfile) or die "Could not open $outfile \n";
	{
		local $/ = ">";	
		open(FASTA, $infile) or die "Could not open $infile\n";
		while(my $record = <FASTA>){
			chomp $record;
			my ($defLine, @seqLines) = split /\n/, $record;		
			my $sequence = join('',@seqLines);
			my $length = length($sequence);
			if ($length >= 10000) {
				print OUTFILE ">".$defLine."\n".$sequence."\n";
			}
		}
		close FASTA;
	}
	close OUTFILE;
