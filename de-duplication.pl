use strict;
use warnings;
use Getopt::Long;

sub usage {
print STDERR <<EOF

usage: perl $0 --file infile

    options:
       --seed    these genomes should be always kept. 
                 multiple genomes should be seperated by comma.
       --aai     the column where to find AAI [default: 6]
       --cutoff  AAI cutoff [default: 99.9] 
       --file    read this file

EOF

}



my %opt=();
GetOptions(\%opt, "seed=s", "aai=i", "cutoff=f", "file=s", "help");

if(defined($opt{help})||!defined($opt{file})){
usage;
exit(1);
}

my $column=defined($opt{aai})? $opt{aai} :6;
my $index=$column-1;

my %seed=();

if(defined($opt{seed})){
	%seed=map {$_=>1} (split /,/,$opt{seed});
}
my $cutoff=99.9;
if(defined($opt{cutoff})){
	$cutoff=$opt{cutoff};
}

my %drop=();
my %kept=();

open my $IN, $opt{file} or die "$!";

while(<$IN>){
chomp;
next if(/^Genome/i);
my @t=split /\t/;

my $genomeA=$t[0];
my $genomeB=$t[2];
next if (defined($drop{$genomeA}) && defined($drop{$genomeA}));
my $aai=$t[$index];
$aai=~s/^\s+//;
$aai=~s/\s+$//;
if($aai=~/[\d\.]{2,}/){
	if($aai >=$cutoff){
		if(defined($seed{$genomeA})){
			if(!defined($seed{$genomeB})){
				$drop{$genomeB}=1;
				next;
			}
		}
		else {
			if(defined ($seed{$genomeB})) {
				$drop{$genomeA}=1;
				next;
			}
			else {
				if(!defined($drop{$genomeA}) && !defined($drop{$genomeB})) {
					$drop{$genomeA}=1;
					next;
				}
			}
		}
	}
}
else {
	print STDERR "line $.: the column specified is not like a number -> $_\n";
}
}
close $IN;


open my $IN2, $opt{file} or die "$!";

while(<$IN2>){
	chomp;
	next if(/^Genome/i);
	my @t=split /\t/;

	my $genomeA=$t[0];
	my $genomeB=$t[2];
	next if (defined($drop{$genomeA}) && defined($drop{$genomeB}));
#	print "$_\n";
	$seed{$genomeA}=1 if (!defined $drop{$genomeA});
	$seed{$genomeB}=1 if (!defined $drop{$genomeB});
}

my $genome=join "\n", (keys %seed);
my $genome2=join "\n", (keys %drop);
print "#Duplicated Genomes\n$genome2\n";
print "#Reserved Genomes\n$genome\n";



