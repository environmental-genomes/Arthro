use strict;
use warnings;

die "\n usage: perl $0 assembly2tax.txt tree\n" unless (scalar(@ARGV)>1);

my $assembly2tax=shift;
my $tree=shift;

open my $FH, $tree or die "$tree:$!\n";
my $tre=<$FH>;

open my $A, $assembly2tax or die "$assembly2tax: $!\n";
while(<$A>){
chomp;
next if(/^#/);
next if(/^$/);
my @t=split /\t/;
$t[1]=~s/\:.*//;
$tre=~s/$t[0]/$t[1]/;
}

print $tre;