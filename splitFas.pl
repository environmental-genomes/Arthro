#!/usr/bin/env perl

# split fasta file to individual entry per file
# 
#
use strict;
use warnings;
use Getopt::Long;

my $INFILE;
my $mode=0;
my $IDFILE;
GetOptions ("fasta=s" => \$INFILE, 
            "id=s"   => \$IDFILE,      
            "mode=i"  => \$mode)   
      or die("Error in command line arguments\n");

sub usage {
print STDERR <<EOF;
  usage: perl $0 --fasta your_fasta_file
	options:
        --fasta [FILE] a fasta file to process
        --id    [FILE] put your target sequence name in this file
        --mode  [INT]  process mode; possible value: 0, or other positive number
                       if set to 0: each sequence would be ouptuted to individual file;
                       if set to other number, sequence would be outputed to file 
                       based in their first INT  number of character in the sequence name
                       eg. a sequence with name "ABCDEF" would be outputed to ABC.fas if --mode set to 3	
	
EOF
exit;
}
usage if (! defined $INFILE);
if ($mode <0){
	die "  --mode must be set to a number greater than 0\n";
}
open my $FAS,'<',$INFILE or die (" failed to open file: $INFILE\n");

my %target=();
if(defined $IDFILE){
	open my $ID, '<',$IDFILE;
	chomp;
	s/\r$//;
	s/\s+$//;
	next if (/^$/);
	next if (/^\s+$/);
	#s/\s+/_/g;
	my $id=(split /\s+/, $_)[0];
	$target{$id}=1;
}

$/='>';
while(<$FAS>){
    chomp;
    next if(/^$/);
    next if(/^\s+$/);
    my $entry=(split /\n/,$_)[0];
    $entry=~s/\s+/_/g;
    $entry=~s/:+/_/g;
    $entry=~s/\W+/_/g;
	my $outfile= $mode==0 ? $entry : substr($entry,0,$mode);
    open my $OUT,">>",$outfile.".fas";
    if(! defined $OUT){
        warn(" failed to creat file for this entry: $entry\n");
        next;
    }
	next if (defined $IDFILE && ! defined $target{$entry});
    print $OUT '>';
    print $OUT $_;
    close $OUT;
}
close $FAS;