#!/usr/bin/perl

#
# add filename to sequence header in fasta
#
use strict;
use warnings;
use File::Basename;

if(scalar(@ARGV)>0){
	foreach my $input (@ARGV){
		my($filename, $dirs, $suffix)=fileparse($input, (".fas",".fasta",".faa"));
		open my $fh,$input or warn "$input: $!\n";
		if($fh){
			while(<$fh>){
				if(/^>/){
					s/^>//;
					print ">$filename","_",$_;
				}
				else {
					print;
				}
			}
			close $fh;
		}
	}
}
