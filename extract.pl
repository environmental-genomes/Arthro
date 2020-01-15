#!/usr/bin/env perl

use strict;
use warnings;
use Getopt::Std;


my %opts=();
getopts('i:f:on', \%opts);

if(! defined ($opts{f}) || !defined ($opts{i})){
	warn "Usage: perl $0 -f in.fas -i id.txt\n";
	warn "  options:\n";
	warn "     -f  input fasta file\n";
	warn "     -i  only ouput sequence whose id found in this file\n";
	warn "     -o  keep sequence order as it in id file specified by -i \n";
	warn "     -n  use whole string rather than first non-blank strings as sequence name \n";
	warn "\n";
	exit(0);
}

my %target=();
open my $ID,'<', $opts{i} or die "failed to open file \"$opts{i}\"\n";
my $num=0;
while(<$ID>){
    chomp;
	s/\r$//g;
	s/\n$//g;
	s/\s+$//g;
    next if (/^\s+$/);
    $target{$_}=++$num;
}
close $ID;

open my $FAS,'<',$opts{f} or die "failed to open file \"$opts{f}\"\n";
my $rand_str='';
if(defined $opts{o}){
	my @str=("A".."Z", "a".."z",0..9);
	$rand_str = join '', map $str[rand @str] , 1 .. 10;
}
$/='>';
while(<$FAS>){
    $_=~s/>$//;
    next if (/^\s+$/);
    next if (/^$/);
    my $id='';
    if (defined $opts{n}) {
	$id=(split /\n/)[0];
	}
	else {
	$id=(split /\s+/)[0];
	}
	$id=~s/\s+$//g;
    if(defined $target{$id} ){
	$_=~s/\s+$//;
	if (defined $opts{o}){
		open my $fh,'>', $rand_str.".".$target{$id}.".fas";
		print $fh ">$_\n";
		close $fh;
	}
	else{
		print ">$_\n";
	}
	}
}

if(defined $opts{o}){
	foreach my $i (1..$num){
		open my $fh2,'<', $rand_str.".".$i.".fas";
		while(<$fh2>){
			print;
		}
		close $fh2;
		unlink $rand_str.".".$i.".fas";
	}
}