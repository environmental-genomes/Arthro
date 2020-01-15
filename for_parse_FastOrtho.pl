#!/usr/bin/env perl

use strict;
use warnings;
use Getopt::Long;
use IO::File;
use Bio::DB::Fasta;
use File::Spec;
use List::MoreUtils qw/ uniq /;;

my %opts;
  GetOptions ("rast=s" => \$opts{rast},    
        "fasta=s"   => \$opts{fasta},      
        "table=s"  => \$opts{table},
	"fastorth=s"   => \$opts{fastorth},
	"outfa=s" =>  \$opts{outfa},
	"annotation=s" => \$opts{anno},
	"mono!" => \$opts{mono},
	"help!" => \$opts{help})   
  or die("Error in command line arguments\n");


sub help {
	print STDERR "\n";
	print STDERR "Usage: perl $0 --rast rast_dir --fasta fas_dir --table out.table --outfa finalseq.fasta\n";
	print STDERR "      options:\n";
	print STDERR "            --rast      directory where RAST annotation results localted\n";
	print STDERR "            --fasta     directory where fasta files located\n";
	print STDERR "            --table     the filename for all ortholog number per strain\n";
	print STDERR "            --mono      keep mono-copy only\n";
	print STDERR "            --fastorth  output file from fastorth\n";
	print STDERR "            --outfa     filename for final fasta\n";
	print STDERR "            --annotation write final annotation to this file\n";
	print STDERR "\n";
	exit(0);
}


if(defined $opts{help} || !defined $opts{rast} || !defined $opts{fasta} || !defined $opts{table}|| !defined $opts{fastorth} || !defined $opts{anno}){
	&help;
}


my $db = Bio::DB::Fasta->new($opts{fasta},'-reindex'=>1);
if(!defined $db){
	die " Failed to open fasta directory. please check -> $opts{fasta}\n";
}

if (! -d $opts{rast}) {
	die " Failed to find the rast directory: $opts{rast}\n";
}

my $fh_table=IO::File->new($opts{table},'w') or die " Failed to create file $opts{table}: $@\n";
my $fh_outfa=IO::File->new($opts{outfa},'w') or die " Failed to create file $opts{outfa}: $@\n";
my $fh_anno=IO::File->new($opts{anno},'w') or die " Failed to create file $opts{anno}: $@\n";

my %result=();
my %select_id=();
my $strain_index=1;
my %counter=();

my $fastorth=IO::File->new($opts{fastorth},'r') or die "$opts{fastorth}: $@\n";

while(<$fastorth>){
	next if(/^#/);
	chomp;
	s/^\s+//;
	s/\s+$//;
	next if(/^$/);
	my ($orth, $gene_list)=split /:/, $_;
	my $orth_id;
	my $num_genes;
	my $num_taxa;
	if($orth=~/(\w+)\s+\((\d+)\s+genes,\s*(\d+)\s+taxa\)/){
		$orth_id=$1;
		$num_genes=$2;
		$num_taxa=$3;
	}
	else {
		warn "$orth\n"; next;}
	$counter{$orth_id}{taxa}=$num_taxa;
	$counter{$orth_id}{genes}=$num_genes;	

	my @genes=split /\)/,$gene_list;
	foreach my $index (0..$#genes) {
		my $gene_str=$genes[$index];
		$gene_str=~s/^\s+//; 
		$gene_str=~s/\s+$//;
		my ($gene_name, $strain)=split /\(/,$gene_str;
		if(!defined $result{$strain}{$orth_id}){
			$result{$strain}{$orth_id}=[$gene_name];
		}
		else{
			push @{$result{$strain}{$orth_id}},$gene_name;
		}
		if($index==0){
			print $fh_anno $orth_id,"\t",$strain;
			my $random_id=int(rand(scalar(@{$result{$strain}{$orth_id}})));
			$select_id{$strain}{$orth_id}=$random_id;
			print $fh_anno "\t",$result{$strain}{$orth_id}->[$random_id];
			my $rast_file=File::Spec->catfile($opts{rast},$strain.'.tsv');
			my $function=query_gene_function_from_file($rast_file,$result{$strain}{$orth_id}->[$random_id]);
			#print "\t",$random_id, "\t", $result{$strain}{$orth_id}->[$random_id], "\n";
			print $fh_anno defined $function? "\t$function": "\tunavailable";
			print $fh_anno "\n";
		}
	}	
}

close $fh_anno;

my @all_orth=keys %counter;
my @all_strain=keys %result;

my %orth_info=();

my $num_shared=0;
my $mono_copy=0;

foreach my $orth (@all_orth){
	if($counter{$orth}{taxa}==scalar(@all_strain)){
		$orth_info{$orth}{shared}=1; $num_shared++;
		if (defined $opts{mono} && $counter{$orth}{taxa}==$counter{$orth}{genes}){
			$orth_info{$orth}{mono}=1; $mono_copy++;
		}
	}
	else {
		$orth_info{$orth}{shared}=0;
	}
}

warn "\n Total ",scalar(@all_orth)," orthlogs; of them $num_shared shared orthologs, and $mono_copy mono-copy orthologs\n\n";

print $fh_table "orth_id\t", join "\t", @all_strain, "\n";
foreach my $orth (sort @all_orth) {
	print $fh_table $orth;
	foreach my $strain (@all_strain) {
		my $num=defined $result{$strain}{$orth} ? scalar(@{$result{$strain}{$orth}}):0;
		print $fh_table "\t$num";
	}
	print $fh_table "\n";
}
close $fh_table;

foreach my $strain (@all_strain) {
	print $fh_outfa '>', $strain,"\n";
	foreach my $orth (sort keys %orth_info) {
		next if ($orth_info{$orth}{shared}==0);
		next if (defined $opts{mono} && !defined $orth_info{$orth}{mono});
		my $seqid= defined $select_id{$strain}{$orth} ? $result{$strain}{$orth}->[$select_id{$strain}{$orth}]: $result{$strain}{$orth}->[int(rand(scalar(@{$result{$strain}{$orth}})))];
		#warn "scalar(@{$result{$strain}{$orth}})==== $counter{$orth}{taxa}=== $counter{$orth}{genes} \n";
		print $fh_outfa query_fasta_by_seqid($db,$seqid);
	}
	print $fh_outfa "\n";
}
close $fh_outfa;



sub query_gene_function_from_file {
	my ($file, $keyword)=@_;
	my $fh=IO::File->new($file,'r');
	if(!defined $fh){
		warn " Failed to read file: $file\n";
		return undef;
	}
	my $res;
	while(<$fh>){
		next if (/^#/);
		next if (/^$/);
		chomp;
		if($_=~/\Q$keyword\E/){
			$res= $_;
			last;
		}
	}
	close $fh;
	return $res;
}

sub query_fasta_by_seqid {
	my ($db_fasta, $seqid)=@_;
	my $seq=$db->get_Seq_by_id($seqid);
    if(!defined $seq){
		warn "failed to fetch sequence for $seqid\n";
		return undef;
     }
     my $seqstr  = $seq->seq;
     $seqstr=~s/\s//g;
    chomp $seqstr;
    return $seqstr;
}

