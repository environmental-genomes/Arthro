use strict;
use warnings;
use Getopt::Long;

my %opts;
  GetOptions ("Category!" => \$opts{Category},    
              "Subcategory!"   => \$opts{Subcategory},      
              "Subsystem!"  => \$opts{Subsystem},  
			  "Role!" =>\$opts{Role},
			  "help!"=> \$opts{help} )
  or die("Error in command line arguments\n");
  
if (defined $opts{Category} || defined $opts{Subcategory} || defined $opts{Role} ||defined $opts{Subsystem}){
	if(scalar(@ARGV)>0){
		my $count={};
		my @files=();
		foreach my $file (@ARGV){
			
			open my $fh,'<', $file;
			if(!defined $fh){
				warn "$file: $@\n";
				next;
			}
			my $filter='';
			while(<$fh>){
				chomp;
				next if(/^$/);
				next if(/^\s+$/);
				next if (/^Category\tSubcategory\tSubsystem/);
				my ($Category, $Subcategory, $Subsystem,$Role,$Features)=split /\t/;
				$Category=~s/\s+$//; $Category=~s/^\s+//;
				$Subcategory=~s/\s+$//; $Subcategory=~s/^\s+//;
				$Subsystem=~s/\s+$//; $Subsystem=~s/^\s+//;
				$Role=~s/\s+$//;  $Role=~s/^\s+//;
				$Features=~s/\s+$//;  $Features=~s/^\s+//;
				
				# next if (defined $opts{Category} && ! grep( /$Category/, @{$opts{Category}} ));
				# next if (defined $opts{Subcategory} && ! grep( /$Subcategory/, @{$opts{Subcategory}} ));
				# next if (defined $opts{Subsystem} && ! grep( /$Subsystem/, @{$opts{Subsystem}} ));
				# next if (defined $opts{Role} && ! grep( /$Role/, @{$opts{Role}} ));
				$filter=$filter."|".$Category if (defined $opts{Category});
				$filter=$filter."|".$Subcategory if (defined $opts{Subcategory});
				$filter=$filter."|".$Subsystem if (defined $opts{Subsystem});
				$filter=$filter."|".$Role if (defined $opts{Role});
				my @genes=split /,/,$Features;
				$count->{$filter}->{$file}+=scalar(@genes);
				$filter='';
			}
			push @files, $file;
			# print "$file\t$total";
			# if(defined $opts{Category}){
				# print "\t\"";
				# if (scalar(@{$opts{Category}})>1) {print join '|', @{$opts{Category}};}
				# else {print $opts{Category}->[0];}
				# print "\"";
			# }
			# if(defined $opts{Subcategory}){
				# print "\t\"";
				# if (scalar(@{$opts{Subcategory}})>1) {print join '|', @{$opts{Subcategory}};}
				# else {print $opts{Subcategory}->[0];}
				# print "\"";
			# }
			# if(defined $opts{Subsystem}){
				# print "\t\"";
				# if (scalar(@{$opts{Subsystem}})>1) {print join '|', @{$opts{Subsystem}};}
				# else {print $opts{Subsystem}->[0];}
				# print "\"";
			# }
			# if(defined $opts{Role}){
				# print "\t\"";
				# if (scalar(@{$opts{Role}})>1) {print join '|', @{$opts{Role}};}
				# else {print $opts{Role}->[0];}
				# print "\"";
			# }			
			# print "\n";
		}
		print "group";
		print "\t\"$_\"" foreach @files;
		print "\n";
		foreach my $FILTER(keys %{$count}){
			my $FILTER2=$FILTER;
			$FILTER2=~s/^\|//;
			print "\"$FILTER2\"";
			foreach my $FILE (@files){
				print "\t";
				print defined ($count->{$FILTER}->{$FILE} ) ? $count->{$FILTER}->{$FILE} :0;
				
			}
			print "\n";
		}
	}
	else {
		die "  at least one file to be specified\n";
	}

}
else {
	warn "\n  At least one of Category, Subcategory or Role must be provided for selecting genes\n\n";
	&help;
}

sub help{
	warn " usage: perl $0 [optiions]\n";
	warn "   options:\n";
	warn "      --Category      count genes number by Category\n";
	warn "      --Subcategory   count genes number by Subcategory\n";
	warn "      --Role          count genes number by Role\n";
	warn "      --Subsystem     count genes number by SUbsystem\n";
	warn "      --help           show this message\n";
	exit(0);
}
