#!/usr/bin/env perl
# ABSTRACT: Read in silva file and remove eukaryota sequences.
# AUTHOR: Arnold
# DATE: Jan 11th, 2024

use strict;
use utf8;
use warnings qw(all);
use Getopt::Long;
use Pod::Usage;


# VERSION 1.0

=head1 SYNOPSIS

    simulate_sequences [-a <alignment_fasta_path> 
                        -e <exclude_pattern>
                        -o <out_file_name>
                       ] 
                        
    alignment_fasta_path: Path to alignment file. The alignment file should be 
    single line fasta format. 
    
    exclude_pattern: The pattern to exclude.  
    
    out_file_name: The output file name / path
    
=head1 DESCRIPTION

=cut

# INPUTS: Input parameters that were used are hard coded  below
my $aln = "/nfs3/Sharpton_Lab/prod/prod_restructure/projects/arnoldhk/2025_HVR_Guide_Phylogenetic_Integration/reference_data/silva.seed_v138_2.align";
my $ex = "Eukaryota;";
my $out = "/nfs3/Sharpton_Lab/prod/prod_restructure/projects/arnoldhk/2025_HVR_Guide_Phylogenetic_Integration/reference_data/silva.seed_v138_2_Eukaryota_excluded.align";

GetOptions(
  "a=s" => \$aln,
  "e=i" => \$ex, 
  "o=s" => \$out, 
  q(help)=> \my $help,
) or pod2usage(q(-verbose) => 1);
pod2usage(q(-verbose) => 1) if $help;

open ALIGN, $aln or die("remove_eukaryotes.pl: Cannot find alignment input file at ${aln}. Please check file exists and try again.\n");
open OUT, ">$out" or die("remove_eukaryotes.pl: Cannot out file at ${out}\n");
  
my $align_count = 0;
my $header = "";
my $align_seq = "";
  
while(<ALIGN>){
  my $align_line =  $_;
  chomp $align_line;
  

  $align_count = $align_count + 1;
  if($align_count % 10000 == 0){
    print "Processing line number: " . $align_count . "\n";
  }
    
  if($align_line =~ />/){ # We are at a header
    $header = $align_line;
    $align_seq = ""; # clear the previous sequence.
      
  }else{ # We are at a sequence.
    $align_seq = $align_line;
    if($header =~ /$ex/){

    }else{
      
      print OUT $header . "\n";
      print OUT $align_seq . "\n";
    }
    $header = "";
  }
       

}
close ALIGN;
