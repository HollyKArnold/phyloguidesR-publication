#!/usr/bin/env perl
# AUTHOR: Arnold
# DATE: Jan 11th, 2024
# PURPOSE: To simulate primer sequences from a given identified hyper conserved primer region. 

use autodie;
use strict;
use utf8;
use warnings qw(all);
use Getopt::Long;
use Pod::Usage;
use Scalar::Util qw(looks_like_number); 

# VERSION 1.0

=head1 SYNOPSIS

    simulate_sequences [-a <alignment_fasta_path>> 
                        -s <start>
                        -t <stop>
                        -o <output_fasta_path>] 
                        
    alignment_fasta_path: Path to the alignment file.

    start: The start position in the alignment file from which to get bases. 
    Note that you should a priori determine this. If searching for how a 
    primer sequences matches in the alginment file, see 
    simulate_sequences_constant_length.pl as an output off that script is the
    starting base of the file.

    stop: The stop position in the alignment file from which to stop getting bases.    
    Note that you should a priori determine this. 

    out: The path to the output file. 
    
    exclude: Exclude sequences of less than this length. Set to 1 if you don't want
    to exclude any sequences. 
=head1 DESCRIPTION

=cut

# INPUTS: Input parameters are described below
 
my $align_in; 
my $start;
my $stop;
my $out;
my $exclude = 1;

GetOptions(
  "a=s" => \$align_in, 
  "s=i" => \$start,
  "t=i" => \$stop,
  "o=s" => \$out,
  "e=i" => \$exclude,
  q(help)=> \my $help,
) or pod2usage(q(-verbose) => 1);
pod2usage(q(-verbose) => 1) if $help;


# Generate sequences. 
_generate_sequences_hvr($align_in, $start, $stop, $out, $exclude);

# FUNCTION: _generate_sequences_hvr
# INPUTS: 
## $start: The start coordinate to start generating sequences from
## $stop: The position in the alignment file to stop generating sequences.
## $aln: The path to the alignment file. Single line fasta format. RNA (Us)
## $out: Path to output file and output. 
## $exclude: The length of sequence that must be generated for it to be output
## to file. Set to 1 if you want all sequences.
# OUTPUTS: 
## Simulates sequences between the start and stop position given for the 
## alignment file. 

sub _generate_sequences_hvr{
  
  # Read in arguments
  my $aln  = shift;
  my $start = shift;
  my $stop = shift;
  my $out = shift;
  my $exclude = shift;
  
  # Print progress message
  print "Generating sequences for the start and stop coordinates provided...\n";
  
  # Check ins
  if(!looks_like_number($start)){
    die("Cannot find start coordinate. Please provide an integer. \n");  
  }
  if(!looks_like_number($stop)){
    die("Cannot find start coordinate. Please provide an integer. \n");  
  }

  # Get string length to return
  my $length = $stop - $start;
  if ($length <= 0){die("Please provide valid start: ${start} and stop: ${stop} coordinates...\n")}

  open ALIGN, $aln or die("simulate_sequences: Cannot find alignment input file at ${aln}. Please check file exists and try again.\n");
  open OUT, ">$out" or die("simulate_sequences: Cannot out file at ${out}\n");
  
  my $align_count = 0;
  my $header = "";
  my $align_seq = "";
  
  while(<ALIGN>){
    
    # Get the line in the alignment file. Could be a header or a sequence. 
    my $align_line =  $_;
    chomp $align_line;

    # Set up a counter for the lines in the alignment file.
    $align_count = $align_count + 1;
    
    if($align_line =~ />/){ # We are at a header
      $header = $align_line;
      $align_seq = ""; # clear the previous sequence.
      
      
    }else{ # We are at a sequence.
      $align_seq = $align_line;
      $align_seq = substr($align_seq, $start, $length);
      $align_seq =~ s/\.//g;
      $align_seq =~ s/-//g;
      my $aln_length = length($align_seq);
      
      if($aln_length >= $exclude){
        print OUT $header . "\n";
        print OUT $align_seq . "\n";
      }
      $header = ""; # clear the previous header.
      
    }
    
  }
  close ALIGN;
  
  # Print finishing method.
  print "Generating sequences for the start and stop coordinates provided...DONE\n";   
}




