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

    simulate_sequences [-p <primer_fasta_path> 
                        -a <alignment_fasta_path>
                        -m <model_16S>
                        -l <simulated_sequence_length>] 
                        
    primer_fasta_path: Path to primer fasta file. Note that special characters 
    are allowed as long as they are standard fasta letters. See here for a list
    of allowed characters: https://en.wikipedia.org/wiki/FASTA_format. The 
    primer file should be in fasta format. Only provide one primer per file.
    
    alignment_fasta_path: The path to the alignment file. Note that the alignment
    file must be a single line fasta format. If converting SILVA files
    which are in multialign fasta format, see multiLineToSingleLineAlign.pl for
    formating help. 
    
    model_16S: The model 16S file is a full length sequence deirved from the 
    alignment fasta. Note that this file should be in single line fasta format. 
    Note that this model 16S should be derrived exactly from the alignment_fasta.
    Do not remove dashes within this file or the starting sequence will not 
    be identified correctly. A subroutine will check that this sequence matches
    to one in the alignment file prior to searching for primers. 
    
    sequence_length: An integer designating the sequence length to 
    simulate.
    
    exclude_length: Exclude sequences that are less than this value. Set to 1 
    as default.
    
    out: out file name    
=head1 DESCRIPTION

=cut

# INPUTS: Input parameters are described below
my $primer_fasta;  
my $sequence_length; 
my $align_in; 
my $model_16S;  
my $out;
my $exclude_length = 1;

GetOptions(
  "p=s" => \$primer_fasta,
  "l=i" => \$sequence_length, 
  "a=s" => \$align_in, 
  "m=s" => \$model_16S,
  "o=s" => \$out,
  "e=i" => \$exclude_length,
  q(help)=> \my $help,
) or pod2usage(q(-verbose) => 1);
pod2usage(q(-verbose) => 1) if $help;


# Check input files
my $input_check = _check_input_files($align_in, $primer_fasta, $model_16S);

# Get the primer search string
my $primer_search_string = _get_primer_string($primer_fasta);

# Get the model sequence to look for search string in.
my $model_sequence = _get_model_string($model_16S);

# Search for start site.
my $start_site = _get_primer_start($model_sequence, $primer_search_string);

# Generate sequences. 
_generate_sequences_length($align_in, $sequence_length, $start_site, $out, $exclude_length);

# FUNCTION: _generate_sequences_length
# INPUTS: 
## $length: Length of sequence to generate
## $start: The start coordinate to start generating sequences from
## $aln: The path to the alignment file. Single line fasta format. RNA (Us)
## $out: Path to output file and output. 
# OUTPUTS: 
## Searches through each alignment in the alignment file. For each sequence, 
## an output sequence is made which is the designated $length. Note that, 
## since there is variation in the number of bases within the RNA sequences
## that using a constant length may generate some sequences which extend further
## into the next hvr.
sub _generate_sequences_length{
  
  # Read in arguments
  my $aln  = shift;
  my $length = shift;
  my $start = shift;
  my $out = shift;
  my $exclude = shift;
  
  # Print progress message
  print "Generating sequences of length ${length}...\n";
  
  # Check ins
  if(!looks_like_number($start)){
    die("Cannot find start coordinate. \n");  
  }
  if(!looks_like_number($length)){
    die("Cannot find start coordinate. \n");  
  }
  open ALIGN, $aln or die("simulate_sequences: Cannot find alignment input file at ${aln}. Please check file exists and try again.\n");
  open OUT, ">$out" or die("simulate_sequences: Cannot out file at ${out}\n");
  
  my $align_count = 0;
  my $header = "";
  my $align_seq = "";
  
  while(<ALIGN>){
    my $align_line =  $_;
    chomp $align_line;
    $align_count = $align_count + 1;
    
    if($align_line =~ />/){ # We are at a header
      $header = $align_line;
      $align_seq = ""; # clear the previous sequence.
      
      
    }else{ # We are at a sequence.
      $align_seq = $align_line;
      $align_seq = substr($align_seq, $start);
      $align_seq =~ s/\.//g;
      $align_seq =~ s/-//g;
      my $aln_length = length($align_seq);
      if($aln_length <= $length){
        $align_seq = substr($align_seq, 0, $aln_length)
      }else{
        $align_seq = substr($align_seq, 0, $length)
      }
      
      
      if($aln_length >= $exclude){
        print OUT $header . "\n";
        print OUT $align_seq . "\n";
      }
      $header = ""; # clear the previous header.
      
    }
    
  }
  close ALIGN;
  
  # Print finishing method.
  print "Generating sequences of length ${length}...DONE\n";
  
}

# FUNCTION: _generate_sequences_length
# INPUTS: 
## $length: Length of sequence to generate
## $start: The start coordinate to start generating sequences from
## $alignment: The path to the alignment file. Single line fasta format. RNA (Us)
## $out_file_path: Path to output file and output. 
# OUTPUTS: 
## Searches through each alignment in the alignment file. For each sequence, 
## an output sequence is made which is the designated $start and $stop. Note that
## sequences generated may not be of the same length.
sub _generate_sequences_region{
  
}

# FUNCTION: _get_primer_start
# INPUTS: 
## $rseq: Model sequence alignment
## $regex: Primer search string 
# OUTPUTS: Returns the character at which the primer search string matches alignment
sub _get_primer_start{
  print "Getting primer start site from model sequence...\n";
  
  my $rseq  = shift;
  my $regex = shift;
  my $start_coord;
  if( $rseq =~ m/($regex)/ ){
    $start_coord = $+[0];
  }
  if(looks_like_number($start_coord)){
    print "My start site is: " . $start_coord . "\n";
    print "Getting primer start site from model sequence...DONE\n";
    return $start_coord;
  }else{
    die("Cannot find start coordinate. Consider looking at _get_primer_string to determine if you should set is_rna to different value. \n");  
  }
}

# FUNCTION: _get_model_string
# INPUTS: 
## $file: Path to model file
# OUTPUTS: Returns the model alignment string.
sub _get_model_string{
  print "Getting model string...\n";
  my $file = $_[0];
  open(FILE, $file ) || die "_get_model_string: cannot find file at ${file}. Check path try again. $!\n";
  my $seq = "";
  while(<FILE>){
    chomp $_;
    next if( $_ =~ m/^\>/ );
    $seq .= $_;
  }
  close FILE;
  return $seq;
  print "Getting model string...DONE\n";
 
}

# FUNCTION: _get_primer_string
# INPUTS: 
## $primer_path: Path to primer file
# OUTPUTS: Primer search string with translated special characters.
sub _get_primer_string{
  print "Getting primer string...\n";
  my ($primer_path) = @_;
  my $sequence = _get_primer($primer_path);
  my $search_string = _get_primer_search_string($sequence, 1);
  print "Getting primer string...DONE\n";
  return $search_string;
}

# FUNCTION _get_primer
# INPUTS:
## $file: Path to primer file
# OUTPUTS: Reads primer sequence.
sub _get_primer{ 
    my $file = $_[0];
    open( FILE, $file ) || die "_get_primer: cannot find file at ${file}. Check path try again. $!\n";
    my $seq = "";
    while(<FILE>){
        chomp $_;
        next if( $_ =~ m/^\>/ );
        $seq .= $_;
    }
    close FILE;
    return $seq;
}

# FUNCTION _get_primer_search_string
# INPUTS:
## seq: Intial primer sequence derrived from _get_primer.
## is_rna: 0 if you would like to change U to T. Note silva has U rather than T.
# OUTPUTS: Reads primer sequence.
sub _get_primer_search_string{
    my $seq = $_[0];
    my $is_rna = $_[1];
    my $string = "";
    my @chars = split //, $seq;
    for( my $i = 0; $i < length( $seq ); $i++ ){
        $string .= $chars[$i] . "(\\-|\\.)*?";
    }
    $string =~ s/R/\[A\|G\]/g;
    $string =~ s/Y/\[C\|T\]/g;
    $string =~ s/S/\[G\|C\]/g;
    $string =~ s/W/\[A\|T\]/g;
    $string =~ s/K/\[G\|T\]/g;
    $string =~ s/M/\[A\|C\]/g;
    $string =~ s/B/\[C\|G\|T\]/g;
    $string =~ s/D/\[A\|G\|T\]/g;
    $string =~ s/H/\[A\|C\|T\]/g;
    $string =~ s/V/\[A\|C\|G\]/g;
    $string =~ s/N/\[A\|C\|G\|T\]/g;
    
    if(!$is_rna){
       $string =~ s/T/U/g;
    }
    my $search_string = $string;
    print "My search string is: " . $string .  "\n";
    return $search_string;
    
}



# FUNCTION: _check_input_files
# INPUTS: 
## aln: Path to alignment file. Single line fasta format.
## primer: Path to primer file. Single line fasta format.
## model: A model sequence in single line fasta format derived from the aln 
## file for which to search for the primer sequence.
# OUTPUTS: 
## 1 if output files assessed to pass checks. 0 otherwise.
sub _check_input_files{
  
  print "Checking input files...\n";
  
  print "Checkin primer file...\n";
  my ($aln, $primer, $model) = @_;
  
  open PRIMER, $primer or die("simulate_sequences: Cannot find alignment input file at ${primer}. Please check file exists and try again.");
  
  my $primer_count = 0;
  my $primer_header = "";
  my $primer_sequence = "";
  my $primer_line = "";
  while(<PRIMER>){
    
    $primer_count = $primer_count + 1;
    $primer_line = $_;
    chomp $primer_line;
    
    if($primer_count == 1){
      $primer_header = $primer_line;
    }
    if($primer_count == 2){
      $primer_sequence = $primer_line;
    }
    
    if($primer_count > 3){
      die("FAIL: Your primer sequence has more than two lines. Consider single line fasta format.");
    }
    
  }
  print "The primer :" . $primer_header . " has sequence: " . $primer_sequence . "\n"; 
  print "Checking primer file...PASS\n";
  close PRIMER;
  
  print "Checking model alignment file...\n";
  open MODEL, $model or die("simulate_sequences: Cannot find alignment input file at ${model}. Please check file exists and try again.");
  
  # Read in the model sequence for which to search for primers.
  my $model_count = 0;
  my $model_line = "";
  my $model_seq = "";  
  my $model_header = "";
  my $model_string = "";
  
  while(<MODEL>){
    
    $model_count = $model_count + 1;
    $model_line = $_;
    chomp $model_line;
    
    if($model_count == 1){
      $model_header = $model_line;
    }
    if($model_count == 2){
      $model_string = $model_line;
    }
    
    $model_seq = $model_seq . "\n" . $model_line;
    
    if($model_count > 2){
      die("STOP: Your model sequence has more than two lines. Consider single line fasta format.");
    }
  }
  $model_seq =~ s/T/U/g;
  $model_string =~ s/T/U/g;
  
  close MODEL;
  print "Checking model alignment file...PASS\n";
  
  print "Checking alignment file...\n";
  open ALIGN, $aln or die("simulate_sequences: Cannot find alignment input file at ${primer}. Please check file exists and try again.");
  my $align_count = 0;
  my $header = "";
  my $align_seq = "";
  my $found_seq_header = 0;
  my $found_sequence = 0;
  my $found = 0;
  
  while(<ALIGN>){
    my $align_line =  $_;
    chomp $align_line;
    $align_count = $align_count + 1;
    
    # We are at a header
    if($align_line =~ />/){ 
      $header = $align_line;
      $align_seq = ""; # clear the previous sequence.
      $align_seq = $header;
      
    }else{ # We are at a sequence.
      $align_seq = $align_seq . "\n" . $align_line;
      #print "For iteration: " . $align_count . "\n\n" . $align_seq . "\n";
      $header = ""; # clear the previous header.
      if ($align_seq =~ /$model_string/){
        $found_seq_header = $align_count;
        print "Found sequence at line number: " . $align_count . "\n";
      }
      if ($align_seq =~ /$model_header/){
        $found_sequence = $align_count;
        print "Found header at line number: " . $align_count . "\n";
      }
    }
    
    if($found_sequence == $found_seq_header && $align_count != 0){
      $found = 1;
      last;
    }
    
  }
  close ALIGN;
  print "Checking alignment file...PASS\n";
  if($found == 1){
    print "Checking model found in alignment files...PASS.\n";
    return 1;
  }else{
    die("Cound not verify model sequence was derrived from alignment file\n");
  }
  
  print "Checking input files...PASS.\n";

}
