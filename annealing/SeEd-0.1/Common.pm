# #!/usr/bin/perl -w

package Common;
use strict;
use vars qw(@ISA @EXPORT @EXPORT_OK %EXPORT_TAGS $VERSION);
use Exporter();

@ISA = qw(Exporter);
@EXPORT = qw(	importGENBANK 
		importFASTA 
		num_to_str 
		str_to_num 
		mean 
		varr 
		stddev 
		chisq
		create_fasta_record 
		powers_of_4 );

# These types of records will have to be created in programs 
# using the importFASTA and importGENBANK routines

# genbank record
# my $genbank_record =
# {
#   name => "STDIN",
#   segment_count => 0,
#   segment_name => [],
#   segment_seq => [],
#   segment_len => []
# };

# FASTA record
# my $fasta_record =
# {
#   name => "STDIN",
#   contig_count => 0,
#   contig_name => [],
#   contig_seq => [],
#   contig_length => []
# };
#
# To use:
#
# open(SEQ_FILE, "<some.file")
#   or die "Could not open input_file some.file for reading: $!\n";
# &importFASTA(*SEQ_FILE, $fasta_record);
# close(SEQ_FILE);
#

my @bases = ('A', 'C', 'G', 'T');
my %base_values = ( "A" => 0, "C" => 1, "G" => 2, "T" => 3 );
my @powers_of_4;
$powers_of_4[0] = 1;
for (my $i = 1; $i < 16; $i++)
{
  push (@powers_of_4, $powers_of_4[$i-1]<<2);
}

# FASTA record to be used implicitly by calling program
my $fr =
{
  name => "",
  contig_count => 0,
  contig_name => [],
  contig_seq => [],
  contig_length => []
};

my %fs;
$fs{$fr->{name}} = $fr;
my $fs_counter = 0;

#%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
# The subrouting create_fasta_record does what its name suggests, basically
# imports a fasta file and saves it in a fasta_record structure, which it returns
# relevant information in a structure
# Arguments:
# 1. Filehandle
# Returns:
# 1. fasta_record (fr) reference
#%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
sub create_fasta_record
{
  my ($file_ptr) = @_;
  $fs{$fs_counter}->{name} = $fs_counter;
  &importFASTA($file_ptr, $fs{$fs_counter});
  return $fs{$fs_counter++};
}

#%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
# The subrouting importGENBANK imports a GENBANK file, storing the contigs and
# relevant information in a structure
# Arguments:
# 1. Filehandle
# 2. Record pointer to store the data
#%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
sub importGENBANK
{
  my $filehandle = $_[0];
  my $genbank_record = $_[1];
  my $segment_count = -1;
  my $seq = "";
  my $definition = "";
  my $len = 0;
  my $mode = "header";
  while(<$filehandle>)
  {
    if (/^\/\//)
    {
      if ($mode eq "sequence")
      {
        $seq =~ y/[U]/[T]/;
        push @{$genbank_record->{segment_seq}}, $seq;
        $seq = "";
        $definition = "";
      }
      $mode = "header";
    }
    elsif ($mode eq "sequence")
    {
      s/[^ACTGU]//gi;
      $seq .= $_;
    }
    elsif (/^LOCUS\s*\S*\s*(\d*)\s.*$/)
    {
      $segment_count++;
      push @{$genbank_record->{segment_len}}, $1 + 0;
    }
    elsif (/^ORIGIN/)
    {
      $mode = "sequence";
    }
    elsif (/^DEFINITION  (.*)$/)
    {
      $mode = "definition";
      $definition .= $1;
    }
    elsif ($mode eq "definition")
    {
      if (/^ {12}(.*)/)
      {
        $definition .= " $1";
      }
      else
      {
        $genbank_record->{segment_name}->[$segment_count] = $definition;
        $mode = "header";
      }
    }
  }
  $genbank_record->{segment_count} = ++$segment_count;
}
#%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

#%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
# The subrouting importFASTA imports a FASTA file, storing the contigs and
# relevant information in a structure
# Arguments:
# 1. Filehandle
#
# A second argument of the title of the record maybe added if multiple file
# input will be allowed.
#%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
sub importFASTA
{
  my $filehandle = $_[0];
  my $fasta_record = $_[1];
  my $contig_count = -1;
  my $seq = "";
  my $len = 0;
  while(<$filehandle>)
  {
    if (!/^>/)
    {
      s/\s//gi;
      $seq .= $_;
    }
    else
    {
      /^>(.*)$/;
      push @{$fasta_record->{contig_name}}, $1;
      if ($contig_count >= 0)
      {
        push @{$fasta_record->{contig_len}}, length($seq);
	$seq =~ tr/a-z/A-Z/;
        push @{$fasta_record->{contig_seq}}, $seq;
      }
      $seq = "";
      $contig_count++;
    }
  }
  push @{$fasta_record->{contig_len}}, length($seq);
  $seq =~ tr/a-z/A-Z/;
  push @{$fasta_record->{contig_seq}}, $seq;
  $fasta_record->{contig_count} = ++$contig_count;
}
#%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%


#%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
# The subroutine 'mean' will calculate the mean of the elements of a list.
# Arguments:
# 1. Pointer to list
#%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
sub mean
{
  my ($arr) = @_;
  my $num_of_elements = @$arr;
  my $sum = 0;
  for (my $i = 0; $i < $num_of_elements; $i++)
  {
    $sum += $arr->[$i];
  }
  return ($sum/$num_of_elements);
}
#%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%


#%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
# The subroutine 'varr' will calculate the variance of the elements of a list.
# Arguments:
# 1. Pointer to list
#%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
sub varr
{
  my ($arr) = @_;
  my $num_of_elements = @$arr;
  my $sum = 0;
  my $mean = &mean($arr);
  for (my $i = 0; $i < $num_of_elements; $i++)
  {
    $sum += ($arr->[$i] - $mean)**2;
  }
  if ($num_of_elements > 1)
  {
    return ($sum/($num_of_elements-1));
  }
  else
  {
    return 0;
  }
}
#%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

#%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
## Given an array of numeric values, return the standard deviation of the values.
## Arguments:
## 1. Pointer to an array of numeric values (pointer)
##%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
sub stddev
{
  return (sqrt(&varr($_[0])));
}

#%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
# Given two values, namely the observed and expected values of some quantity
# calculates the chi square statistic.
# Arguments:
# 1. Observed value
# 2. Expected value
##%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
sub chisq
{
  return (($_[0] - $_[1])**2)/$_[1];
}


#%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
# The next subroutine takes as input a number and a size and returns a
# corresponding string of that size composed by DNA bases.
# Arguments:
#   1. Number to be converted to DNA string
#   2. Desired size of string
# Note: 1st_argument <= 2nd_argument**4
# Definitions:
#%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
sub num_to_str
{
  my $num;

  $num = $_[0];
  my $return_value = '';
  for (my $j = 0; $j < $_[1]; $j++)
  {
    $return_value = $bases[$num % 4].$return_value;
    $num = $num>>2;
  }
  $return_value;
}
#%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%


#%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
# The str_to_num subroutine takes as input a string returns a corresponding
# number, such that there is a one-to-one mapping from the strings based on
# the {A, C, G, T} alphabet and the natural numbers.
# Arguments:
#   1. String to be converted/indexed to a number
# Definitions:
#%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
sub str_to_num
{
  my $str = $_[0];
  my $return_value = 0;
  for (my $j = 0; $j < length($str); $j++)
  {
    my $char = substr($str, $j, 1);
    if (exists($base_values{$char}))
    {
      $return_value += $base_values{substr($str, $j, 1)}*$powers_of_4[length($str) - $j - 1];
    }
    else
    {
      return -1;
    }
  }
  return $return_value;
}
#%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

1;
