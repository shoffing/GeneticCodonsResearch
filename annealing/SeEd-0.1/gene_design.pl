#!/usr/bin/perl -w

# Codon pair bias sequence design program
#
# Version 0.1
#
# Creation: 2005-07-18
# 
# This utility will be used to design coding sequences using
# specified codon distributions, maximizing or minimizing
# their codon pair score and eliminating specified restriction
# sites.
# The simulated annealing technique is employed for the 
# optimization of the codon pair score.
#
#
#  Copyright (c) 2008 Research Foundation of the State University of
#  New York. All rights reserved.
#
#  Redistribution and use of the source code, with or without
#  modification, are permitted provided that the following conditions are met:
#
#  1.      Redistributions must reproduce the above copyright notice, this
#  list of conditions and the following disclaimer in the  documentation
#  and/or other materials provided with the distribution.  Redistributions of
#  source code must also reproduce this information in the source code itself.
#
#  2.      If the program is modified, redistributions must include a notice
#  (in the same places as above) indicating that the redistributed program is
#  not identical to the version distributed by the original creator.
#
#  3.      The name of the original creator may not be used to endorse or
#  promote products derived from this software without specific prior written
#  permission.
#
#  We also request that use of this software be cited in publications as
#
#     J. R. Coleman, D. Papamichail, S. Skiena, B. Futcher, E. Wimmer
#     and S. Mueller, "Virus attenuation by genome-scale changes in
#     codon-pair bias: a general method for developing viral vaccines"
#     Science, ...
#     Code available at http://www.cs.miami.edu/~dimitris/SeEd-perl/.
#
#  THIS SOFTWARE IS PROVIDED BY THE ORIGINAL AUTHOR ``AS IS'' AND  ANY
#  EXPRESS OR IMPLIED WARRANTIES, INCLUDING, BUT NOT LIMITED TO, THE  IMPLIED
#  WARRANTIES OF MERCHANTABILITY AND FITNESS FOR A PARTICULAR PURPOSE  ARE
#  DISCLAIMED. IN NO EVENT SHALL THE ORIGINAL AUTHOR BE LIABLE  FOR ANY
#  DIRECT, INDIRECT, INCIDENTAL, SPECIAL, EXEMPLARY, OR CONSEQUENTIAL  DAMAGES
#  (INCLUDING, BUT NOT LIMITED TO, PROCUREMENT OF SUBSTITUTE GOODS  OR
#  SERVICES; LOSS OF USE, DATA, OR PROFITS; OR BUSINESS INTERRUPTION)  HOWEVER
#  CAUSED AND ON ANY THEORY OF LIABILITY, WHETHER IN CONTRACT, STRICT
#  LIABILITY, OR TORT (INCLUDING NEGLIGENCE OR OTHERWISE) ARISING IN ANY WAY
#  OUT OF THE USE OF THIS SOFTWARE, EVEN IF ADVISED OF THE POSSIBILITY OF
#  SUCH DAMAGE.
#

use Common;
use Getopt::Std;

$MAX_SEQUENCE_LENGTH = 1000000;

&getopts("o:l:c:r:d:p:i:hs", \%args); # -v, -D, -o ARG, sets $args{v}, $args{D}, $args{o}

if ($#ARGV == -1)
{
  print "\nUsage: gene_design.pl [-l lock_file] [-c coding_regions_file] [-r restriction_site_file]\n";
  print "                      [-d distribution_file] [-h] [-i iterations] [-s] [-o output_file] <input_file> \n\n";
  print "  <input_file>: input file in fasta format (Put '-' for STDIN)\n";
  print "  -l lock_file: import locks from lock file (Default: no locks)\n";
  print "  -c coding_regions_file: Import coding regions from file (Default: whole input)\n";
  print "  -r restriction_site_file: Specify restriction sites to eliminate (Default: none)\n";
  print "  -d distribution_file: Import codon distribution for codon shuffling (Default: input data codon distribution)\n";
  print "  -p codon pair preference file: Import codon pair preferences\n";
  print "  -i iterations: Codon pair min/max algorithm iterations\n";
  print "  -s Apply simulated annealing method for min/max calculation\n";
  print "  -h Shuffle codons achieving maximum hamming distance\n";
#  print "  -i input_file: use input_file as fasta input (Put '-' for STDIN)\n";
  print "  -o output_file: use output_file as fasta output (Put '-' for STDOUT)\n";
  print "\n";
  exit;
}

$args{f} = $ARGV[0];
#$args{i} = "-" unless defined $args{i};
$args{o} = "-" unless defined $args{o};
$args{i} = 10000 unless defined $args{i};

# Define DEBUG variables

$DEBUG_CODON_IMPORT = 0;
$DEBUG_CODON_PAIR_IMPORT = 0;
$DEBUG_WORKSPACE = 0;
$DEBUG_WORKSPACE_LOCK = 0;
$DEBUG_WORKSPACE_TOTAL = 0;
$DEBUG_IMPORT_DIST = 0;
$DEBUG_SOURCE_DIST = 0;
$DEBUG_PERMUTATION = 0;
$DEBUG_RESTRICTION_SITES = 0;
$DEBUG_PATTERN_MATCH = 0;
$DEBUG_MODULO = 0;
$DEBUG_PATTERN = 0;
$DEBUG_ELIMINATE_RS = 0;

$CREATE_CODON_PAIR_STATISTICS = 1;
$CREATE_CODON_PAIR_SCORE_HISTOGRAM = 0;
$CREATE_RANDOM_CODON_PAIR_BIASES = 0;
$EVALUATE_AVERAGE_SCORES_FOR_AA_PAIRS = 0;
$RANDOM_CODON_PAIR_BIAS_CREATION_REPEATS = 10000;

# Simulated Annealing variables

my $sa_t = 1;		# Temperature
my $sa_k = 0.4;		# The k (c) constant
my $sa_iter = 100;	# Number of iteration at a specific temperature
my $sa_i = 0;		# Current iteration index
my $sa_a = 0.999;	# Temperature decrement factor
my $sa_acc_tran = 0;	# Indicated if a transition was accepted during the last iteration
my $sa_exit = 0;	# Indicated if the method should be exited, since no transition
			# was made in the last iteration.

# FASTA record
$fasta_record =
{
  contig_count => 0,
  contig_name => [],
  contig_seq => [],
  contig_len => []
};

# AA record
$aa_record =
{
  aa => "na",		# Amino acid symbol
  num_of_codons => 0,	# Number of codons for this AA
  num_per_type => [],	# how many of each type of codon we have in our sequence
  location_list => [],	# locations of the codons for this AA in our sequence
  location_type => [],	# type of codon for each location
  at_least_two_different => 0,	# greater than 0 if the AA has at least two different codons in our sequence
  codon_map => [],
  perm => [],
  rev_perm => [],
  new_type => []
};

$aa_byname{$aa_record->{aa}} = $aa_record; # Create a hash of aa records

open (LOG_FILE, ">output/design.log");
my $file_name = $args{f};
if ($file_name =~ /\/([^\/\.]+)\.fasta$/)
{
  $file_name = $1;
#  print $file_name, "\n";
}
open (STATISTICS_FILE, ">output/codon_pair_bias_statistics.".$file_name.".".$args{i})
  or die "Could not open file >output/codon_pair_bias_statistics.".$file_name.".$RANDOM_CODON_PAIR_BIAS_CREATION_REPEATS for writing!\n";

open(SEQ_FILE, $args{f})
  or die "Could not open input_file ", $args{f}, " for reading: $!\n";
&importFASTA(*SEQ_FILE, $fasta_record);
close(SEQ_FILE);

if ($fasta_record->{contig_count} > 1)
{
  die "\nCurrent version does not support more than one fasta input sequence\n\n";
}

if ($fasta_record->{contig_len}->[0] > $MAX_SEQUENCE_LENGTH)
{
  die "\nMaximum sequence length supported is $MAX_SEQUENCE_LENGTH. Limit exceeded.\n\n";
}

$num_of_locks = 0;
for (my $i = 0; $i < $fasta_record->{contig_len}->[0]; $i++)
{
  $lock_mask[$i] = 1;	# lock mask is the table that shows which locations are unlocked
}
%lock_start = ();
if ($args{l})
{
  &import_locks;
}

$num_of_coding_regions = 0;
if ($args{c})
{
  &import_coding_regions;
}
else
{
  $num_of_coding_regions = 1;
  $coding_regions[0][0] = 0;
  $coding_regions[0][1] = $fasta_record->{contig_len}->[0];
}

&import_codon_info;
&create_work_space;
&calculate_source_distribution;

if ($args{d})
{
  &import_distribution;
}
else
{
  $target_dist = $source_dist;
}

$num_of_restriction_sites = 0;
if ($args{r})
{
  &import_restriction_sites;
}

if ($args{p})
{
  &import_codon_pair_info;
}

&eliminate_restriction_sites;
=for other purposes
&shuffle_codons;
&create_new_distribution_sequence;
&output_new_sequence;
&output_log_statistics;
&output_log_distributions;
&output_log_protein_equivalence;
&output_log_hamming_distance;
=cut
close (LOG_FILE);
close (STATISTICS_FILE);

# Codon pair structure
$codon_pair_record =
{
  bases => "na",	# The six bases that consist the codon pair
  AAs => "na",		# Two Amino acid symbols
  num => 0,	 	# Number of codon pairs appearing
  exp_num => 0,		# Number of codon pairs expected
  value => 0,		# value for this codon pair (log of ration observed/predicted with discounting)
  chisq => 0		# chi square value for observed and expected
};

$codon_pair{$codon_pair_record->{bases}} = $codon_pair_record; # Create a hash of codon pair records

#%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
# Function that imports the codon pair information from the codon pair file and 
# stores it in a hash.
#%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
sub import_codon_pair_info
{
# print STATISTICS_FILE "---------------------------------------------------------\n";
  open (CODON_PAIR_FILE, $args{p})
    or die "Could not open codon pair info file ", $args{p}, " for reading: $!\n";
  my $garbage = <CODON_PAIR_FILE>;
  while(<CODON_PAIR_FILE>)
  {
    chop($_);
    /^(\S+)\s+(\S+)\s+(\S+)\s+(\S+)/;
#    print "$1 $2 $3 $4\n";
    if (!$aa_pair_freq{$1})
    {
      $aa_pair_freq{$1} = $4;
    }
    else
    {
      $aa_pair_freq{$1} += $4;
    }
    my $first_codon = substr($2, 0, 3);
    my $second_codon = substr($2, 3, 3);
    if (!$codon_freq{$first_codon})
    {
      $codon_freq{$first_codon} = $4/2;
    }
    else
    {
      $codon_freq{$first_codon} += $4/2;
    }
    if (!$codon_freq{$second_codon})
    {
      $codon_freq{$second_codon} = $4/2;
    }
    else
    {
      $codon_freq{$second_codon} += $4/2;
    }
    my $first_aa = substr($1, 0, 1);
    my $second_aa = substr($1, 1, 1);
    $codon_to_aa_code{$first_codon} = $first_aa;
    $codon_to_aa_code{$second_codon} = $second_aa;
    if (!$aa_freq{$first_aa})
    {
      $aa_freq{$first_aa} = $4/2;
    }
    else
    {
      $aa_freq{$first_aa} += $4/2;
    }
    if (!$aa_freq{$second_aa})
    {
      $aa_freq{$second_aa} = $4/2;
    }
    else
    {
      $aa_freq{$second_aa} += $4/2;
    }
    $codon_pair{$2}->{AAs} = $1;
    my $discounted_num = &discount($4 + 0);
    $codon_pair{$2}->{num} = $discounted_num;
    my $exp_num = $3;
  }
  # Calculate the expected number of codon pairs and the codon pair values
  @codon_pair_scores = ();
  $min_codon_pair_score = 1000000;
  $max_codon_pair_score = -1000000;
  @codon_pair_chisq_values = ();
  $min_codon_pair_chisq_value = 1000000;
  $max_codon_pair_chisq_value = -1000000;
  foreach my $key (keys %codon_pair)
  {
    my $first_codon = substr($key, 0, 3);
    my $second_codon = substr($key, 3, 3);
    my $first_aa = substr($codon_pair{$key}->{AAs}, 0, 1);
    my $second_aa = substr($codon_pair{$key}->{AAs}, 1, 1);
    my $aa_pair = $codon_pair{$key}->{AAs};
    $codon_pair{$key}->{exp_num} = ($codon_freq{$first_codon}/$aa_freq{$first_aa}) *
      ($codon_freq{$second_codon}/$aa_freq{$second_aa}) * $aa_pair_freq{$codon_pair{$key}->{AAs}};
    my $value = log($codon_pair{$key}->{num}/$codon_pair{$key}->{exp_num});
    $codon_pair{$key}->{value} = $value;
    my $value_chisq = &chisq($codon_pair{$key}->{num}, $codon_pair{$key}->{exp_num});
    $codon_pair{$key}->{chisq} = $value_chisq;
    push (@codon_pair_scores, $value);
    push (@codon_pair_chisq_values, $value_chisq);
    if ($value > $max_codon_pair_score) {$max_codon_pair_score = $value;}
    if ($value < $min_codon_pair_score) {$min_codon_pair_score = $value;}
    if ($value_chisq > $max_codon_pair_chisq_value) {$max_codon_pair_chisq_value = $value_chisq;}
    if ($value_chisq < $min_codon_pair_chisq_value) {$min_codon_pair_chisq_value = $value_chisq;}
    if ($aa_pair_value_sum{$aa_pair})
    {
      $aa_pair_value_sum{$aa_pair} += $codon_pair{$key}->{value};
    }
    else
    {
      $aa_pair_value_sum{$aa_pair} = $codon_pair{$key}->{value};
    }
  }
  if ($EVALUATE_AVERAGE_SCORES_FOR_AA_PAIRS)
  {
    # Evaluate the average scores for each AA pair
    foreach my $key (keys %aa_pair_value_sum)
    {
      my $mult1 = $aa_codon_num{$code_to_aa{substr($key, 0, 1)}};
      my $mult2 = $aa_codon_num{$code_to_aa{substr($key, 1, 1)}};
      $aa_pair_average_score{$key} = $aa_pair_value_sum{$key}/($mult1*$mult2);
      print STATISTICS_FILE "$key\t$mult1\t$mult2\t", $aa_pair_average_score{$key}, "\n";
    }
  }
  if ($CREATE_CODON_PAIR_SCORE_HISTOGRAM)
  {
    open (HIST_FILE, ">output/".$file_name."_score_histogram")
      or die "Cannot open file output/".$file_name."_score_histogram for writing!\n";
    @sorted_codon_pair_scores = sort {$a <=> $b} @codon_pair_scores;
    my $step = 0.02;
    my $current_level = int(($min_codon_pair_score - $step)*(1/$step))*$step;
    my $counter = 0;
    foreach my $current_value (@sorted_codon_pair_scores)
    {
      if ($current_value > $current_level)
      {
        while ($current_value > $current_level)
        {
          print HIST_FILE "$current_level\t$counter\n";
          $counter = 0;
          $current_level += $step;
        }
        $counter = 1;
      }
      else
      {
        $counter++;
      }
    }
    close(HIST_FILE);
  }

  if ($CREATE_CODON_PAIR_STATISTICS)
  {
    my $total_codons = 0;
    my $total_exp_num = 0;
    my $total_value_sum = 0;
    foreach $cp (keys %codon_pair)
    {
      if ($DEBUG_CODON_PAIR_IMPORT)
      {
        print STATISTICS_FILE "$cp\t"; 
        print STATISTICS_FILE $codon_pair{$cp}->{num}, "\t"; 
        print STATISTICS_FILE $codon_pair{$cp}->{exp_num}, "\t";
        print STATISTICS_FILE $codon_pair{$cp}->{value}, "\n";
      }
      $total_codons += $codon_pair{$cp}->{num};
      $total_exp_num += $codon_pair{$cp}->{exp_num};
      $total_value_sum += $codon_pair{$cp}->{value};
      $total_codons += 1;
    }
    $min_score = &find_cp_opt(0, 1);
    $max_score = &find_cp_opt(1, 1);

    print STATISTICS_FILE "\n---------------- STATISTICS FOR HOMO SAPIENS CODON PAIR SCORES ----------------\n\n";
    print STATISTICS_FILE "log(observed/expected) scores:\n";
    print STATISTICS_FILE "  Mean = ", &mean(\@codon_pair_scores), "\n";
    print STATISTICS_FILE "  Standard Deviation = ", &stddev(\@codon_pair_scores), "\n";
    print STATISTICS_FILE "  Min value = ", $min_codon_pair_score, "\n";
    print STATISTICS_FILE "  Max value = ", $max_codon_pair_score, "\n";
    if ($args{s})
    {
      print STATISTICS_FILE "\n---------------- SIMULATED ANNEALING METHOD APPROXIMATIONS ----------------\n\n";
      print STATISTICS_FILE "SA parameters: iter = $sa_iter, k = $sa_k, a = $sa_a\n\n";
    }
    else
    {
      print STATISTICS_FILE "\n---------------- GRADIENT DESCENT METHOD APPROXIMATIONS ----------------\n\n";
    }
    print STATISTICS_FILE "Minimum score = $min_score\n";
    print STATISTICS_FILE "Maximum score = $max_score\n";
  }
  if ($CREATE_RANDOM_CODON_PAIR_BIASES)
  {
    my $min_score = 1000000000;
    my $max_score = -1000000000;
    my @score_table = ();
    my @hamming_dist_table = ();
    my @codon_diff_table = ();
    my $num_of_seq_below_initial = 0;
    my $num_of_seq_above_initial = 0;
    my $total_num_of_codons = 0;
    for (my $r = 0; $r < $RANDOM_CODON_PAIR_BIAS_CREATION_REPEATS; $r++)
    {
      $new_seq = $fasta_record->{contig_seq}->[0];
      foreach my $aa (@aas)
      {
        if ($aa_codon_num{$aa} > 1)
        {
          my $num_of_codons = $aa_byname{$aa}->{num_of_codons};
          my ($permuted_index, $rp) = &permutation($num_of_codons);
          for (my $i = 0; $i < $num_of_codons; $i++)
          {
            my $codon_at_i = $type_to_codon{$aa}->[$aa_byname{$aa}->{location_type}->[$i]];
  	    my $new_location_for_codon_i = $aa_byname{$aa}->{location_list}->[$permuted_index->[$i]];
  	    substr($new_seq, $new_location_for_codon_i, 3) = $codon_at_i;
          }
        }
      }
      my $new_score = 0;
      for (my $i= 0; $i < $num_of_coding_regions; $i++)
      {
        for (my $j = $coding_regions[$i][0]; $j < $coding_regions[$i][1] - 3; $j += 3)
        {
          my $cp = substr($new_seq, $j, 6);
          if (exists $codon_pair{$cp})
	  {
            $new_score += $codon_pair{$cp}->{value};
	  }
        }
      }
      if ($new_score > $max_score) {$max_score = $new_score;}
      if ($new_score < $min_score) {$min_score = $new_score;}
      if ($new_score < $cp_initial_score) 
      {
        $num_of_seq_below_initial++;
      }
      else
      {
        $num_of_seq_above_initial++;
      }
      push (@score_table, $new_score);
      # Check if new sequence encodes the same protein as the old one
      # what is the hamming distance and the codon distance and what distributions
      # they are encoding (if they are the same)
      my $aa_that_differ = 0;
      my $codons_that_differ = 0;
      my $hamming_distance = 0;
      my %codon_dist1 = ();
      my %codon_dist2 = ();
      for (my $i= 0; $i < $num_of_coding_regions; $i++)
      {
	for (my $j = $coding_regions[$i][0]; $j < $coding_regions[$i][1]; $j += 3)
	{
	  my $in_codon = substr($fasta_record->{contig_seq}->[0], $j, 3);
	  my $out_codon = substr($new_seq, $j, 3);
	  if (exists $codon_dist1{$in_codon}) { $codon_dist1{$in_codon}++; }
	    else { $codon_dist1{$in_codon} = 1; }
	  if (exists $codon_dist2{$out_codon}) { $codon_dist2{$out_codon}++; }
	    else { $codon_dist2{$out_codon} = 1; }
	  if ($codon_to_aa{$in_codon} ne $codon_to_aa{$out_codon})
	  {
	    $aa_that_differ++;
	  }
	  elsif ($in_codon ne $out_codon)
	  {
	    $codons_that_differ++;
	  }
	  $total_num_of_codons++;
	}
	for (my $j = $coding_regions[$i][0]; $j < $coding_regions[$i][1]; $j += 1)
	{
	  my $in_base = substr($fasta_record->{contig_seq}->[0], $j, 1);
	  my $out_base = substr($new_seq, $j, 1);
	  if ($in_base ne $out_base)
	  {
	    $hamming_distance++;
	  }
	}
      }
      push (@hamming_dist_table, $hamming_distance);
      push (@codon_diff_table, $codons_that_differ);
    }
    $total_num_of_codons /= $RANDOM_CODON_PAIR_BIAS_CREATION_REPEATS;
    # Now create a histogram file
    open (HIST_FILE, ">output/".$file_name."_random_score_histogram.$RANDOM_CODON_PAIR_BIAS_CREATION_REPEATS")
      or die "Cound not open file output/".$file_name."_random_score_histogram.$RANDOM_CODON_PAIR_BIAS_CREATION_REPEATS for writing!\n";
    my @sorted_scores = sort {$a <=> $b} @score_table;
    my $step = 0.1;
    my $current_level = int(($min_score - $step)*(1/$step))*$step;
    my $counter = 0;
    foreach my $current_value (@sorted_scores)
    {
      if ($current_value > $current_level)
      {
        while ($current_value > $current_level)
        {
          print HIST_FILE "$current_level\t$counter\n";
          $counter = 0;
          $current_level += $step;
        }
        $counter = 1;
      }
      else
      {
        $counter++;
      }
    }
    print HIST_FILE "$current_level\t$counter\n";
    close(HIST_FILE);
    print STATISTICS_FILE "\n------------ RANDOM POLIO SCORE GENERATION for $RANDOM_CODON_PAIR_BIAS_CREATION_REPEATS repeats -------------\n";
    print STATISTICS_FILE "         (KEEPING CODON BIAS CONSTANT AND CHANGING CODON PAIR BIAS)\n\n";
    my $sc_mean = &mean(\@score_table);
    my $sc_stddev = &stddev(\@score_table);
    print STATISTICS_FILE "Minimum Score = $min_score\n";
    print STATISTICS_FILE "Maximum Score = $max_score\n";
    print STATISTICS_FILE "Mean = ", $sc_mean, "\n";
    print STATISTICS_FILE "Standard Deviation = ", $sc_stddev, "\n";
    print STATISTICS_FILE "\nStandard deviations from the mean the sequence score appears: ";
    print STATISTICS_FILE abs($cp_initial_score - $sc_mean)/$sc_stddev;
    print STATISTICS_FILE "\nNumber of random sequences with score above (or equal) original: $num_of_seq_above_initial/$RANDOM_CODON_PAIR_BIAS_CREATION_REPEATS";
    print STATISTICS_FILE "\nNumber of random sequences with score below original: $num_of_seq_below_initial/$RANDOM_CODON_PAIR_BIAS_CREATION_REPEATS";
    print STATISTICS_FILE "\n\nMean and standard deviation of random sample hamming distance from original sequence: ";
    print STATISTICS_FILE &mean(\@hamming_dist_table), ", ", &stddev(\@hamming_dist_table), "\n";
    print STATISTICS_FILE "Mean and standard deviation of random sample codon difference from original sequence: ";
    print STATISTICS_FILE &mean(\@codon_diff_table), ", ", &stddev(\@codon_diff_table), "\n";
    print STATISTICS_FILE "\nTotal number of codons: $total_num_of_codons\n";
  }
}

sub discount
{
  return $_[0];
}

#%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
# Procedure to calculate the current, minimum and maximum codon pair scores. It
# will do so with gradient descent, exchanging codons, whenever there is a score
# improvement, in random.
#%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
sub find_cp_opt
{
  my $max_condition = $_[0];
  my $sa_on = $_[1];
  $opt_seq = $fasta_record->{contig_seq}->[0];
  # Calculate the score of the sequence to start with
  $cp_initial_score = 0;
  for (my $i= 0; $i < $num_of_coding_regions; $i++)
  {
    for (my $j = $coding_regions[$i][0]; $j < $coding_regions[$i][1] - 3; $j += 3)
    {
      my $cp = substr($fasta_record->{contig_seq}->[0], $j, 6);
      if (exists $codon_pair{$cp})
      {
        $cp_initial_score += $codon_pair{$cp}->{value};
      }
    }
  }
  unless ($max_condition)
  {
    print STATISTICS_FILE "\nInitial codon pair score for our sequence = $cp_initial_score\n\n";
  }
  $opt_score = $cp_initial_score;
  # Create two hashes of the border locations
  for (my $i= 0; $i < $num_of_coding_regions; $i++)
  {
    $start_border{$coding_regions[$i][0]} = 1;
    $end_border{$coding_regions[$i][1]-2} = 1;
  }
  # Now do a randomized gradient descent
  if ($args{i})
  {
    $repeats = $args{i} + 0;
  }
  else
  {
    $repeats = 100000;
  }
  my $cond = "min";
  my $method = "grad_desc";
  if ($max_condition)
  {
    $cond = "max";
  }
  if ($args{s})
  {
    $method = "sim_anneal";
  }
  open(OPT_FILE, ">output/$cond\_score_$method\_".$file_name.".$repeats")
    or die "Could not open output file output/min_score_gradient_descent_".$file_name.".$repeats for writing: $!\n";
  for (my $i = 0; $i < $repeats; $i++)
  {
    my $aa = "";
    do
    {
      my $random_seq_pos = int(rand $all_region_seq_length);
      my $random_pos_region = &coding_region_of_position($random_seq_pos);
      my $random_pos_codon_start = $coding_regions[$random_pos_region][0] + $random_seq_pos - $region_start[$random_pos_region] - (($random_seq_pos - $region_start[$random_pos_region]) % 3);
      my $random_pos_codon = substr($fasta_record->{contig_seq}->[0], $random_pos_codon_start, 3);
      $aa = $codon_to_aa{$random_pos_codon};
    } until (($aa_codon_num{$aa} > 1) && ($aa ne "Ter") && ($aa_byname{$aa}->{at_least_two_different}));
    my $num_of_codon_types_of_random_aa = $aa_codon_num{$aa};
    if ($num_of_codon_types_of_random_aa < 2)
    {
      die "Number of codon types for AA $aa less than 2. Critical error!\n";
    }
    my $num_of_codons_for_this_aa = $aa_byname{$aa}->{num_of_codons};
    my $rand_pos1 = int(rand $num_of_codons_for_this_aa);
    my $rand_pos2 = int(rand $num_of_codons_for_this_aa);
    while (substr($opt_seq, $aa_byname{$aa}->{location_list}->[$rand_pos1], 3) eq
      substr($opt_seq, $aa_byname{$aa}->{location_list}->[$rand_pos2], 3))
    { 
      $rand_pos2 = int(rand $num_of_codons_for_this_aa);
    } 
    my $location1 = $aa_byname{$aa}->{location_list}->[$rand_pos1];
    my $location2 = $aa_byname{$aa}->{location_list}->[$rand_pos2];
    my $codon1 = substr($opt_seq, $location1, 3);
    my $codon2 = substr($opt_seq, $location2, 3);
    # now check the scores now around the two codons and the scores if they are exchanged
    my $total_score_before = 0;
    my $total_score_after = 0;
    my ($previous_codon1, $previous_codon2, $next_codon1, $next_codon2) = ("", "", "", "");
    if (!$start_border{$location1})
    {
      $previous_codon1 = substr($opt_seq, $location1-3, 3);
      $total_score_before += (exists $codon_pair{$previous_codon1.$codon1}) ? $codon_pair{$previous_codon1.$codon1}->{value} : 0;
      $total_score_after += (exists $codon_pair{$previous_codon1.$codon2}) ? $codon_pair{$previous_codon1.$codon2}->{value} : 0;
    }
    if (!$end_border{$location1})
    {
      $next_codon1 = substr($opt_seq, $location1+3, 3);
      $total_score_before += (exists $codon_pair{$codon1.$next_codon1}) ? $codon_pair{$codon1.$next_codon1}->{value} : 0;
      $total_score_after += (exists $codon_pair{$codon2.$next_codon1}) ? $codon_pair{$codon2.$next_codon1}->{value} : 0;
    }
    if (!$start_border{$location2})
    {
      $previous_codon2 = substr($opt_seq, $location2-3, 3);
      $total_score_before += (exists $codon_pair{$previous_codon2.$codon2}) ? $codon_pair{$previous_codon2.$codon2}->{value} : 0;
      $total_score_after += (exists $codon_pair{$previous_codon2.$codon1}) ? $codon_pair{$previous_codon2.$codon1}->{value} : 0;
    }
    if (!$end_border{$location2})
    {
      $next_codon2 = substr($opt_seq, $location2+3, 3);
      $total_score_before += (exists $codon_pair{$codon2.$next_codon2}) ? $codon_pair{$codon2.$next_codon2}->{value} : 0;
      $total_score_after += (exists $codon_pair{$codon1.$next_codon2}) ? $codon_pair{$codon1.$next_codon2}->{value} : 0;
    }
    if ($max_condition)
    {
      if (($total_score_after > $total_score_before) || &sa($total_score_before - $total_score_after))
      {
        # Exchange
	substr($opt_seq, $location1, 3) = $codon2;
	substr($opt_seq, $location2, 3) = $codon1;
	my $score_diff = $total_score_after - $total_score_before;
	$opt_score += $score_diff;
	print OPT_FILE "$i\t$opt_score\n";
	$sa_acc_tran = 1;
      }
    }
    else
    {
      if (($total_score_after < $total_score_before) || &sa($total_score_after - $total_score_before))
      {
        # Exchange
	substr($opt_seq, $location1, 3) = $codon2;
	substr($opt_seq, $location2, 3) = $codon1;
	my $score_diff = $total_score_after - $total_score_before;
	$opt_score += $score_diff;
	print OPT_FILE "$i\t$opt_score\n";
	$sa_acc_tran = 1;
      }
    }
  }
  close(OPT_FILE);
  print ">seq.$opt_score\n";
  print $opt_seq, "\n";
  return $opt_score;
}

#%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
# This function returns yes if we are to accept a "negative" transition
# (opposite to our goals) and no if not. It makes the decision on the 
# value of the simulated annealing acceptance criteria
#%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
sub sa
{
  if ($args{s})
  {
    my ($sa_diff) = @_;
    my $sa_prob = rand;
    my $sa_val = exp(-$sa_diff/($sa_k*$sa_t));
    $sa_i++;
    if ($sa_i == $sa_iter)
    {
      $sa_i = 0;
      $sa_t *= $sa_a;
      if (!$sa_acc_tran)
      {
        $sa_exit = 1;
      }
    }
    if ($sa_val >= $sa_prob) {return 1;}
      else {return 0;}
  }
  else {return 0};
}

#%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
# Function that imports the codon information from the codons file and stores
# it in a hash.
#%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
sub import_codon_info
{
  open (CODON_FILE, "<data/codons")
    or die "Could not open codon info file data/codons for reading!\n";
  %aa = ();
  while(<CODON_FILE>)
  {
    chop($_);
    /^([^\s]+)\s+([^\s]+)\s+([^\s]+)\s+([^\s]+)/;
    push (@{$aa{$2}}, $1);
    $codon_to_aa{$1} = $2;
    $aa_code{$2} = $3;
    $code_to_aa{$3} = $2;
  }
  close(CODON_FILE);
  $variable_aas = 0;	# AAs that have more than one codon representations
  foreach $key (keys %aa)
  {
    push (@aas, $key);
    $aa_codon_num{$key} = 0;
    foreach $codon (@{$aa{$key}})
    {
      $codon_type{$codon} = $aa_codon_num{$key};
      $type_to_codon{$key}->[$aa_codon_num{$key}] = $codon;
      $aa_codon_num{$key}++;
    }
    if (($aa_codon_num{$key} > 1) && ($key ne "Ter"))
    {
      $variable_aas++;
    }
  }
  if ($DEBUG_CODON_IMPORT)
  {
    foreach $key (keys %aa)
    {
      print $key, "\t", $aa_code{$key}, "\n";
      foreach $codon (@{$aa{$key}})
      {
         print "\t$codon";
	 print "\t", $codon_to_aa{$codon};
	 print "\t", $codon_type{$codon}, "\n";
      }
      print "\n";
      print "Codons of $key: ";
      for (my $i = 0; $i < $aa_codon_num{$key}; $i++)
      {
        print $type_to_codon{$key}->[$i], " ";
      }
      print "\n";
    }
  }
}

#%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
# Function that imports the locks from the lock file $args{l}
#%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
sub import_locks
{
  open (LOCK_FILE, $args{l})
    or die "Could not open lock_file ", $args{l}, " for reading!\n";
  while(<LOCK_FILE>)
  {
    chop($_);
    /([^-]+)-([^-]+)$/;
    if (($1+0) > ($2+0))
    {
      die "\nERROR: Lock $1-$2 is reversely indicated. Please correct.\n\n";
    }
    ($locks[$num_of_locks][0], $locks[$num_of_locks++][1]) = ($1-1, $2-1);
  }
  close(LOCK_FILE);
  # Now create a mask table for all locations of the genome that are locked.
  for (my $i = 0; $i < $num_of_locks; $i++)
  {
    $lock_start{$locks[$i][0]} = 1;
    for ($j = $locks[$i][0]; $j <= $locks[$i][1]; $j++)
    {
      $lock_mask[$j] = 0;
    }
  }
}

#%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
# Function that imports the coding regions from the coding region file $args{c}
#%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
sub import_coding_regions
{
  open (CODING_REGION_FILE, $args{c})
    or die "Could not open coding regions file ", $args{c}, " for reading!\n";
  $all_region_seq_length = 0;
  while(<CODING_REGION_FILE>)
  {
    chop($_);
    /([^-]+)-([^-]+)$/;
    if ((($2 + 1 - $1) % 3) != 0)
    {
      die "\n ERROR: Coding region $1-$2 not a multiple of 3!\n\n";
    }
    if (($1+0) > ($2+0))
    {
      die "\nERROR: Coding region $1-$2 is reversely indicated. Please correct.\n\n";
    }
    ($coding_regions[$num_of_coding_regions][0], $coding_regions[$num_of_coding_regions][1]) = ($1-1, $2-1);
    my $reg_length = $coding_regions[$num_of_coding_regions][1] - $coding_regions[$num_of_coding_regions][0] + 1;
    $region_start[$num_of_coding_regions] = $all_region_seq_length;
    $all_region_seq_length += $reg_length;
    $region_end[$num_of_coding_regions] = $all_region_seq_length - 1;
    $region_length[$num_of_coding_regions++] = $region_length;
  }
  close(CODING_REGION_FILE);
}

#%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
# Function that returns the coding region a random position belongs to
# Arguments:
# 1. Random position (integer)
#%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
sub coding_region_of_position
{
  my $pos = $_[0];
  my $region = 0;
  while ($pos > $region_end[$region])
  {
    $region++;
  }
  return ($region);
}

#%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
# Function that imports the restriction sites
#%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
sub import_restriction_sites
{
  open (RESTRICTION_SITE_FILE, $args{r})
    or die "Could not open restriction site file ", $args{r}, " for reading!\n";
  while (<RESTRICTION_SITE_FILE>)
  {
    chop($_);
    if (/([^\s]+)/)
    {
      $restriction_site[$num_of_restriction_sites++] = $1;
    }
  }
  if ($DEBUG_RESTRICTION_SITES)
  {
    print "\nRestriction sites imported\n";
    for (my $i = 0; $i < $num_of_restriction_sites; $i++)
    {
      print $restriction_site[$i], "\n";
    }
    print "\n";
  }
}

#%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
# Function that isolates the codons we will work with. Basically creates
# a table with pointers to all codons we can change in the coding regions
# specified, after excluding the locked regions.
#%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
sub create_work_space
{
  $locked_codons_count = 0;
  for (my $i = 0; $i < $num_of_coding_regions; $i++)
  {
    for (my $j = $coding_regions[$i][0]; $j < $coding_regions[$i][1]; $j +=3)
    {
      my $codon = substr($fasta_record->{contig_seq}->[0], $j, 3);
      my $codon_t = $codon_type{$codon};
      if ($DEBUG_WORKSPACE) {print $codon, "\n";}
      if ($lock_mask[$j] && $lock_mask[$j+1] && $lock_mask[$j+2])
      {
        my $aa = $codon_to_aa{$codon};
        $aa_byname{$aa}->{num_of_codons}++;
        push (@{$aa_byname{$aa}->{location_list}}, $j);
        push (@{$aa_byname{$aa}->{location_type}}, $codon_t);
      }
      else
      {
        $locked_codons_count++;
	if ($DEBUG_WORKSPACE_LOCK)
	{
	  print "Codon $codon at position $j is locked\n";
	}
      }
    }
  }
  foreach my $aa (@aas)
  {
    if (!$aa_byname{$aa}->{num_of_codons})
    {
      $aa_byname{$aa}->{num_of_codons} = 0;
    }
  }
  if ($DEBUG_WORKSPACE_TOTAL)
  {
    foreach my $aa (@aas)
    {
      print "AminoAcid $key: ", $aa_byname{$aa}->{num_of_codons}, "\n";
    }
  }
}

#%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
# Function that calculates the codon distribution for the input sequence
#%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
sub calculate_source_distribution
{
  foreach my $a (@aas)
  {
    $number_of_codons = $aa_byname{$a}->{num_of_codons};
    $number_of_types = $aa_codon_num{$a};
    for (my $i = 0; $i < $number_of_types; $i++)
    {
      $type_mult[$i] = 0;
    }
    foreach my $t (@{$aa_byname{$a}->{location_type}})
    {
      $type_mult[$t]++;
    }
    my $how_many_different = 0;
    for (my $i = 0; $i < $number_of_types; $i++)
    {
      $aa_byname{$a}->{num_per_type}->[$i] = $type_mult[$i];
      if ($type_mult[$i] > 0) { $how_many_different++; }
    }
    $aa_byname{$a}->{at_least_two_different} = ($how_many_different) ? $how_many_different - 1 : 0;
    foreach my $codon (@{$aa{$a}})
    {
      my $type = $codon_type{$codon};
      if ($number_of_codons >0)
      {
        $source_dist->{$codon} = $type_mult[$type]/$number_of_codons;
      }
      else
      {
        $source_dist->{$codon} = 0;
      }
    }
    if ($DEBUG_SOURCE_DIST)
    {
      print "Amino Acid $a\t";
      foreach my $codon (@{$aa{$a}})
      {
        print $source_dist->{$codon}, " ";
      }
      print "\n";
    }
  }
}

#%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
# Function that imports a codon distribution from the file $args{d}
#%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
sub import_distribution
{
  open (DIST_FILE, $args{d})
    or die "Could not open codon info file ", $args{l}, " for reading: $!\n";
  $target_dist = {};
  while(<DIST_FILE>)
  {
    chop($_);
    if (/^([^\s]+)\s+([^\s]+)\s+([^\s]+)\s+([^\s]+)\s+([^\s]+)/)
    {
      $target_dist->{$2} = $5+0.0;
    }
  }
  close(DIST_FILE);
  if ($DEBUG_IMPORT_DIST)
  {
    foreach $key (keys %$target_dist)
    {
      print $key, "\t", $target_dist->{$key}, "\n";
    }
  }
}

#%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
# Search for each individual restriction site and eliminate it, with a synonymous
# change that does not fall in a locked region.
# This is an ad-hoc routine that does not do anything clever, so it should
# probably be substituted later for a more sophisticated one. For the time
# being it serves its purpose...
#%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
sub eliminate_restriction_sites
{
  print LOG_FILE "\n-----------------------------------------------------\n";
  print LOG_FILE "| Entering restriction site elimination phase       |\n";
  print LOG_FILE "-----------------------------------------------------\n\n";

  $base{"A"} = "A";
  $base{"C"} = "C";
  $base{"G"} = "G";
  $base{"T"} = "T";
  $base{"Y"} = "(C|T)";
  $base{"R"} = "(A|G)";
  $base{"N"} = "(A|C|G|T)";
  $sol{"A"} = ["C", "G", "T"];
  $sol{"C"} = ["A", "G", "T"];
  $sol{"G"} = ["A", "C", "T"];
  $sol{"T"} = ["A", "C", "G"];
  $sol{"Y"} = ["A", "G"];
  $sol{"R"} = ["C", "T"];
  $sol{"N"} = [];

  my @pattern = ();
  for (my $i = 0; $i < $num_of_restriction_sites; $i++)
  {
    $pattern[$i] = "";
    for (my $j = 0; $j < length($restriction_site[$i]); $j++)
    {
      $pattern[$i] .= $base{substr($restriction_site[$i], $j, 1)};
    }
    if ($DEBUG_PATTERN)
    {
      print "Restriction site ", $restriction_site[$i], " pattern: ", $pattern[$i], "\n";
    }
  }
  
  if ($DEBUG_PATTERN_MATCH)
  {
    my $seq = "ATTTCTAGGTCAG";
    print "sequence = $seq\n";
    my $pattern = $base{"A"}.$base{"R"}.$base{"N"}.$base{"Y"};
    print "pattern = $pattern\n";
    if ($seq =~ m/$pattern/g)
    {
      print "pattern found at position ", pos($seq), "\n";
      print "$` $& $'\n";
    }
    die;
  }
  
#  $output_sequence = $shuffled_sequence;
  $output_sequence = $fasta_record->{contig_seq}->[0];;
  for (my $r = 0; $r < $num_of_coding_regions; $r++)
  {
    my $modulo = $coding_regions[$r][0] % 3;
    if ($DEBUG_MODULO)
    {
      print "\nCoding region $r: Start position ", $coding_regions[$r][0], " mod 3 = $modulo\n\n";
    }
    my $satisfied = 0;
    # Until we do not encounter another restriction site that is not locked
    $elimination_rounds = 0;
    until ($satisfied)
    {
      $satisfied = 1;
      for (my $i = 0; $i < $num_of_restriction_sites; $i++)
      {
	my $rs_len = length($restriction_site[$i]);
	if ($DEBUG_ELIMINATE_RS) {print "Examining restriction site ", $restriction_site[$i], " with length $rs_len\n"}
        my $p = $pattern[$i];
        pos($output_sequence) = $coding_regions[$r][0];
        while (($output_sequence =~ /$p/g) && (pos($output_sequence) <= $coding_regions[$r][1]))
	{
	  my $rs_start = pos($output_sequence) - $rs_len;
	  if ($DEBUG_ELIMINATE_RS) {print "--Found at position $rs_start\n"}
	  print LOG_FILE "Detected site ", $restriction_site[$i], " at position ", $rs_start, "\n";
	  # Check if this is a lock we want to preserve. If not, continue
	  if (!exists $lock_start{$rs_start})
	  {
	    # Oops, found another restriction site.
	    $satisfied = 0;
	    # For each of the codons that intersect the restriction site area, determine
	    # if they can be changed with a synonymous one that would work. Three conditions:
	    # 1. Not a single amino-acid codon
	    # 2. Change in area outside the intersection
	    # 3. Change affects a polycharacter base of the restriction site that will not make a difference.
	    # $pos_diff is the number of bases before the start of the restriction site that a codon starts
	    my $pos_diff = ($rs_start - $coding_regions[$r][0]) % 3;
	    my $count = $rs_start - $pos_diff;
	    my $found = 0;
	    if ($DEBUG_ELIMINATE_RS) {print "  First intersecting codon at position $count ($pos_diff difference)\n"}
	    EXIT1: while (($count + 3 < pos($output_sequence)) && (!$found))
	    {
	      my $cdn = substr($output_sequence, $count, 3);
	      my $a = $codon_to_aa{$cdn};
	      if ($DEBUG_ELIMINATE_RS) {print "    Examining to change codon $cdn (amino acid $a)\n"}
	      # if this amino acid has more than one codons...
	      if ($aa_codon_num{$a} > 1)
	      {
	        EXIT2: for (my $j = 0; $j < $aa_codon_num{$a}; $j++)
		{
		  $cdn2 = $aa{$a}->[$j];
	          if ($DEBUG_ELIMINATE_RS) {print "      Maybe substitute it with $cdn2?\n"}
		  # for each codon of this amino acid that is not the one we have
		  if ($cdn2 ne $cdn)
		  {
		    EXIT3: for (my $k = 0; $k < 3; $k++)
		    {
		      # if they differ in some character, check if that character is inside the restriction
		      # site area and whether is is a legal character to substitute (in case the restriction
		      # site contains multicharacters
		      my $current_pos = $count+$k;
		      if ($DEBUG_ELIMINATE_RS) {print "Checking whether ", substr($cdn, $k, 1), " ne ", substr($cdn2, $k, 1), " AND $current_pos >= $rs_start AND $current_pos < ", pos($output_sequence), "\n"}
		      if ((substr($cdn, $k, 1) ne substr($cdn2, $k, 1)) && \
		        ($current_pos >= $rs_start) && ($current_pos < pos($output_sequence)))
		      {
			my $current_char = substr($restriction_site[$i], $current_pos - $rs_start, 1);
	                if ($DEBUG_ELIMINATE_RS) {print "        Examining character $current_char at position $current_pos\n"}
			my $proposed_char = substr($cdn2, $k, 1);
	                if ($DEBUG_ELIMINATE_RS) {print "          Proposed substitution: $proposed_char\n"}
			EXIT4: for (my $l = 0; $l <= $#{$sol{$current_char}}; $l++)
			{
	                  if ($DEBUG_ELIMINATE_RS) {print "            Checking against: ", $sol{$current_char}->[$l], "\n"}
			  if ($proposed_char eq $sol{$current_char}->[$l])
			  {
			    substr($output_sequence, $current_pos, 1) = $proposed_char;
			    $found = 1;
			    print LOG_FILE "  Eliminated by substituting $current_char with $proposed_char at position $current_pos\n";
	                    if ($DEBUG_ELIMINATE_RS) {print "  Eliminated by substituting $current_char with $proposed_char at position $current_pos\n"}
			    last EXIT1;
			    last EXIT2;
			    last EXIT3;
			    last EXIT4;
			  }
			}
		      }
		    }
		  }
		}
	      }
	      $count += 3;
	    }
	    if (!$found)
	    {
	      print LOG_FILE "  Cannot eliminate restriction site ", $restriction_site[$i], " at position $rs_start\n\n";
	      if ($DEBUG_ELIMINATE_RS) {print "  Cannot eliminate restriction site ", $restriction_site[$i], " at position $rs_start\n"};
	    }
	  }
	  else
	  {
	    print LOG_FILE "  Site locked\n";
	    if ($DEBUG_ELIMINATE_RS) {print "  Site Locked!\n"};
	  }
	  pos($output_sequence) = $rs_start + 1;
	}
      }
      $elimination_rounds++;
      print LOG_FILE "**************** End of round $elimination_rounds\n";
      if ($DEBUG_ELIMINATE_RS) {print "One more elimination round? ($elimination_rounds)\n"};
    }
  }
  print LOG_FILE "Total elimination rounds = $elimination_rounds\n";
}

#%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
# Create an index permutation. Given the index size, return a pointer to an 
# array that has the indices permuted.
# Arguments:
# 1. Size of index (integer)
#%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
sub permutation
{
  my @index = ();
  my @temp_index = ();
  my ($size) = @_;
  my @rand_index = ();
  for (my $i = 0; $i < $size; $i++)
  {
    push (@temp_index, $i);
  }
  for ($i = 0; $i < $size; $i++) 
  {
    $rand_index[$i] = rand;
  }
  @index = sort {$rand_index[$a] <=> $rand_index[$b]} @temp_index;
  for ($i = 0; $i < ($size - 1); $i++) 
  {
    my $random = int(rand($size-$i)+$i);
    ($index[$i], $index[$random]) = ($index[$random], $index[$i]);
  }
  my @rev_index = ();
  for ($i = 0; $i < $size; $i++)
  {
    $rev_index[$index[$i]] = $i;
  }
  if ($DEBUG_PERMUTATION)
  {
    for ($i = 0; $i < $size; $i++)
    {
      print $index[$i], " ", $rev_index[$i], "\n";
    }
  }
  return (\@index, \@rev_index);
}
