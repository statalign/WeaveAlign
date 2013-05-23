#!/usr/bin/env perl
#
# Runs ClustalW on an alignment file. File is transformed to Fasta
# format first, see convert.pl.
#
# Input: [-keeptree] file.ext [reftree]
#  output is file.clustal.fsa
#
# Output:
#  _infile       input file (after conversion)
#  [scores]      score.pl's each output field
#  ztime         elapsed time (seconds)
#
##unused
##  outfsa        Output Fasta file name
#

use warnings;

$path = $0 =~ m%(.*[/\\])% ? $1 : "";
require "${path}module.pl";

sub clustal {

  my $keeptree = 0;

  my $in = shift;
  $in = shift, $keeptree = 1 if $in eq '-keeptree';
  my $reftree = shift;

  require "${path}convert.pl";
  $in = ${convert($in,'fsa')}{'out'} or return{};
  (my $out=$in) =~ s/fsa$/clustal.fsa/;
  require "${path}rmgaps.pl";
  my ($in_gaprm,$in_conv) = @{rmgaps($in,'fasta')}{'out','converted'};

  my $time = -time;
  $_=`clustalw2 -infile=$in_gaprm -align -output=FASTA -outfile=$out`;
  print(),return{} if $?;
  $time += time;

  my $dnd;
  ($dnd = $in) =~ s/fsa$/dnd/;
  unlink $dnd if -e $dnd and !$keeptree;
  unlink $in_gaprm if $in_conv;

  require "${path}scoreall.pl";
  my $score = scoreall($in,$out,$reftree//());

  { _infile=>$in,
    %$score,
    ztime=>$time };
}

1;
