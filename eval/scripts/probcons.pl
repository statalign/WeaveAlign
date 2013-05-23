#!/usr/bin/env perl
#
# Runs ProbConsRNA on nucleotide sequences. File is transformed into Fasta
# format first, see convert.pl.
#
# Input: file.ext [reftree]
#  output is file.probcons.fsa
#
# Output:
#  _infile       input file (after conversion)
#  [scores]      score.pl's each output field
#  ztime         elapsed time (seconds)
#

use warnings;

$path = $0 =~ m%(.*[/\\])% ? $1 : "";
require "${path}module.pl";

sub probcons {

  my $in = shift;
  my $reftree = shift;

  require "${path}convert.pl";
  $in = ${convert($in,'fsa')}{'out'} or return{};
  (my $out=$in) =~ s/fsa$/probcons.fsa/;
  require "${path}rmgaps.pl";
  my ($in_gaprm,$in_conv) = @{rmgaps($in,'fasta')}{'out','converted'};

  my $time = -time;
  $_=`probconsRNA $in_gaprm 2>&1 >$out`;
  print(),return{} if $?;
  $time += time;

  unlink $in_gaprm if $in_conv;

  require "${path}scoreall.pl";
  my $score = scoreall($in,$out,$reftree//());

  { _infile=>$in,
    %$score,
    ztime=>$time };
}

1;
