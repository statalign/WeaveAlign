#!/usr/bin/env perl
#
# Runs PCMA on an alignment file. File is transformed to Fasta
# format first, see convert.pl.
#
# Input: file.ext [reftree]
#  output is file.pcma.fsa
#
# Output:
#  _infile       input file (after conversion)
#  [scores]      score.pl's each output field
#  ztime         elapsed time (seconds)
#

use warnings;

$path = $0 =~ m%(.*[/\\])% ? $1 : "";
require "${path}module.pl";

sub pcma {

  my $in = shift;
  my $reftree = shift;

  require "${path}convert.pl";
  $in = ${convert($in,'fsa')}{'out'} or return{};
  (my $out=$in) =~ s/fsa$/pcma.fsa/;
  (my $pout=$out) =~ s/fsa$/aln/;
  require "${path}rmgaps.pl";
  my ($in_gaprm,$in_conv) = @{rmgaps($in,'fasta')}{'out','converted'};
  (my $dnd = $in_gaprm) =~ s/.fasta$/.dnd/;

  my $time = -time;
  $_=`pcma $in_gaprm -outfile=$pout 2>&1`;
  print(),return{} if $?;
  $time += time;

  unlink $dnd if -e $dnd;

  `clustalw2 -convert -infile=$pout -output=fasta -outfile=$out`;
  -e $out or print(STDERR "pcma.pl: converting MSA output to Fasta failed!\n"), return{};
  unlink $pout;

  unlink $in_gaprm if $in_conv;

  require "${path}scoreall.pl";
  my $score = scoreall($in,$out,$reftree//());

  { _infile=>$in,
    %$score,
    ztime=>$time };
}

1;
