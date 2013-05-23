#!/usr/bin/env perl
#
# Runs DIALIGN-TX in DNA mode on an alignment file. File is transformed to
# Fasta format first, see convert.pl.
#
# Input: file.ext [reftree]
#  output is written into file.dialign.fsa
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

sub dialign {

  my $in = shift or die "dialign.pl: required parameter missing";
  my $reftree = shift;

  require "${path}convert.pl";
  $in = ${convert($in,'fsa')}{'out'} or return{};
  (my $out = $in) =~ s/fsa$/dialign.fsa/;
  require "${path}rmgaps.pl";
  my ($in_gaprm,$in_conv) = @{rmgaps($in,'fasta')}{'out','converted'};

  (my $dout = $in_gaprm) =~ s/fasta$/fa/;
  (my $dali = $in_gaprm) =~ s/fasta$/ali/;

  my $time = -time;
  $_=`dialign-tx -n -fa $in_gaprm`;
  print(),return{} if $?;
  $time += time;

  my $ali = slurp($dout);
  $ali =~ s/\r//g;
  open my $f, ">$out";
  print $f $ali;
  close $f;
  unlink $dout;
  unlink $dali if -e $dali;

  unlink $in_gaprm if $in_conv;

  require "${path}scoreall.pl";
  my $score = scoreall($in,$out,$reftree//());

  { _infile=>$in,
    %$score,
    ztime=>$time };#, treeout=>$treeout };

}

1;
