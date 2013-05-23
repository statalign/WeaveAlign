#!/usr/bin/env perl
#
# Runs T-Coffee on an alignment file. File is transformed to Fasta
# format first, see convert.pl.
#
# Input: [-multicore] [-keeptree] file.ext [reftree]
#  output is file.tcoffee.fsa
#  multicore mode is activated when -multicore is used
#  guide tree is written into file.tcoffee.tree when -keeptree used
#
# Output:
#  _infile       input file (after conversion)
#  [scores]      score.pl's each output field
#  ztime         elapsed time (seconds)
#
##unused
##  treeout       Tree file created or empty
##  outfsa        Output Fasta file name
#

use strict;
use warnings;

my $path = $0 =~ m%(.*[/\\])% ? $1 : "";
require "${path}module.pl";

sub tcoffee {

  my $temptree="tcoffee.tmp.tree";
  my $keeptree=0;
  my $singlecore='-multi_core=no';

  my $in = shift;
  $in = shift, $singlecore = '' if $in eq '-multicore';
  $in = shift, $keeptree = 1 if $in eq '-keeptree';
  my $reftree = shift;

  require "${path}convert.pl";
  $in = ${convert($in,'fsa')}{'out'} or return{};
  (my $out=$in) =~ s/fsa$/tcoffee.fsa/;
  require "${path}rmgaps.pl";
  my ($in_gaprm,$in_conv) = @{rmgaps($in,'fasta')}{'out','converted'};

  my $time = -time;
  $_=`t_coffee $in_gaprm -output=fasta -newtree=$temptree $singlecore -outfile=$out 2>&1`;
  print(),return{} if $?;
  $time += time;

  my $treeout = '';
  if($keeptree) {
    require "${path}treenorm.pl";
    ($treeout = $in) =~ s/fsa$/tcoffee.tree/;
    ${treenorm($temptree,$treeout)}{'converted'} or return{};
  }
  unlink $temptree;
  unlink $in_gaprm if $in_conv;

  require "${path}scoreall.pl";
  my $score = scoreall($in,$out,$reftree//());

  { _infile=>$in,
    %$score,
    ztime=>$time };#, treeout=>$treeout };

}

1;
