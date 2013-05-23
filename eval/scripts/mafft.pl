#!/usr/bin/env perl
#
# Runs MAFFT on an alignment file. File is transformed to Fasta
# format first, see convert.pl.
#
# Input: [-keeptree] file.ext [reftree]
#  output is written into file.mafft.fsa
#  guide tree is written into file.mafft.tree when -keeptree used
#
# Output:
#  _infile       input file (after conversion)
#  [scores]      score.pl's each output field
#  ztime         elapsed time (seconds)
#
##unused
##  treeout       tree file created or empty
##  outfsa        output Fasta file name
#

use warnings;

$path = $0 =~ m%(.*[/\\])% ? $1 : "";
require "${path}module.pl";

sub mafft {

  $keeptree='';

  my $in = shift;
  $in = shift, $keeptree = ' --treeout' if $in eq '-keeptree';
  #$keeptree = ' --treeout';
  my $reftree = shift;

  require "${path}convert.pl";
  $in = ${convert($in,'fsa')}{'out'} or return{};
  (my $out = $in) =~ s/fsa$/mafft.fsa/;
  require "${path}rmgaps.pl";
  my ($in_gaprm,$in_conv) = @{rmgaps($in,'fasta')}{'out','converted'};

  my $time = -time;
  $_=`mafft --auto$keeptree $in_gaprm 2>&1 >$out`;
  print(),return{} if $?;
  $time += time;

  my $treeout = '';
  if($keeptree) {
    require "${path}treenorm.pl";
    ($treeout = $in) =~ s/fsa$/mafft.tree/;
    ${treenorm("$in_gaprm.tree",$treeout)}{'converted'} or return{};
    unlink "$in_gaprm.tree";
  }
  unlink $in_gaprm if $in_conv;

  require "${path}scoreall.pl";
  my $score = scoreall($in,$out,$reftree//());

  { _infile=>$in,
    %$score,
    ztime=>$time };#, treeout=>$treeout };

}

1;
