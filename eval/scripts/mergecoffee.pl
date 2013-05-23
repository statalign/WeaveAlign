#!/usr/bin/env perl
#
# Runs MergeAlign on alignments created by several other methods.
# Alignments must be in Fasta format.
#
# Methods currently used are: clustal, dialign, mafft, muscle, pcma, poa,
# probcons, tcoffee. Methods are assumed to have been pre-run.
#
# Input: [-multicore] refalign.ext [reftree]
#  alignment created by method mmm must be in refalign.mmm.fsa
#  output is refalign.mergecoffee.fsa
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

use strict;
use warnings;

my $path = $0 =~ m%(.*[/\\])% ? $1 : "";
require "${path}module.pl";

sub mergecoffee {

  my @methods = qw(clustal dialign mafft muscle pcma poa probcons tcoffee);

  my $ref = shift;
  my $reftree = shift;

  $ref =~ /(.*)\./, my $base = $1;
  my @flist = grep -e $_, map "$base.$_.fsa", @methods;
  @flist > 1 or print(STDERR "mergecoffee.pl: did not find alignments to combine!\n"), return{};
  my $out = "$base.mergecoffee.fsa";

  require "${path}mergecons.pl";
  mergecons($ref, @flist, $reftree?$reftree:(), "-out=$out");

}

1;
