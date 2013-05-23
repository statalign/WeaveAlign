#!/usr/bin/env perl
#
# Runs MPD to combine several alignments created by other methods into
# one consensus alignment. Inputs must be in Fasta format.
#
# Methods currently used are: clustal, dialign, mafft, muscle, pcma, poa,
# probcons, tcoffee. Methods are assumed to have been pre-run.
#
# Input: refalign.ext [reftree] [-out=outfile.fsa] [params...]
#  alignment created by method mmm must be in refalign.mmm.fsa
#  default output is refalign.mpdcomb.fsa
#  optional params... are transferred to MPD directly
#
# Output:
#  _infile       reference alignment
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

sub mpdcomb {

  my @methods = qw(clustal dialign mafft muscle pcma poa probcons tcoffee);

  my $ref = shift;
  my $reftree;
  $reftree = shift if @_ && $_[0] !~ /^-/;

  $ref =~ /(.*)\./, my $base = $1;
  my @flist = grep -e $_, map "$base.$_.fsa", @methods;
  @flist > 1 or print(STDERR "mpdcomb.pl: did not find alignments to combine!\n"), return{};
  grep /^-out=/, @_ or push @_, "-out=$base.mpdcomb.fsa";

  require "${path}mpd.pl";
  mpd($ref, $reftree//(), @flist, @_);

}

1;
