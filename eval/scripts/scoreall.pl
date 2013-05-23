#!/usr/bin/env perl
#
# Tests accuracy of alignment using all available scores, including scores
# weighted using a tree and ML tree accuracy.
#
# Input: refalign.ext testalign.ext [reftree] [-key=value...]
#  falls back to simply calling score.pl when reftree is missing
#  optional key=value pairs are added to the output
#
# Output:
#  [aliscores]   score.pl's output fields
#  [waliscores]  weighted scores from score.pl, prefixed with a 'w'
#  [treescores]  tscore.pl's output fields, prefixed with a 'x'
#

use warnings;

$path = $0 =~ m%(.*[/\\])% ? $1 : "";
require "${path}module.pl";

sub scoreall {

  my $ref = shift;
  my $test = shift or print(STDERR "scoreall.pl: required parameters missing"), return{};
  my $tree; $tree = shift if @_ && $_[0] !~ /^-/;

  my %out = map split('=',substr $_,1), @_;

  require "${path}score.pl";
  my $scores = score($ref, $test);
  if($tree) {
    my $wscores = score($ref, $test, $tree);
    $scores->{"w$_"} = $wscores->{$_} for grep !/^[as_]/, keys %$wscores;

    require "${path}tscore.pl";
    my $tscores = tscore($test, $tree);
    unlink $tscores->{_outfile};  # do not keep tree
    $scores->{"x$_"} = $tscores->{$_} for grep !/^_/, keys %$tscores;
  }
  $scores->{$_} = $out{$_} for keys %out;

  $scores;
}

1;
