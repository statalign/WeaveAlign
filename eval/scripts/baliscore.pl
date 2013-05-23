#!/usr/bin/env perl
#
# Gives score of an alignment with respect to a reference alignment by
# either of baliscore2.pl or baliscore3.pl depending on whether an
# annotation file is present.
#
# Converts both alignments into MSF format first (temporarily) if
# they're in a different format, see convert.pl. If conversion took
# place temporary files are removed afterwards.
#
# Input: reffile.ext infile.ext
#  annotation file is assumed to be the first of reffile*.ann
#
# Output:
#  spscore    SP score
#  tcscore    TC score
#  fullsp     full SP score
#  fulltc     full TC score
#

use warnings;

$path = $0 =~ m%(.*[/\\])% ? $1 : "";
require "${path}module.pl";

sub baliscore {

  my $ref = shift;
  my $in = shift;

  my $base = $ref;
  $base =~ s/.[^.]*$//;
  my $ann = (<${base}*.ann>)[0] || "";

  require "${path}baliscore2.pl";
  require "${path}baliscore3.pl";
  -e $ann ? baliscore2($ref,$in) : baliscore3($ref,$in);
  #baliscore3($ref,$in);
}

1;
