#!/usr/bin/env perl
#
# Gives score of an alignment by bali_score with respect to a reference
# alignment. Converts both into MSF format first (temporarily) if they're
# in a different format, see convert.pl.
#
# If conversion took place temporary files are removed afterwards.
#
# Input: reffile infile [basename.ext]
#  basename.ext defaults to reffile
#  annotation file is assumed to be the first of basename*.ann
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

sub baliscore2 {

  my $ref = shift;
  my $in = shift;
  my $base = shift;

  $base ||= $ref;
  $base =~ s/.[^.]*$//;
  my $ann = (<${base}*.ann>)[0] || "";
  -e $ann or print(STDERR "baliscore2.pl: annotation file $ann does not exist\n"), $ann='';

  require "${path}convert.pl";
  my ($rconv,$iconv);
  ($ref,$rconv) = @{convert($ref,'msf')}{'out','converted'}, $ref or return{};
  ($in, $iconv) = @{convert($in, 'msf')}{'out','converted'}, $in  or return{};


  $_=`bali_score $ref $in $ann`;
  print(), return{} if $?;
  my ($sp,$tc) = /score= (.*)/g;

  $_=`bali_score $ref $in`;
  print(), return{} if $?;
  my ($fullsp,$fulltc) = /score= (.*)/g;

  unlink $ref if $rconv;
  unlink $in if $iconv;

  {spscore=>$sp, tcscore=>$tc, fullsp=>$fullsp, fulltc=>$fulltc};
}

1;
