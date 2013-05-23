#!/usr/bin/env perl
#
# Gives score of an alignment by bali_score3 with respect to a reference
# alignment. Converts both into MSF format first (temporarily) if they're
# in a different format, see convert.pl.
#
# If conversion took place temporary files are removed afterwards.
#
# Input: reffile.ext infile
#  annotation file is assumed to be reffile.xml
#
# Output:
#  spscore    SP score (0 if annotation file is not found)
#  tcscore    TC score (0 if annotation file is not found)
#  fullsp     full SP score
#  fulltc     full TC score
#

use warnings;

$path = $0 =~ m%(.*[/\\])% ? $1 : "";
require "${path}module.pl";

sub baliscore3 {

  my $ref = shift;
  my $in = shift;
  my $base = shift;

  require "${path}convert.pl";
  my ($rconv,$iconv);
  ($ref,$rconv) = @{convert($ref,'msf')}{'out','converted'}, $ref or return{};
  ($in, $iconv) = @{convert($in, 'msf')}{'out','converted'}, $in  or return{};

  $_=`bali_score3 $ref $in`;
  print(), return{} if $?;
  my ($fullsp,$fulltc) = /score= (.*)/g;

  my ($sp,$tc) = (0,0);
  my $ann = $ref;
  $ann =~ s/msf$/xml/;
  if(-e $ann) {
    $_=`bali_score3 $ann $in`;
    print(), return{} if $?;
    ($sp,$tc) = /score= (.*)/g;
  }

#  unlink $ref if $rconv;
#  unlink $in if $iconv;

  {spscore=>$sp, tcscore=>$tc, fullsp=>$fullsp, fulltc=>$fulltc};
}

1;
