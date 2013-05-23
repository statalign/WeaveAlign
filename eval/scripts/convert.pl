#!/usr/bin/perl
#
# Converts file to given format if it's in a different format.
# Doesn't do anything if target file exists.
#
# If EXTA and EXTB are the format extensions the module
# named EXTA2EXTB.pl is employed to carry out the conversion.
#
# Input: file.EXTA EXTB
#
# Output:
#  out        output file name (file.EXTB)
#  converted  1 if conversion took place 0 otherwise
#

use warnings;

$path = $0 =~ m%(.*[/\\])% ? $1 : "";
require "${path}module.pl";

sub convert {

  my $in = shift;
  my $extb = lc shift or print(STDERR "convert.pl: target extension not specified\n"), return{};
  my $out = $in;
  -e $in or print(STDERR "convert.pl: input file $in does not exist\n"), return{};
  $in =~ /([^.]*)$/;
  my $exta = lc $1;
  $out =~ s/$exta$/$extb/;
  my $conv=(-e $out)?0:1;
  if($conv) {
    my $mod = "${exta}2${extb}";
    require "${path}$mod.pl";
    ${&$mod($in,$out)}{'converted'} or return{};
  }

  {out=>$out,converted=>$conv};
}

1;