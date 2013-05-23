#!/usr/bin/perl
#
# Removes gaps from a fasta file containing alignments.
# Doesn't do anything if target file exists.
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

sub rmgaps {

  my $in = shift;
  my $extb = lc shift;
  my $out = $in;
  -e $in or print(STDERR "convert.pl: input file $in does not exist\n"), return{};
  $in =~ /([^.]*)$/;
  my $exta = lc $1;
  $out =~ s/$exta$/$extb/;
  my $conv=(-e $out)?0:1;
  if($conv) {
    open IN, "$in";
    open OUT, ">$out";
    while(<IN>) {
      y/-//d if !/^>/;
      print OUT;
    }
    close IN;
    close OUT;
  }

  return {out=>$out,converted=>$conv};
}

1;