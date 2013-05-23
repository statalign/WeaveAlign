#!/usr/bin/perl
#
# Converts MPD file to FASTA format. No checks are performed. As internal alignment
# format it supports both Fasta and StatAlign.
#
# Input: infile.mpd outfile.fsa
#
# Output:
#  converted  1
#

use warnings;

$path = $0 =~ m%(.*[/\\])% ? $1 : "";
require "${path}module.pl";

sub mpd2fsa {

  my $in = shift or print(STDERR "mpd2fsa.pl: input file not specified\n"), return{};
  my $out = shift or print(STDERR "mpd2fsa.pl: output file not specified\n"), return{};

  open IN, $in;
  open OUT, ">$out";

  while(<IN>) {
    last if /^#scores/;
    @A = split "\t";
    print(OUT), next if @A != 2;
    print OUT ">$A[0]\n$A[1]";
  }

  close IN;
  close OUT;

  {'converted' => 1};
}

1;