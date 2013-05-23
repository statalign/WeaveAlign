#!/usr/bin/perl
#
# Converts MSF file to FASTA format. No checks are performed.
#
# Input: infile.msf outfile.fsa
#
# Output:
#  converted  1
#

use warnings;

$path = $0 =~ m%(.*[/\\])% ? $1 : "";
require "${path}module.pl";

sub msf2fsa {

  my $in = shift or print(STDERR "msf2fsa.pl: input file not specified\n"), return{};
  my $out = shift or print(STDERR "msf2fsa.pl: output file not specified\n"), return{};

  open IN, $in;
  open OUT, ">$out";

  while(<IN>) {
    m#^//# and $seqstart=1, next;
    $seqstart or next;

    my ($name,$seq) = split / {3,}/;
    $seq or next;
    $seq =~ tr/. /-/d;

    # removing gaps
    # $seq =~ tr/-//d;

    $seqs{$name} .= $seq;
  }

  print OUT ">$_\n$seqs{$_}" for keys %seqs;

  close IN;
  close OUT;

  {'converted' => 1};
}

1;