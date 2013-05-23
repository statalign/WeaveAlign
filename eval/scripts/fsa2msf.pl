#!/usr/bin/perl
#
# Converts FASTA file to MSF format. No checks are performed.
#
# Input: infile.fsa outfile.msf
#
# Output:
#  converted  1
#

use warnings;

$path = $0 =~ m%(.*[/\\])% ? $1 : "";
require "${path}module.pl";

sub fsa2msf_max {
 $_[0] > $_[1] ? $_[0] : $_[1];
}

sub fsa2msf {

  my $in = shift or print(STDERR "fsa2msf.pl: input file not specified\n"), return{};
  my $out = shift or print(STDERR "fsa2msf.pl: output file not specified\n"), return{};

  open IN, $in;
  open OUT, ">$out";

  my $name;
  my %seq;

  while(<IN>) {
   chomp;
   s/\r//g;
   /^> *(.*)/ and $name=$1, next;
   y/-/./;
   $seq{$name} .= $_;
  }

  $len=0;
  $len = fsa2msf_max($len,length($_)) for keys %seq;
  %seq = map {$_.' 'x($len-length($_)) => $seq{$_} } keys %seq;
  $len=0;
  $len = fsa2msf_max($len,length($seq{$_})) for keys %seq;
  $seq{$_} .= '.'x($len-length($seq{$_})) for keys %seq;

  print OUT "PileUp\n\n\n\n   MSF:  $len  Type: P    Check:  1   ..\n\n";
  print OUT " Name: $_ oo  Len:   $len  Check:  1  Weight:  1.00\n" for keys %seq;
  print OUT "\n//\n\n\n\n";

  for($i = 0; $i < $len; $i += 50) {
   for $name (keys %seq) {
    print OUT "$name     ";
    print OUT "$_ " for (substr $seq{$name}, $i, 50) =~ /.{1,10}/g;
    print OUT "\n";
   }
   print OUT "\n\n";
  }

  close IN;
  close OUT;

  {'converted' => 1};
}

1;
