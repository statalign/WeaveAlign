#!/usr/bin/perl
#
# Generates aligned sequences consisting of fast-slow blocks for testing
# footprinting methods.
#
# (C) novadam, 2011.
#

use strict;
use warnings;

#@ARGV == 1 or print("use: eval_dawg_bonphy.pl comp_size"), exit;

my $seqno = 10;
my $first = 91;
my $n = 10;

#my $segments = 9;  # number of fast-slow segments total
my $segments = 1;  # number of fast-slow segments total
#my @fastlenrange = (20, 80);  # range of fast segment lengths (uniform)
my @fastlenrange = (800, 1000);  # range of fast segment lengths (uniform)
my @slowlenrange = (10, 18);  # range of slow segment lengths (uniform)

#my $fastscale = 1.0;  # edge scaling ratio in slow segments
my $fastscale = 1;  # edge scaling ratio in slow segments
my $slowscale = 0.3;  # edge scaling ratio in slow segments

#my $fastscale = 2.0;  # edge scaling ratio in slow segments
#my $slowscale = 0.3;  # edge scaling ratio in slow segments
#my $fastscale = 5.0;  # edge scaling ratio in slow segments
#my $slowscale = 0.05;  # edge scaling ratio in slow segments

my $templatef = "dawg.templ";
my $dawgin = "dawg_in.cfg";
my $dawgout = "dawg_out.fasta";
my $alignfile = "data.fsa";
my $refmotifs = "refmotifs";

for my $i ($first..($first+$n-1)) {
  my $dir = sprintf "t%03d", $i;
  -e $dir and `rm -r $dir`;
  mkdir $dir or die "create dir failed";

#  `perl randtree.pl $seqno >$dir/data.tree`;
  `cp data.tree $dir/data.tree`;

  my $dawg = `perl subst.pl $templatef *TREE* $dir/data.tree`;

  my %seq;
  my @refmots;

  for my $seg (1..$segments) {
    my @lenrange = (($seg & 1) ? @fastlenrange : @slowlenrange);
    my $scale = (($seg & 1) ? $fastscale : $slowscale);

    (my $dawgf = $dawg) =~ s/\*SCALE\*/$scale/g;
    my $len = int(rand($lenrange[1]-$lenrange[0]+1))+$lenrange[0];
    $dawgf =~ s/\*LEN\*/$len/;
    save($dawgin, $dawgf);

    `./dawg $dawgin`;

    my $flen;
    open F, $dawgout;
    my $name;
    while(<F>) {
      chomp;
      /^> *(.*)/ and $name=$1, $flen = 0, next;
      $seq{$name} .= $_;
      $flen += length($_);
    }
    close F;

    push @refmots, (1-($seg&1)) x $flen;
  }
  unlink $dawgin;
  unlink $dawgout;

  open F, ">$dir/$alignfile";
  print F ">$_\n$seq{$_}\n" for keys %seq;
  close F;

  open F, ">$dir/refmotifs";
  print F "$_\n" for @refmots;
  close F;
}

sub save {
  open F, ">$_[0]";
  print F $_[1];
  close F;
}

#sub load {
#  open F, $_[0];
#  my $lines = join '', <F>;
#  close F;
#  $lines;
#}