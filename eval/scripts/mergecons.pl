#!/usr/bin/env perl
#
# Runs MergeAlign on a set of alignments. Files must all be
# FASTA files or a single StatAlign log file with samples.
#
# Input: refalign.ext input_1.fsa ... [reftree] [-out=out] OR
#        refalign.ext [input.log] [reftree] [-out=out]
#  default output is refalign.mergecons.fsa but can be overridden with -out
#  input.log is assumed to be refalign.ext.log when omitted
#
# Output:
#  _infile       input file (after conversion)
#  [scores]      score.pl's each output field
#  ztime         elapsed time (seconds)
#

use strict;
use warnings;
use File::Copy;

my $path = $0 =~ m%(.*[/\\])% ? $1 : "";
require "${path}module.pl";

sub mergecons_xlog {
  my $tempfiles = "mctemp.%d.fsa";

  my $in = shift;
  my $dir = shift;
  open my $f, $in or return;
  my $cur = -1;
  my @flist = ();
  my $g;
  while(<$f>) {
    if(/Sample (\d+)\tAlignment:\t(.*)/s) {
      if($cur != $1) {
        close $g if $cur != -1;
        open $g, ">".($_="$dir/".(sprintf $tempfiles, $1)) or return;
        push @flist, $_;
        $cur = $1;
      }
      print $g $2;
    }
  }
  close $g if $cur != -1;
  close $f;

  # removing gap-only cols
  for my $fn (@flist) {
    open $f, $fn or print(STDERR "mergecons.pl: cannot open $fn\n"), return;
    my (%seq,$name,@nlist);
    while(<$f>) {
      chomp;
      s/\r//g;
      /^> *(.*)/ and $name=$1, push(@nlist,$1), next;
      $seq{$name} .= $_;
    }
    close $f;
    my $len = length($seq{$name});
    my %seq2;
    for my $i (0..($len-1)) {
      if(grep substr($seq{$_},$i,1) ne '-', @nlist) {
         $seq2{$_} .= substr($seq{$_},$i,1) for @nlist;
      }
    }
    open $f, ">$fn" or print(STDERR "mergecons.pl: cannot create $fn\n"), return;
    print $f ">$_\n$seq2{$_}\n" for @nlist;
    close $f;
  }

  @flist;
}

sub mergecons {


  my $out; $out = $1, pop if @_ && $_[-1] =~ /-out=(.*)/;
  my $reftree = pop if @_ && $_[-1] !~ /\.(log|fsa)$/;
  my $ref = shift or print(STDERR "mergecons.pl: missing required arguments\n"), return{};

  my $tmpdir = ($ref =~ m#(.*)/# ? $1 : ".")."/mergetmp";
  mkdir $tmpdir or print(STDERR "mergecons.pl: creating $tmpdir failed\n"), return{};

  my $in = $_[0] || "$ref.log";
  if($in =~ /\.log$/) {
    @_ > 1 and print(STDERR "mergecons.pl: you cannot specify more than one StatAlign log file\n"), return{};
    @_ = mergecons_xlog($in,$tmpdir) or print(STDERR "mergecons.pl: error reading log file\n"), return{};
  } else {
    for my $i (@_) {
      (my $j = $i) =~ s#.*/##;
      copy $i, "$tmpdir/$j" or print(STDERR "mergecons.pl: copying $i to $tmpdir/$j failed\n"), return{};
    }
  }

  $out or ($out=$ref) =~ s/\.[^.]+$/.mergecons.fsa/;

  my $time = -time;
  $_=`mergealign -a $tmpdir -f $out 2>&1`;
  print(),return{} if $?;
  $time += time;

  unlink <$tmpdir/*> or print(STDERR "mergecons.pl: error removing files");
  rmdir $tmpdir or print(STDERR "mergecons.pl: error removing temp dir");

  require "${path}scoreall.pl";
  my $score = scoreall($ref,$out,$reftree//());

  { _infile=>$in,
    %$score,
    ztime=>$time };

}

1;
