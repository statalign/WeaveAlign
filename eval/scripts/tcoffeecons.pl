#!/usr/bin/env perl
#
# Runs T-Coffee in consensus alignment mode on a set of alignments. Files
# must all be Fasta files or a single StatAlign log file with samples.
#
# Input: [-multicore] refalign.ext input_1.fsa ... [reftree] [-out=out] OR
#        [-multicore] refalign.ext [input.log] [reftree] [-out=out]
#  default output is refalign.tcoffeecons.fsa but can be overridden with -out
#  input.log is assumed to be refalign.ext.log when omitted
#  multicore mode is activated when -multicore is used
#
# Output:
#  _infile       input file (after conversion)
#  [scores]      score.pl's each output field
#  ztime         elapsed time (seconds)
#
##unused
##  treeout       Tree file created or empty
##  outfsa        Output Fasta file name
#

use strict;
use warnings;

my $path = $0 =~ m%(.*[/\\])% ? $1 : "";
require "${path}module.pl";

sub tcoffeecons_xlog {
  my $tempfiles = "%s.tcctemp.%d.fsa";

  my $in = shift;
  open my $f, $in or return;
  my $cur = -1;
  my @flist = ();
  my $g;
  while(<$f>) {
    if(/Sample (\d+)\tAlignment:\t(.*)/s) {
      if($cur != $1) {
        close $g if $cur != -1;
        open $g, ">".($_=sprintf $tempfiles, $in, $1);
        push @flist, $_;
        $cur = $1;
      }
      print $g $2;
    }
  }
  close $g if $cur != -1;
  close $f;
  @flist;
}

sub tcoffeecons {

  my $singlecore='-multi_core=no';

  $singlecore = '', shift if $_[0] eq '-multicore';
  my $out; $out = $1, pop if @_ && $_[-1] =~ /-out=(.*)/;
  my $reftree = pop if @_ && $_[-1] !~ /\.(log|fsa)$/;
  my $ref = shift or print(STDERR "tcoffeecons.pl: missing required arguments\n"), return{};

  my $in = $_[0] || "$ref.log";
  if($in =~ /\.log$/) {
    @_ > 1 and print(STDERR "tcoffeecons.pl: you cannot specify more than one StatAlign log file\n"), return{};
    @_ = tcoffeecons_xlog($in) or print(STDERR "tcoffeecons.pl: error reading log file\n"), return{};
  }

  $out or ($out=$ref) =~ s/\.[^.]+$/.tcoffeecons.fsa/;
  (my $temptree=$ref) =~ s/\.[^.]+$/.tcoffee.tmp.tree/;

  my $l = join ' ', @_;
  my $time = -time;
  $_=`t_coffee -aln $l -output=fasta -newtree=$temptree $singlecore -outfile=$out 2>&1`;
  print(),return{} if $?;
  $time += time;

  unlink @_ if $in =~ /\.log$/;
  unlink $temptree;

  require "${path}scoreall.pl";
  my $score = scoreall($ref,$out,$reftree//());

  { _infile=>$in,
    %$score,
    ztime=>$time };

}

1;
