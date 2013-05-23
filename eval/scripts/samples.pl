#!/usr/bin/env perl
#
# Scores each sample from a given StatAlign log file using scoreall.pl.
# Prints one line in the output file per sample, columns are scores.
#
# Input: [-h] refalign.ext [reftree] [input.log]
#  default input is refalign.ext.log
#  output is written into input.log.samples
#  -h prints a header in the first line
#
# Output:
#  none
#

use warnings;

$path = $0 =~ m%(.*[/\\])% ? $1 : "";
require "${path}module.pl";

sub samples_score {
  my ($ref,$test,$tree,$h,$n) = @_;
  reformat($h//(), 'scoreall', $ref, $test, $tree//(), "-_infile=$n");
}

sub samples {

  require "${path}reformat.pl";

  my $h; $h = shift if $_[0] =~ /^-h/;
  my $ref = shift or print(STDERR "samples.pl: reference alignment is missing\n"), return{};
  my $reftree; $reftree = shift if @_ && $_[0] !~ /.log$/;
  my $in = shift // "$ref.log";
  my $out = "$in.samples";

  my $tmp = "$in.sampletmp.fsa";

  open my $f, $in or print(STDERR "samples.pl: error opening input file\n"), return{};
  open my $g, ">$out" or print(STDERR "samples.pl: error creating output file\n"), return{};
  select $g;

  my $cur = -1;
  my $hh;
  while(<$f>) {
    chomp;
    if(/Sample (\d+)\tAlignment:\t(.*)/) {
      my ($n,$ali) = ($1,$2);
      if($cur != $n) {
        close($hh), samples_score($ref,$tmp,$reftree,$h,$cur), undef($h) if $cur != -1;
        open $hh, ">$tmp";
        $cur = $n;
      }
      print $hh "$ali\n";
    }
  }
  close($hh), samples_score($ref,$tmp,$reftree,$h,$cur) if $cur != -1;

  close $f;
  close $g;
  unlink $tmp if -e $tmp;

  {};
}

1;
