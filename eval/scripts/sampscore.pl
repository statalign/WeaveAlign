#!/usr/bin/env perl
#
# Scores each sample from a given StatAlign log file based on
# 'internal' scores of the MPD procedure, i.e. the sum of column
# marginals for the alignment divided by the reference length
# (developer) or the alignment sample length (modeler).
#
# Prints one line in the output file per sample, columns are tab-
# separated scores (sample no, developer score, modeler score)
#
# Input: [-h] refalign.ext [input.log]
#  default input is refalign.ext.log
#  output is written into input.log.sampscores
#  -h prints a header in the first line
#
# Output:
#  none
#

use warnings;

$path = $0 =~ m%(.*[/\\])% ? $1 : "";
require "${path}module.pl";

sub sampscore_score {
  my ($ref,$test,$tree,$h,$n) = @_;
  reformat($h//(), 'scoreall', $ref, $test, $tree//(), "-_infile=$n");
}

sub sampscore_readfsa {
  open my $in, $_[0] or return ();

  my $name;
  my %seq;

  while(<$in>) {
    chomp;
    s/[\s\r]//g;
    /^> *(.*?) *$/ and $name=uc($1), next;
    s/\s+//g;
    next if length == 0;
    $seq{$name} .= uc;
  }

  close $in;

  %seq;
}


sub sampscore {

  require "${path}reformat.pl";

  my $h; $h = shift if $_[0] =~ /^-h/;
  my $ref = shift or print(STDERR "sampscore.pl: reference alignment is missing\n"), return{};
  my $in = shift // "$ref.log";
  my $out = "$in.sampscores";
  my $tmp = "$out.tmp";

  $_=`java -Xmx512m -jar ${path}mpd_score.jar $in -out $tmp 2>&1`;
  print(),return{} if $?;

  my %ali = sampscore_readfsa($ref) or print(STDERR "sampscore.pl: failed to read ref\n"), return{};
  my $alilen = length $ali{(keys %ali)[0]};

  open my $f, $tmp or print(STDERR "sampscore.pl: mpd_score.jar's output ($tmp) not found\n"), return{};
  open my $g, ">$out" or print(STDERR "sampscore.pl: failed to create output file\n"), return{};
  while(<$f>) {
    my @A = split;
    print $g $A[0],"\t",$A[1]/$alilen,"\t",$A[1]/$A[2],"\n";
  }
  close $f;
  close $g;

  unlink $tmp if -e $tmp;

  {};
}

1;
