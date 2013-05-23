#!/usr/bin/perl
#
# Runs MPD on a StatAlign log file or a set of alignments in Fasta format and
# assesses resulting MPD alignment quality using a reference alignment.
#
# Input: ref.ext [reftree] input_1.[log|fsa],... [-out=output.fsa] [params...]
#  ref.ext is the reference alignment
#  default output is input_1.ext.fsa
#  optional params... are transferred to MPD directly
#  MPD process is logged into output.log
#
# Output:
#  _infile       reference alignment
#  [scores]      scores from score(all).pl
#  ztime         elapsed time (seconds)
#

use warnings;

$path = $0 =~ m%(.*[/\\])% ? $1 : "";
require "${path}module.pl";

sub mpd {

  my $ref = shift or print(STDERR "mpd.pl: reference file not specified\n"), return{};
  my $reftree;
  $reftree = shift if @_ && $_[0] !~ /\.(log|fsa)$/;

  my $in0 = $_[0] or print(STDERR "mpd.pl: input file(s) not specified\n"), return{};
  my ($out) = grep /^-out/, @_;
  $out ? $out =~ s/^-out=// : ($out = "$in0.fsa");
  (my $log = $out) =~ s/\.[^.]*$/.log/;

  # canonicalisation and -out fix
  my @temps;
  for(@_) {
    s/^-out=/-out / and next;
    next if !/.fsa$/;
    `${path}canonicalise -i $_ 2>&1`;
    $_ .= ".canon";
    rename $_, "$_.fsa";
    $_ .= ".fsa";
    push @temps, $_;
  }

  my $pars = join ' ', @_;
  my $time = -time;
  `java -Xmx512m -jar ${path}mpd.jar -outgi $pars >$log`, $_=slurp($log);
  print(),return{} if $?;
  $time += time;
#  unlink $log;

  unlink @temps;

  require "${path}convert.pl";
  $outfsa = ${convert($out,'fsa')}{'out'} or return{};

  require "${path}scoreall.pl";
  my $score = scoreall($ref,$outfsa,$reftree//());

#  unlink $outfsa;

  { _infile=>$ref,
    %$score,
    ztime=>$time };
}

1;
