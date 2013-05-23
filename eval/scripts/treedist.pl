#!/usr/bin/env perl
#
# Compares two trees using Phylip's treedist with two methods.
#
# Input: treefile1 treefile2
#
# Output:
#  brsdist       Branch Score Distance of Kuhner & Felsenstein (1994)
#  symdiff       Symmetric Difference of Robinson & Foulds (1981)
#

use warnings;
use Cwd;

$path = $0 =~ m%(.*[/\\])% ? $1 : "";
require "${path}module.pl";

sub treedist {

  my $phout = "outfile";
  my $nul = $^O =~ /win/i ? "nul" : "/dev/null";

  my $t1 = shift; $t1 && -e $t1 or die "missing treefile1";
  my $t2 = shift; $t2 && -e $t2 or die "missing treefile2";

  my ($t1path,$t1name) = $t1 =~ m%(.*[/\\])(.*)% ? ($1,$2) : (".",$t1);
  my ($t2path,$t2name) = $t2 =~ m%(.*[/\\])(.*)% ? ($1,$2) : (".",$t2);

  my $dir = getcwd, chdir $t1path, $t1 = $t1name, $t2 = $t2name if $t1path eq $t2path;

  unlink $phout if -e $phout;
  open my $f, "|treedist >$nul";
  print $f <<EOL;
$t1
2
C
S
Y
$t2

EOL
  close $f;

  my $brsdist = sprintf "%.3f", (split " ", slurp($phout))[1];

  unlink $phout;
  open $f, "|treedist >$nul";
  print $f <<EOL;
$t1
D
2
C
S
Y
$t2

EOL
  close $f;

  my $symdiff = sprintf "%.3f", (split " ", slurp($phout))[1];
  unlink $phout;

  chdir $dir if $dir;

  { brsdist  => $brsdist,
    symdiff => $symdiff };

}

1;
