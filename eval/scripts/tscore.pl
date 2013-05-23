#!/usr/bin/env perl
#
# Creates a ML tree from the given alignment using Phylip's DNAML and
# checks its accuracy to a reference tree if provided (Newick format).
#
# Alignment is converted to Fasta first if it is in a different format
# in which case the temporary file is removed afterwards.
#
# Input: align.ext [reftree]
#  output file name is align.dnaml.tree
#
# Output:
#  _infile       input file (after conversion)
#  _outfile      output tree file
#  [scores]      treedist.pl's output fields (if reftree given)
#

use warnings;

use Cwd;

$path = $0 =~ m%(.*[/\\])% ? $1 : "";
require "${path}module.pl";

sub tscore {

  my $phout = "outfile";
  my $phtree = "outtree";
  my $nul = $^O =~ /win/i ? "nul" : "/dev/null";

  my $infile = shift or printf(STDERR "tscore.pl: input alignment not specified\n"), return{};
  my $reftree = shift;

  require "${path}convert.pl";
  my $iconv;
  ($infile,$iconv) = @{convert($infile,'fsa')}{'out','converted'}, $infile or return{};

  my ($tpath,$tname) = $infile =~ m%(.*[/\\])(.*)% ? ($1,$2) : (".",$infile);
  my $dir = getcwd;
  chdir $tpath;

  `clustalw2 -convert -output=phylip -outfile=infile -infile=$tname 2>$nul`;
  -e "infile" or print(STDERR "tscore.pl: error converting input file to Phylip format\n"), return{};

  unlink $phout if -e $phout;
  unlink $phtree if -e $phtree;
  open my $f, "|dnaml >$nul";
  print $f <<EOL;
Y

EOL
  close $f;

  (my $outtree = $tname) =~ s/.fsa$/.dnaml.tree/;
  rename $phtree, $outtree or printf(STDERR "treedist.pl: DNAML stopped with an error\n"), return{};
  unlink $phout if -e $phout;
  unlink "infile";

  chdir $dir;
  unlink $infile if $iconv;
  ($outtree = $infile) =~ s/.fsa$/.dnaml.tree/;

  my $scores = {};
  if($reftree) {
    require "${path}treedist.pl";
    $scores = treedist($reftree,$outtree);
  }

  { _infile  => $infile,
    _outfile => $outtree,
    %$scores };

}

1;
