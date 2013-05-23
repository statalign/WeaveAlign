#!/usr/bin/perl
#
# Scores test alignment by comparing it against a reference alignment. Some
# checks are performed on the alignments, currently the individual
# characters of the sequences are also matched for the two alignments.
#
# Converts both alignments into FSA format first (temporarily) if
# they're in a different format, see convert.pl. If conversion took
# place temporary files are removed afterwards.
#
# You can specify a tree (in Newick format) to assign a weight of 1/d_ij to
# each sequence pair when evaluating pairwise scores where d_ij is the
# distance between sequence i and j in the tree.
#
# Input: refalign testalign [scorefile] [treefile]
#
# Output:
#  _infile      reference file
#  baliscore    ratio of correct residue-residue pairs
#  fsascore     ratio of all correct pairs (residue pairs double-weighted)
#  nhomscore    ratio of correct residue-gap pairs
#  scolscore    ratio of correct full columns where gaps are not distinguished
#  smpdscore    ratio of correct full columns where gaps are distinguished
#
use lib "/data/genome/lib/perl5";

use warnings;
use strict;

use Hash::MultiKey;

my $path = $0 =~ m%(.*[/\\])% ? $1 : "";
require "${path}module.pl";

sub score_readfsa {
  open my $in, $_[0] or return ();

  my $name;
  my %seq;

  while(<$in>) {
    chomp;
    s/[\s\r]//g;
    /^> *(.*?) *$/ and $name=uc($1), next;
    next if length == 0;
    $seq{$name} .= uc;
  }

  close $in;

  %seq;
}

sub score_readtree {
  open my $in, $_[0] or return;
  my $tree = join '', <$in>;
  close $in;
  $tree =~ s/\s//g;
  $tree =~ s/([^:(),;]+)/"'".uc($1)."'"/ge;
  $tree =~ y/:()/,[]/;
  $tree = eval $tree;
  if(@$tree > 4) {
    $$tree[1] /= 2;
    @$tree[2..5] = ([@$tree[2..5]], $$tree[1]);
  }
  $tree;
}

sub tree2dist {
  my $tree = shift;
  my $nameh = shift;
  my $dist = [];

  t2d_rec($tree, $nameh, $dist);
  return $dist;

  sub t2d_rec {
    my $subt = shift;
    my $nh = shift;
    my $d = shift;

    if(ref $subt) {
      my $l = t2d_rec($$subt[0],$nh,$d);
      $_ += $$subt[1] for values %$l;
      my $r = t2d_rec($$subt[2],$nh,$d);
      $_ += $$subt[3] for values %$r;
      for my $i (keys %$l) {
        for my $j (keys %$r) {
          $d->[$i][$j] = $d->[$j][$i] = $l->{$i}+$r->{$j};
        }
      }
      return {%$l, %$r};
    } else {
      defined(my $id = $nh->{$subt}) or print(STDERR "score.pl: tree incompatible with sequences\n"), exit 1;
      return {$id => 0};
    }
  }
}

sub score {

  my $checkchars = 1;   # true if individual characters must be checked
  my $output_correct = 0; # true if we want to output a vector of which columns are correct with scores alongside

  my $ref = shift or print(STDERR "score.pl: ref file not specified\n"), return{};
  my $test = shift or print(STDERR "score.pl: test file not specified\n"), return{};
  my $score_file = "";
  #my $score_file = shift; 
  #if ($score_file) {
  #    $output_correct = 1;
  #}

  my $treef = shift;

  my $correct_file = $test.".correct";

  require "${path}convert.pl";
  my ($rconv,$tconv);
  ($ref,$rconv) = @{convert($ref,'fsa')}{'out','converted'}, $ref or return{};
  ($test,$tconv) = @{convert($test, 'fsa')}{'out','converted'}, $test or return{};

  my %rali = score_readfsa($ref) or print(STDERR "score.pl: error reading ref file\n"), return{};
  my %tali = score_readfsa($test) or print(STDERR "score.pl: error reading test file\n"), return{};
  !$treef or my $tree = score_readtree($treef) or print(STDERR "score.pl: error reading tree file\n"), return{};

  unlink $ref if $rconv;
  unlink $test if $tconv;

  # sort sequences
  my @rnames = sort keys(%rali);
  my @tnames = sort keys(%tali);

  # perform checks
  @rnames == @tnames or print(STDERR "score.pl: number of sequences in ref and test do not match\n"), return{};
  my ($rlen, $tlen);
  for my $i (0..$#rnames) {
    $rnames[$i] lt $tnames[$i] and print(STDERR "score.pl: seq $rnames[$i] not found in test\n"), return{};
    $rnames[$i] gt $tnames[$i] and print(STDERR "score.pl: seq $tnames[$i] not found in ref\n"), return{};
    my $name = $rnames[$i];
    my $rseq = $rali{$name};
    my $tseq = $tali{$name};
    if($rlen) {
      length($rseq) == $rlen or print(STDERR "score.pl: ref file not aligned, seq $name is of different length\n"), return{};
    } else {
      $rlen = length($rseq);
    }
    if($tlen) {
      length($tseq) == $tlen or print(STDERR "score.pl: test file not aligned, seq $name is of different length\n"), return{};
    } else {
      $tlen = length($tseq);
    }
    $rseq =~ y/-//d;
    $tseq =~ y/-//d;
    $checkchars and $rseq ne $tseq and print(STDERR "score.pl: seqs $name in ref and test do not match\n"), return{};
    length($rseq) == length($tseq) or print(STDERR "score.pl: seqs $name in ref and test do not match\n"), return{};
  }

  # preprocess ref alignment
  my ($rali, $rcol);
  tie my %hmpd, 'Hash::MultiKey';
  tie my %hcol, 'Hash::MultiKey';
  my $mpd = [(0) x @rnames];
  my $col = [(0) x @rnames];
  for my $c (0..($rlen-1)) {
    my $allg = 1;
    for my $i (0..$#rnames) {
      $_ = substr $rali{$rnames[$i]}, $c, 1;
      my $j = $mpd->[$i];
      $j += $j & 1;
      $j++, $allg = 0 if $_ ne '-';
      $mpd->[$i] = $j;
      my $ch = ($j & 1) ? $j >> 1 : -1;
      $col->[$i] = $ch;
      $rali->[$i][$c] = $ch;
      $rcol->[$i][$ch] = $c if $ch >= 0;
    }
    $hmpd{$mpd} = 1, $hcol{$col} = 1 if !$allg;
  }
  keys(%hmpd) == keys(%hcol) or die "internal error";

  # preprocess test alignment + calc column scores
  my ($tali, $tcol);
  my ($mpdok, $colok)=(0,0);
  my $colt = 0;
  $mpd = [(0) x @rnames];
  $col = [(0) x @rnames];

  my $correct = [(0) x $tlen]; 
  my $mpd_correct = [(0) x $tlen]; 

  # Binary vector yielding presence/absence of column
  # in test alignment, using MPD (S minus) score.
  for my $c (0..($tlen-1)) {
    my $allg = 1;
    for my $i (0..$#tnames) {
      $_ = substr $tali{$tnames[$i]}, $c, 1;
      my $j = $mpd->[$i];
      $j += $j & 1;
      $j++, $allg = 0 if $_ ne '-';
      $mpd->[$i] = $j;
      my $ch = ($j & 1) ? $j >> 1 : -1;
      $col->[$i] = $ch;
      $tali->[$i][$c] = $ch;
      $tcol->[$i][$ch] = $c if $ch >= 0;
    }
    $colt++ if !$allg;
    #$mpdok++ if $hmpd{$mpd};
    if ($hmpd{$mpd}) {
	$mpdok++;
	$mpd_correct->[$c] = 1;
    }
    else {
	$mpd_correct->[$c] = 0;
    }
    if ($hcol{$col}) {
	$colok++;
	$correct->[$c] = 1;
    }
    else {
	$correct->[$c] = 0;
    }
  }

  

  # get pairwise distances
  my $d;
  if($tree) {
    my $nh = {};
    $nh->{$rnames[$_]} = $_ for 0..$#rnames;
    $d = tree2dist($tree,$nh);
    @$d == @rnames or print(STDERR "score.pl: tree incompatible with sequences\n"), return{};
  } else {
    $d->[$_] = [(1)x@rnames] for 0..$#rnames;
  }

  # calc pairwise scores
  my ($homok,$nhomok) = (0,0);
  my ($homn,$nhomn) = (0,0);
  my ($homnt,$nhomnt) = (0,0);
  for my $i (0..$#rnames) {
    for my $j (0..$#{$rcol->[$i]}) {
      my $rc = $rcol->[$i][$j];
      my $tc = $tcol->[$i][$j];
      for my $k (0..$#rnames) {
        next if $i == $k;
        my $w = 1/$d->[$i][$k];
        my $hom = $rali->[$k][$rc] >= 0;
        $hom ? $homn += $w : ($nhomn += $w);
        $hom ? $homok += $w : ($nhomok += $w) if $rali->[$k][$rc] == $tali->[$k][$tc];
        $tali->[$k][$tc] >= 0 ? $homnt += $w : ($nhomnt += $w);
      }
    }
  }

  if ($output_correct) {
      open(SCORES, "< $score_file") or die "Cannot open test file: $!";
      print $score_file."\n";
      my @scores;
      #my $scores_started = 0;
      while(<SCORES>) {
	  #if (/#scores/) {
	  #    $scores_started = 1;
	  #    next;
	  #}
	  #next unless $scores_started;
	  if (/(\d\.\d+)/) {
	      push(@scores,$1);
	  }
      }
      close SCORES;
      open(COR,"> $correct_file");
      for (my $i=0; $i<$tlen; $i++) {
	  print COR $mpd_correct->[$i]."\t".$correct->[$i]."\t".$scores[$i]."\n";
      }
      close COR;
  }	  
  
  {_infile => $ref,
   alignlen => $colt,

   baliscore => sprintf("%.3f", $homok/$homn),
   fsascore => sprintf("%.3f", ($homok+$nhomok)/($homn+$nhomn)),
   nhomscore => sprintf("%.3f", $nhomok/$nhomn),
   scolscore => sprintf("%.3f", $colok/keys(%hcol)),
   smpdscore => sprintf("%.3f", $mpdok/keys(%hmpd)),

   baliscorem => sprintf("%.3f", $homok/$homnt),
   fsascorem => sprintf("%.3f", ($homok+$nhomok)/($homnt+$nhomnt)),
   nhomscorem => sprintf("%.3f", $nhomok/$nhomnt),
   scolscorem => sprintf("%.3f", $colok/$colt),
   smpdscorem => sprintf("%.3f", $mpdok/$colt)};
}

1;
