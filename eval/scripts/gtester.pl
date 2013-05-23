#!/usr/bin/perl
#
# Tests MPD with different g values.
#
# Input: file.ext
#
# Output:
#  _input    Input file
#  mpd.pl's each non-underscored output field suffixed by g value
#

use warnings;

$path = $0 =~ m%(.*[/\\])% ? $1 : "";
require "${path}module.pl";

sub gtester {

  my $in = shift;
  (my $log = $in) =~ s/[^.]+$/fasta.log/;

#  my @test_vals = (0,0.5,1,1.5);
  my @test_vals = map $_/20.0, 0..10;

  my %result = (_input=>$in);
  require "${path}mpd.pl";


  (my $mpd = $log) =~ s/log$/mpd/;
  require "${path}convert.pl";
  $outfsa = ${convert($mpd,'fsa')}{'out'} or return{};

  require "${path}score.pl";
  my $oret = score($in,$outfsa);
  $result{"${_}_ori"}=$$oret{$_} for grep !/^_/, keys %$oret;

  unlink $outfsa;


  for (@test_vals) {
    my $g = sprintf "%.02f", $_;
    (my $out = $log) =~ s/log$/g$g.mpd/;
    my $ret = mpd($in,$log,$out,"-g=$g");
    $result{"${_}_g$g"}=$$ret{$_} for grep !/^_/, keys %$ret;
  }

  \%result;
}

1;
