use warnings;

@ARGV >= 3 && @ARGV % 2 == 1 or print('Use: subst.pl infile pattern substfile [pattern substfile]...'), exit;
-e $ARGV[$_*2] or die "file ".$ARGV[$_*2]." doesn't exist" for (1..@ARGV/2);

open F, shift;
my %patts = @ARGV;
lab: while(defined(my $l = <F>)) {
  $l =~ /^(.*)\Q$_\E(.*)$/ and print($1), printfile($patts{$_}), print($2), next lab for keys %patts;
  print $l;
}
close F;

sub printfile {
  open G, $_[0];
  print while(<G>);
  close G;
}