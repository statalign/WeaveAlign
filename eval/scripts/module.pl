#!/usr/bin/env perl
#
# Calls the method named after the main file if it exists. Expects a
# hash reference as the return value the contents of which are then
# sorted and printed.
#
# Provides a way of creating a set of modules that can either be built
# on top of each other or run separately.
#

#use warnings;

if(!@ARGV) {
  open my $f, "$0";
  <$f>; #<$f>;
  while(<$f>) {
    last unless /^#/;
    print substr $_,1;
  }
  close $f;
  exit 1;
}
$0 =~ m%.*?([^/\\]*?)(.pl)?$%;
#print "$0\n$1\n";
if(exists &$1) {
  my $map = &$1(@ARGV);
  print "$_\t$$map{$_}\n" for sort keys %$map;
} else {
#  print STDERR "module.pl: method $1 does not exist\n";
}

sub slurp {
    local $/ = wantarray ? $/ : undef;
    open my $f, $_[0];
    my $cont = <$f>;
    close $f;
    $cont;
}

1;
