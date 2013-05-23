#!/usr/bin/perl

use warnings;
use strict;

use feature 'say';
use Hash::MultiKey;

my $x;

sub test {
  say "$_:$_[$_]" for 0..$#_;
}

test(3, $x//(), 4);
exit;

tie my %h, 'Hash::MultiKey';

my $l = [1,2,3,6];
my $ll = [1,2,3,6];

$h{$l} = 1;

$l->[2] = 7;
#say for @$l;

print "yes" if $h{$ll};