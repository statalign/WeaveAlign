#!/usr/bin/perl
#
# a
#
# Input: a
#
# Output:
#  a
#

use warnings;

$path = $0 =~ m%(.*[/\\])% ? $1 : "";
require "${path}module.pl";

sub ?? {

  {};
}

1;