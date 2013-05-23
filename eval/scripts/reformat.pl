#!/usr/bin/perl
#
# Runs a module with an optional list of arguments and prints a tab-separated
# list of output values to STDOUT (sorted by keys)
#
# Input: [-h[eader]] module [params...]
#  -header   prints a line of output keys as header
#
# Output:
#  none
#

use warnings;

$path = $0 =~ m%(.*[/\\])% ? $1 : "";
require "${path}module.pl";

sub reformat {

  $_ = shift;
  my $head = /^-h/ and ($_=shift);
  my $module = $_ or die "reformat.pl: missing module name";

  require "${path}${module}.pl";

  my $result = &$module(@_);

  local ($\,$,)=("\n","\t");
  print sort keys %$result if $head;
  print map $$result{$_}, sort keys %$result;

  +{};
}

1;