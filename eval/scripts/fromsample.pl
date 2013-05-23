#!/usr/bin/perl

$from = shift;
$to = shift;
$in = shift;

open IN, $in;

while(<IN>) {
 last if /^Sample $from/;
}

s/^Sample ([0-9]+)/"Sample ".($1-$from)/e;
print;

while(<IN>) {
 last if /^Sample $to/;
 s/^Sample ([0-9]+)/"Sample ".($1-$from)/e;
 print;
}
