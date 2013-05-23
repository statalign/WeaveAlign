#!/bin/bash
#
# Usage: ./run.sh [module] subdir/*.msf
# 
# module defaults to runall
#

if [ -z "$1" ]; then echo "Usage: ./run.sh [module] [subdir/*/*.rsf]"; exit 1; fi

out="results/out.xls"

export PATH=$PATH:scripts

if [ -e "scripts/$1.pl" ]; then cmd=$1; shift; else cmd=runall; fi

list=$@
if [ -z "$list" ]; then list="*/*.msf"; fi

if [ -e $out ]; then rm $out; fi
head="-h"
for i in $list; do
 echo $i
 /usr/bin/perl scripts/reformat.pl $head $cmd $i >>$out
 head=
done
