#!/bin/bash

perl generate_data.pl
mafft t091/data.fsa >t091/data.mafft.fsa
perl ../scripts/score.pl t091/data.fsa t091/data.mafft.fsa