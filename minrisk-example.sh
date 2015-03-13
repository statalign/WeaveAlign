#!/bin/bash

# Before running, read the README, and run the "build-all.sh" script

# This script will compute the minimum-risk summary alignment
# for the alignments contained in a set of log files, for the OXBench
# datasets described in 

# Herman JL, Novák Á, Lyngsø R, Szabó A, Miklós I and Hein J (2014)
#  "Efficient representation of uncertainty in multiple sequence alignments 
# using directed acyclic graphs."

thisDir=$( cd "$( dirname "${BASH_SOURCE[0]}" )" && pwd )

echo $thisDir

pushd $thisDir > /dev/null

  g=0.5
  burnin=300
  nSamples=2000
  refAliDir=$thisDir/testdata/oxbench/reference_alignments
  logFileDir=$thisDir/testdata/oxbench/sample_alignments
  outFile=$thisDir/testdata/oxbench/wvaAnalysis.dat

  cd $logFileDir
  unzip -n ali_samples.zip > /dev/null
  cd $thisDir

  # start WeaveAlign
  echo "Running WeaveAlign..."
  for m in 15 33 60 122; do 
      java -jar $thisDir/WeaveAlign.jar -g=$g -f=$burnin -n=$nSamples $logFileDir/oxbench_12_size_$m.log
  done
    
popd > /dev/null
