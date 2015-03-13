#!/bin/bash

# Before running, read the README, and run the "build-all.sh" script.

# This script will compute the rank score for the minimum-risk
# summary alignment, relative to the set of alignments contained
# within a log file, with respect to a reference alignment, 
# for the OXBench datasets described in 

# Herman JL, Novák Á, Lyngsø R, Szabó A, Miklós I and Hein J (2014)
#  "Efficient representation of uncertainty in multiple sequence alignments 
# using directed acyclic graphs."

nSamples=2000
# Running with nSamples=2000 will reproduce the results in the paper,
# but may take a couple of minutes to execute.

thisDir=$( cd "$( dirname "${BASH_SOURCE[0]}" )" && pwd )

pushd $thisDir > /dev/null

  g=0.5
  burnin=300
  refAliDir=$thisDir/testdata/oxbench/reference_alignments
  logFileDir=$thisDir/testdata/oxbench/sample_alignments
  outFile=$thisDir/testdata/oxbench/wvaAnalysis.dat
  
  # evaluate the results
  for m in 15 33 60 122; do 
      logFile=$logFileDir/oxbench_12_size_$m.log
      if [ ! -f "$logFile.fsa" ]; then	  
	  cd $logFileDir
	  unzip -n ali_samples.zip > /dev/null
	  cd $thisDir
	  printf "Computing minimum-risk summary alignment..."
	  java -jar $thisDir/WeaveAlign.jar -g=$g -f=$burnin -n=$nSamples $logFileDir/oxbench_12_size_$m.log > /dev/null
	  printf "done.\n"
      fi
      java -cp $thisDir/WeaveAlign.jar wvalign.eval.AlignSetEval $refAliDir/oxbench_12_size_$m.fasta $logFile.fsa $logFile $outFile $burnin $(( nSamples + burnin - 1 ))
      echo ""
      echo ""
  done
    
popd > /dev/null
