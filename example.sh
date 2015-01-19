#!/bin/bash

# Before running, read the README, and run the "build-all.sh" script

thisDir=$(dirname $0) || false

pushd $thisDir > /dev/null

  g=0.5
  burnin=300
  maxSample=2000 # Must be larger than $burnin
  refAliDir=$thisDir/testdata/oxbench/reference_alignments
  logFileDir=$thisDir/testdata/oxbench/sample_alignments
  outFile=$thisDir/testdata/oxbench/wvaAnalysis.dat

  cd $logFileDir
  unzip ali_samples.zip > /dev/null
  cd $thisDir

  # start WeaveAlign
  echo "Running WeaveAlign..."
  for m in 15 33 60 122; do 
      java -jar WeaveAlign.jar -g=$g -f=$burnin -n=$maxSample $logFileDir/oxbench_12_size_$m.log
  done

  # evaluate the results
  echo ""
  echo "Summary of results:"
  echo "-------------------"
  for m in 15 33 60 122; do 
      logFile=$logFileDir/oxbench_12_size_$m.log
      java -cp WeaveAlign.jar wvalign.eval.AlignSetEval $refAliDir/oxbench_12_size_$m.fasta $logFile.fsa $logFile $outFile $burnin $maxSample
      echo ""
      echo ""
  done
    
popd > /dev/null
