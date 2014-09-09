#!/bin/bash

# At first read the README, run the "build-all.sh" script!

thisDir=$(dirname $0) || false

pushd ${thisDir} > /dev/null

  inputFsaDir=${thisDir}/testdata/exampleInputDir
  logFileForWvaInput=${thisDir}/testdata/exampleAlignments.log
  resultFile=${thisDir}/testdata/wva-result-example500.fsa
  
  # clean up (maybe prev runs)
  rm -rf ${inputFsaDir}
  # recreate directory of input alignments
  mkdir -p ${inputFsaDir}
  
  # unzip example input (500 randomized alignments)
  echo "Unzipping alignment samples..."
  unzip testdata/sampleAligns.zip -d ${inputFsaDir} > /dev/null
  echo "...done"
  echo "Concatenating alignments into log file..."
  # (WeaveAlign will run much faster this way than requiring it to read in many small fsa files)
  ./create-wva-input.sh ${inputFsaDir} 500 ${logFileForWvaInput}
  echo "...done"
  
  # start WeaveAlign
  echo "Running WeaveAlign..."
  java -jar runnable-jars/WeaveAlign-1.2-all.jar -out ${resultFile} -g=0.5 ${logFileForWvaInput}

  # evaluate the results
  echo ""
  echo "Summary of results:"
  echo "-------------------"
  java -jar runnable-jars/WeaveAlignEval-1.2-all.jar ${thisDir}/testdata/testrun500.properties
    
popd > /dev/null
