#!/bin/bash

# At first read the README, run the "build.sh" script!

thisDir=$(dirname $0) || false

pushd ${thisDir}

  inputFsaDir=${thisDir}/testdata/exampleInputDir
  logFileForWvaInput=${thisDir}/testdata/exampleAlignments.log
  resultFile=${thisDir}/testdata/wva-result-example500.fsa
  
  # clean up (maybe prev runs)
  rm -rf ${inputFsaDir}
  # recreate directory of input alignments
  mkdir -p ${inputFsaDir}
  
  # unzip example input (500 randomized alignments)
  unzip testdata/sampleAligns.zip -d ${inputFsaDir}
  ./create-wva-input.sh ${inputFsaDir} 500 ${logFileForWvaInput}
  
  # start WeaveAlign

  java -cp $(ls lib/*.jar | awk '{printf $0":"}'):build/libs/WeaveAlign-1.2.jar wvalign.WeaveMain -out ${resultFile} -g=0.5 ${logFileForWvaInput}

  # evaluate the results
  
  java -cp $(ls lib/*.jar | awk '{printf $0":"}'):build/libs/WeaveAlign-1.2.jar wvalign.eval.AlignEvaluator ${thisDir}/testdata/testrun500.properties
    
popd
