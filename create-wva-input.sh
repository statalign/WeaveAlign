#!/bin/bash

thisDir=$(dirname $0) || false

inputDirOfAligns=$1
numOfSamples=$2
logFileForWva=$3


pushd ${thisDir} > /dev/null

echo -n "" > ${logFileForWva}

for i in $(seq 1 ${numOfSamples})
do 
  fsaFile="${inputDirOfAligns}/sampleAligns/m15-BL62-rndmx-${i}.fsa"
  cat ${fsaFile} | awk -v num=${i} '{print "Sample " (num-1) "\tAlignment:\t" $0}' >> ${logFileForWva}
done

popd > /dev/null
