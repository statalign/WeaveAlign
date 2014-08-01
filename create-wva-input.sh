#!/bin/bash

thisDir=$(dirname $0) || false

inputDirOfAligns=$1
numOfSamples=$2
logFileForWva=$3

fakeSummaryLine1="Acceptances	Alignment	0.5555	Edge	0.9000	Topology	0.1000	Indel	0.5000	Substitution	0.0"
fakeSummaryLine2="Report	LogLikelihood	-1000.5555	R	0.8888	Lamda	0.0111	Mu	0.0555	"

pushd ${thisDir}

echo -n "" > ${logFileForWva}

for i in $(seq 1 ${numOfSamples})
do 
  fsaFile="${inputDirOfAligns}/sampleAligns/m15-BL62-rndmx-${i}.fsa"
  cat ${fsaFile} | awk -v num=${i} '{print "Sample " (num-1) "\tAlignment:\t" $0}' >> ${logFileForWva}
  echo ${fakeSummaryLine1} >> ${logFileForWva}
  echo ${fakeSummaryLine2} >> ${logFileForWva}
done

popd
