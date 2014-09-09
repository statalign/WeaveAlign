#!/bin/bash

thisDir=$(dirname $0) || false

pushd ${thisDir} > /dev/null
  gradle -b buildEval.gradle clean build shadowJar
popd > /dev/null

