#!/bin/bash

thisDir=$(dirname $0) || false

pushd ${thisDir} > /dev/null
  mkdir -p runnable-jars
  ./build-wva.sh
  cp build/libs/*-all.jar runnable-jars/
  ./build-eval.sh
  cp build/libs/*-all.jar runnable-jars/
  cp runnable-jars/WeaveAlign-*-all.jar WeaveAlign.jar
popd > /dev/null

