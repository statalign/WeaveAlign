#!/bin/bash

thisDir=$(dirname $0) || false

pushd ${thisDir} > /dev/null
  gradle clean build jar
popd > /dev/null

