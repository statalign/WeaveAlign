#!/bin/bash

thisDir=$(dirname $0) || false

pushd ${thisDir} > /dev/null
  #./build-all.sh
  gradle distZip
  echo "Created zip archive: build/distributions/WeaveAlign-1.2.zip"
popd > /dev/null

