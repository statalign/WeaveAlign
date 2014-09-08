#!/bin/bash

thisDir=$(dirname $0) || false

pushd ${thisDir} > /dev/null
  gradle distZip
popd > /dev/null

