#!/bin/bash

thisDir=$(dirname $0) || false

pushd ${thisDir}
  gradle distZip
popd

