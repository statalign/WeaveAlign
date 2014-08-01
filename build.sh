#!/bin/bash

thisDir=$(dirname $0) || false

pushd ${thisDir}
  gradle clean build jar copytolib
popd

