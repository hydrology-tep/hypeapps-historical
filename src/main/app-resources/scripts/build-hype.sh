#!/bin/bash

# compile HYPE model code and copy to /application/util/bin
mkdir ../util/bin
pushd ../util/fortran/HYPE-5.x.0
make -fmakefile_assimilation
cp ./hype_assimilation ../../bin/hype_assimilation-5.x.0.exe
make -fmakefile_assimilation clean
rm hype_assimilation
popd