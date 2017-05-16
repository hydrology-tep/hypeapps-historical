# build hype_assimilation  and copy to ../bin-file directory
cd HypeDataAssimilation
make -fmakefile_assimilation
cp ./hype_assimilation ../../bin/hype_assimilation
make .fmakefile_assimilation clean
rm hype_assimilation
cd ..
