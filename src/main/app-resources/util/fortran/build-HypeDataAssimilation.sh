# build hype_assimilation  and copy to ../bin-file directory
unzip HypeDataAssimilation.zip
make -fmakefile_assimilation
cp ./hype_assimilation ../bin/hype_assimilation
make .fmakefile_assimilation clean
rm *.f90
rm steplength_dependent.txt
rm MAK_HYPE.BAT
rm hype_assimilation
rm makefile_assimilation
rm makefile
ls