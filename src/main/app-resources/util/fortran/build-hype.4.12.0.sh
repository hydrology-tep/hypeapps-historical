# build hype v4.12.0 and copy to ../bin-file directory
unzip src-open-4.12.0.zip
make
cp ./hype ../bin/hype-4.12.0
make clean
rm *.f90
rm steplength_dependent.txt
rm MAK_HYPE.BAT
rm hype
rm makefile
ls