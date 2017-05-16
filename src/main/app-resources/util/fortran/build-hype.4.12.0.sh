# build hype v4.12.0 and copy to ../bin-file directory
cd HYPE-4.12.0
make
cp ./hype ../../bin/hype-4.12.0
make clean
rm hype
cd ..
