all:
	g++ -o test.exe  test.cpp LyReco.cpp FullReco.cpp EnergyReader.cpp Illustration.cpp -g `root-config --cflags --libs --glibs` -lMinuit -I./ 

