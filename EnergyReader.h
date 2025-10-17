#ifndef EnergyReader_h 
#define EnergyReader_h 1
#include "DetGeo.h"
#include <TRandom3.h>
#include <vector>
#include <fstream>
#include <iostream>
#include <array>
#include <algorithm>
#include <string>

using namespace std;

class EnergyReader {
public:
  EnergyReader();
  void ReadEnergy(int evID, vector<int>* vADCData, vector<int>* vTrueHitCubeN, vector<int>* vTrueHitZ, vector<double>* vTrueHitEnergy);

private:
};

#endif
