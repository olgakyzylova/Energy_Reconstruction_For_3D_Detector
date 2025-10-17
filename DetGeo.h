#ifndef DetGeo_h
#define DetGeo_h

const int N_PMTS = 16;
const int N_CUBES = N_PMTS*N_PMTS;
const int N_CUBE_LAYERS = 41;
const int N_PMT_LAYERS = (N_CUBE_LAYERS+1)/2;
const int N_SIDES = 4;
const int N_TOTAL_PMTS = N_PMTS * N_SIDES/2;
const int N_ITER = 7;
const double PE_TO_ADC = 6.05;

struct CubeE {
public:
  int x;
  int y;
  int z;
  int cubeN;
  double e;
  CubeE(int x, int y, int z, int cubeN, double e) {
    this -> x = x;
    this -> y = y;
    this -> z = z;
    this -> cubeN = cubeN;
    this -> e = e;
  }
};

#endif
