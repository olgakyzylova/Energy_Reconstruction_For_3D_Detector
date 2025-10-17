#ifndef FullReco_h 
#define FullReco_h 1
#include "DetGeo.h"
#include <TFile.h>
#include <TTree.h>
#include <TBranch.h>
#include <TCanvas.h>
#include <TPad.h>
#include <TH1D.h>
#include <TH2D.h>
#include <vector>
#include <fstream>
#include <iostream>
#include <cmath>
#include <TMinuit.h>
#include <string>
#include <TStopwatch.h>

using namespace std;

class FullReco {
public:
  FullReco(); // Constructor
  double fullreco_total_error;
  struct selection {
    int layer;
    string type;
    double power;
    int position;
    int alt_position;
  };
  int initialize(string matrix_name, string prime_matrix_name); // takes two matrices to initialize
  void getSelection(vector<selection>& sel, vector<CubeE> hits);//TODO returns an order list of selections
  void processBlock(vector<CubeE>& hits, vector<int> layer, double ADC[2*N_PMT_LAYERS][N_TOTAL_PMTS]); // process block of layers
  int processFull(double ADC[2*N_PMT_LAYERS][N_TOTAL_PMTS], vector<CubeE>& vFullRecoHits); //TODO add energy starting from the block matrix definition
  int processPredefined(double ADC[2*N_PMT_LAYERS][N_TOTAL_PMTS], int nhits, vector<int> layer, vector<int> cubePosition, vector<double> initEnergy); //minimizes energies in locations provided from the layer definition

private:
  static vector<double>* TransferMatrixVec; // formerly xVecPtr
  static vector<double>* PrimeMatrixVec; // tracks the prime transfer matrix separately
  static vector<double>* PMTVec; // formerly xVecPtr2
  static vector<int>* hitCubeLayers; // tracks the cube layer that has an energy deposit
  static vector<int>* cubeIndex; // formerly cubeNumber
  static void fcn(int& npar, double* deriv, double& f, double par[], int flag);
};

#endif
