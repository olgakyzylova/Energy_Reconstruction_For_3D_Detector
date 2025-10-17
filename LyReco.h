#ifndef LyReco_h 
#define LyReco_h 1
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

using namespace std;

class LyReco {
public:
  LyReco();
  int initialize(string matrix_name);
  int processLayer(double ADC[], int recoPosition[], double recoEnergy[]);

private:
  static vector<double>* xVecPtr;
  static vector<double>* xVecPtr2;
  static vector<int>* cubeNumber;
  static void fcn(int& npar, double* deriv, double& f, double par[], int flag);
};

#endif

