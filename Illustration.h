#ifndef Illustration_h 
#define Illustration_h 1
#include "DetGeo.h"
#include <TCanvas.h>
#include <TPad.h>
#include <TAttFill.h>
#include <TLine.h>
#include <TLatex.h>
#include <TMarker.h>
#include <TBox.h>
#include <vector>
#include <fstream>
#include <iostream>
#include <cmath>
#include <vector>
#include <string>
#include <algorithm>

using namespace std;

class Illustration {
public:
  Illustration();
  void DrawTable();
  void Draw(vector<CubeE> trueHits, /*vector<CubeE> leftRecoHits, vector<CubeE> rightRecoHits, vector<CubeE> fullRecoHits,*/ vector<CubeE> vFullRecoHits, int evID);

private:
};

#endif
