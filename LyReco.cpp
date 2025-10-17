#include "LyReco.h"
#include "DetGeo.h"
#include <fstream>
#include <iostream>
#include <cmath>
#include <vector>
#include <TMinuit.h>
#include <string>

vector<double>* LyReco::xVecPtr = new vector<double>();
vector<double>* LyReco::xVecPtr2 = new vector<double>();
vector<int>* LyReco::cubeNumber = new vector<int>();

//++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++

void LyReco::fcn(int& npar, double* deriv, double& f, double par[], int flag){
  
  double sum = 0;
  vector<double> pos_map = *xVecPtr;
  vector<double> PMT_value = *xVecPtr2;
  vector<int> cubeN = *cubeNumber;

  Double_t model[N_TOTAL_PMTS];
  memset(model, 0, sizeof(model));

  for (int j = 0; j < N_TOTAL_PMTS; j++) {
    for (int k = 0; k < npar; k++) {
      model[j] += par[k] * pos_map[cubeN[k] + N_CUBES * j];
    }
  }

  for (int j = 0; j < N_TOTAL_PMTS; j++) {
    sum += PE_TO_ADC * model[j] - PMT_value[j];
//    cout << "o:" << PMT_value[j] << "  e:" << model[j] << endl;    
    if (PMT_value[j] != 0 && model[j] == 0) sum += pow(PMT_value[j], 5);
    if (PMT_value[j] != 0) {
      sum += PMT_value[j] * log(PMT_value[j]) - PMT_value[j] * log(PE_TO_ADC * model[j]);
    }
  }
//  cout << "loglikelihood:" << sum << endl;
  f = sum;
}

//++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++

LyReco::LyReco() { }

int LyReco::initialize(string matrix_name) {

  ifstream f_matrix(matrix_name);
  xVecPtr->clear();
  if (!f_matrix.good()) cout << "File failed! Quitting..." << endl;

  double M_ij;
  for (int j = 0; j < N_TOTAL_PMTS; j++) {
    for (int i = 0; i < N_CUBES; i++) {
      f_matrix >> M_ij;
      xVecPtr->push_back(M_ij); 
    }
  }
  f_matrix.close();
  return 0;
}

//++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++

int LyReco::processLayer(double ADC[], int recoPosition[], double recoEnergy[]) {

  //Initialize Minuit
  Double_t amin, edm, errdef;
  Int_t nvpar, nparx, icstat;

  double totPH = 0;
  double ADC_sum_x = 0;
  double ADC_sum_y = 0;

  for(int i = 0; i < N_PMTS; i++) {
    ADC_sum_x += ADC[i];
    totPH += ADC[i];
  }

  for (int i = N_PMTS; i < N_TOTAL_PMTS; i++) {
    ADC_sum_y += ADC[i];
    totPH += ADC[i];
  }

  if (ADC_sum_x == 0 || ADC_sum_y == 0) return 0; //skip zero ADC layers for x and y

  double amin_min, par_temp, par_temp_err;  
  double resultPosition[N_ITER]; // resulted position
  double result[N_ITER];
  double logL[N_ITER]; // resulted logL
  memset(result, 0, sizeof(result));
  memset(resultPosition, 0, sizeof(resultPosition));
  memset(logL, 0, sizeof(logL));

  for (int i = 0; i < N_ITER; i++) {
    resultPosition[i] = -1;
    result[i] = 0.0;
    logL[i] = 0.0;
    cubeNumber->push_back(-1);
  }
  xVecPtr2->clear();
  
  for (int j = 0; j < N_PMTS; j++) {
    xVecPtr2->push_back(ADC[j]);
//    cout << ADC[j] << "  ";
  }
  
  for (int j = N_PMTS; j < N_TOTAL_PMTS; j++) {
    xVecPtr2->push_back(ADC[j]);
//    cout << ADC[j] << "  ";
  }
//  cout << endl;

//+++++++++++++++++++++++ MINIMIZATION ++++++++++++++++++++++++++++++++++++++

  for (int k = 0; k < N_ITER; k++) {
    amin_min = 1e13;
    for (int i = 0; i < N_CUBES; i++) {
      int a = 1;
      for (int m = 0; m < k; m++) {
        if (i == resultPosition[m]) a = a * 0;
      }
      if (a == 0) continue;
      cubeNumber->at(k) = i;
      for (int s = 0; s < k; s++) {
        cubeNumber->at(s) = resultPosition[s];
      }
      
      TMinuit* minuit = 0;
      minuit = new TMinuit(k+1);
      minuit->SetFCN(LyReco::fcn);
      minuit->SetPrintLevel(-1);
     
      for (int s = 0; s <= k; s++) {
        minuit->DefineParameter(s, Form("AAA%d",s), 0.0, 0.01, 0, 8.333);
      }
      minuit->Migrad();
      
      minuit->mnstat(amin, edm, errdef, nvpar, nparx, icstat);
      if (amin < amin_min) {
        amin_min = amin;
        for (int s = 0; s<= k; s++) {
          minuit->GetParameter(s, result[s], par_temp_err);
          resultPosition[k] = i;
          logL[k] = amin;
        }
      }
      delete minuit;
    }
    cout << k << "	" << logL[k] << endl;
    if (logL[k] < 1 || result[k] < 0.016) {
      result[k] = 0;
      resultPosition[k] = -1;
      break;
    }
/*    if (result[k] < 0.016) {
      result[k] = 0;
      resultPosition[k] = -1;
      break;
    }*/
  } // end of iterations loop

  for (int k = 0; k < N_ITER; k++) {
    recoPosition[k] = resultPosition[k];
    recoEnergy[k] = result[k];
  }
  
//++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++
  return 1;
} // end of layer loop
