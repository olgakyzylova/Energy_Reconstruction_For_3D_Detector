#include "DetGeo.h"
#include "LyReco.h"
#include "FullReco.h"
#include "EnergyReader.h"
#include "Illustration.h"
#include <iostream>
using namespace std;

bool comparator(const CubeE& lhs, const CubeE& rhs) {
  return lhs.z < rhs.z;
}

int main(int argc, char* argv[]) {

  string eventNumber = "value";
  eventNumber = argv[1];
  cout << eventNumber << endl;
  int evtID = stoi(eventNumber);
  vector<int> vADCData, vTrueCubeN, vTrueZ;
  vector<double> vTrueEnergy;
  EnergyReader *eRead = new EnergyReader();
  eRead->ReadEnergy(evtID, &vADCData, &vTrueCubeN, &vTrueZ, &vTrueEnergy);

  // Getting simulated ADC values from EnergyReader.cpp
  int adcPMTData[N_PMT_LAYERS][N_SIDES][N_PMTS];
  for (int m = 0; m < N_PMT_LAYERS; m++) {
    for (int s = 0; s < N_SIDES; s++) {
      for (int j = 0; j < N_PMTS; j++) {
        adcPMTData[m][s][j] = vADCData.at(m*N_SIDES*N_PMTS + s*N_PMTS + j);
      }
    }
  }
  
  string matrix = "./Matrix 32x256 scaled.dat";
  string matrixPrime = "./Matrix 32x256 primed scaled.dat";

  vector<int> vResultZ, vResultCubeN, vPrimeResultZ, vPrimeResultCubeN, vFullResultZ, vFullResultCubeN;
  vector<double> vResultEnergy, vPrimeResultEnergy, vFullResultEnergy;

// ++++++++++++++++++++++++ FILLING OUT PMT DATA ARRAYS +++++++++++++++++++++++++++++++++++++  
  double PMTDataArray[N_PMT_LAYERS][N_TOTAL_PMTS];
  double PMTDataArrayPrime[N_PMT_LAYERS][N_TOTAL_PMTS];
  for (int k = 0; k < N_PMT_LAYERS; k++) {
    for (int j = 0; j < N_PMTS; j++) {
      for (int s = 0; s < N_SIDES; s++) {
        if (s == 0) PMTDataArray[k][j] = adcPMTData[k][s][j]; // x values
        if (s == 2) PMTDataArray[k][N_PMTS + j] = adcPMTData[k][s][j]; // y values
        if (s == 1) PMTDataArrayPrime[k][j] = adcPMTData[k][s][j]; // x' values
        if (s == 3) PMTDataArrayPrime[k][N_PMTS + j] = adcPMTData[k][s][j]; // y' values
      }
    }
  }
// ++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++
/*
  // Looping over the PMT layers
  for (int k = 0; k < N_PMT_LAYERS; k++) {

    int positionPrime[N_ITER];// = {0, 0, 0, 0, 0, 0, 0};
    double energyPrime[N_ITER];// = {0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0};
    for (int i = 0; i < N_ITER; i++) {
      positionPrime[i] = 0;
      energyPrime[i] = 0.0;
    }
    LyReco *recoPrime = new LyReco();
    recoPrime->initialize(matrixPrime);
    recoPrime->processLayer(PMTDataArrayPrime[k], positionPrime, energyPrime);
    if (energyPrime[0] > 0) cout << "++++++ Layer " << k << ", Right Side:" << endl;
    for (int m = 0; m < N_ITER; m++) {
      if (energyPrime[m] > 0) {
        cout << "Cube #        " << positionPrime[m] << ",  E = " << 1000 * energyPrime[m] << " keV" << endl;
        vPrimeResultZ.push_back(k);
        vPrimeResultCubeN.push_back(positionPrime[m]);
        vPrimeResultEnergy.push_back(energyPrime[m]);
      }
    }
    delete recoPrime;

    if (energyPrime[0] > 0) { // if there is energy in a layer
      int positionA[N_ITER] = {0, 0, 0, 0, 0, 0, 0};
      double energyA[N_ITER] = {0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0};
      LyReco *recoA = new LyReco();
      recoA->initialize(matrix);
      recoA->processLayer(PMTDataArray[k], positionA, energyA);
      delete recoA;
      if (energyA[0] > 0) cout << "++++++ Layer " << k << ", Left Side:" << endl;
      for (int m = 0; m < N_ITER; m++) {
        if (energyA[m] > 0) {
          cout << "Cube #      " << positionA[m] << ",  E = " << 1000 * energyA[m] << " keV" << endl;
          vResultZ.push_back(k);
          vResultCubeN.push_back(positionA[m]);
          vResultEnergy.push_back(energyA[m]);
          for (int m2 = 0; m2 < N_ITER; m2++) {
            if (positionA[m] == positionPrime[m2]) {
              vFullResultZ.push_back(k+k);
              vFullResultCubeN.push_back(positionA[m]);
              vFullResultEnergy.push_back((energyPrime[m2]+energyA[m])/2);
            }
          }
        }
      }

      int positionB[N_ITER] = {0, 0, 0, 0, 0, 0, 0};
      double energyB[N_ITER] = {0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0};
      LyReco *recoB = new LyReco();
      recoB->initialize(matrix);
      recoB->processLayer(PMTDataArray[k+1], positionB, energyB);
      delete recoB;
      if (energyB[0] > 0) cout << "++++++ Layer " << k+1 << ", Left Side:" << endl;
      for (int m = 0; m < N_ITER; m++) {
        if (energyB[m] > 0) {
          cout << "Cube #      " << positionB[m] << ",  E = " << 1000 * energyB[m] << " keV" << endl;
          vResultZ.push_back(k+1);
          vResultCubeN.push_back(positionB[m]);
          vResultEnergy.push_back(energyB[m]);
          for (int m2 = 0; m2 < N_ITER; m2++) {
            if (positionB[m] == positionPrime[m2]) {
              vFullResultZ.push_back(k+k+1);
              vFullResultCubeN.push_back(positionB[m]);
              vFullResultEnergy.push_back((energyPrime[m2]+energyB[m])/2);
            }
          }
        }
      }
    }
  } // End of the loop over layers
*/
// ++++++++++++++++++++++++ MEGA MATRIX FULL RECONSTRUCTION +++++++++++++++++++++++++++++++++

  cout << "TESTING FULL RECONSTRUCTION" << endl;
  FullReco *fullreco = new FullReco();
  fullreco->initialize(matrix, matrixPrime);
  std::vector<CubeE> vFullRecoHits;
  double fullPMTDataArray[2*N_PMT_LAYERS][N_TOTAL_PMTS];
  for (int k = 0; k < N_PMT_LAYERS; k++) { 
    for (int j = 0; j < N_PMTS; j++) {
      for (int s = 0; s < N_SIDES; s++) {
        if (s == 0) fullPMTDataArray[2*k][j] = adcPMTData[k][s][j];
        if (s == 2) fullPMTDataArray[2*k][N_PMTS + j] = adcPMTData[k][s][j];
        if (s == 1) fullPMTDataArray[2*k+1][j] = adcPMTData[k][s][j];
        if (s == 3) fullPMTDataArray[2*k+1][N_PMTS + j] = adcPMTData[k][s][j];
      }
    }
  }
//  for (int i = 0; i < vFullResultZ.size(); i++) {
//    cout << "AAA " << vFullResultZ[i] << "	" << vFullResultCubeN[i] << "	" << vFullResultEnergy[i] << endl;
//  }

//  fullreco->processPredefined(fullPMTDataArray, vFullResultZ.size(), vFullResultZ, vFullResultCubeN, vFullResultEnergy);
  fullreco->processFull(fullPMTDataArray, vFullRecoHits);
  sort(vFullRecoHits.begin(), vFullRecoHits.end(), &comparator);

// ++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++

  // Vector of objects of true hits
  std::vector<CubeE> trueHits;
  for (int w = 0; w < vTrueZ.size(); w++) {
    trueHits.push_back(CubeE(vTrueCubeN.at(w)%N_PMTS, vTrueCubeN.at(w)/N_PMTS, vTrueZ.at(w), vTrueCubeN.at(w), vTrueEnergy.at(w)));
  }

  // Printing true information
  cout << "TRUE ENERGY INFORMATION: " << endl;
  for (int w = 0; w < vTrueZ.size(); w++) {
    cout << "Event ID = " << evtID << ",	z = " << vTrueZ.at(w) << ",	Cube # = " << vTrueCubeN.at(w) << ",	x = ";
    cout << vTrueCubeN.at(w)%N_PMTS << ",	y = " << vTrueCubeN.at(w)/N_PMTS << ",	Energy = " << vTrueEnergy.at(w)*1000 << " keV" << endl;
  }
  sort(trueHits.begin(), trueHits.end(), &comparator);
//  fullreco->processPredefined(fullPMTDataArray, vTrueZ.size(), vTrueZ, vTrueCubeN, vTrueEnergy);
/*
  // Vector of objects of left reconstructed hits
  std::vector<CubeE> leftRecoHits;
  for (int w = 0; w < vResultZ.size(); w++) {
    leftRecoHits.push_back(CubeE(vResultCubeN.at(w)%N_PMTS, vResultCubeN.at(w)/N_PMTS, vResultZ.at(w), vResultCubeN.at(w), vResultEnergy.at(w)));
  }
  sort(leftRecoHits.begin(), leftRecoHits.end(), &comparator);

  // Vector of objects of right reconstructed hits
  std::vector<CubeE> rightRecoHits;
  for (int w = 0; w < vPrimeResultZ.size(); w++) {
    rightRecoHits.push_back(CubeE(vPrimeResultCubeN.at(w)%N_PMTS, vPrimeResultCubeN.at(w)/N_PMTS, vPrimeResultZ.at(w), vPrimeResultCubeN.at(w), vPrimeResultEnergy.at(w)));
  }
  sort(rightRecoHits.begin(), rightRecoHits.end(), &comparator);

  // Vector of objects of full reconstructed hits
  std::vector<CubeE> fullRecoHits;
  for (int w = 0; w < vFullResultZ.size(); w++) {
    fullRecoHits.push_back(CubeE(vFullResultCubeN.at(w)%N_PMTS, vFullResultCubeN.at(w)/N_PMTS, vFullResultZ.at(w), vFullResultCubeN.at(w), vFullResultEnergy.at(w)));
  }
  sort(fullRecoHits.begin(), fullRecoHits.end(), &comparator);

  // Printing reconstructed information
  cout << "RECONSTRUCTED ENERGY: " << endl;
  for (int d = 0; d < vFullResultZ.size(); d++) {
    cout << "z = " << vFullResultZ.at(d) << ",	Cube # = " << vFullResultCubeN.at(d) << ",	x = ";
    cout << vFullResultCubeN.at(d)%N_PMTS << ",	y = " << vFullResultCubeN.at(d)/N_PMTS << ",	Energy = " << vFullResultEnergy.at(d)*1000 << " keV" << endl;
  }
*/  
  Illustration *illustr = new Illustration();
/*  bool dodraw = false;
  if(vFullRecoHits.size() != trueHits.size()) dodraw=true;
  for(int t = 0; t < trueHits.size(); t++) {
    CubeE truth = trueHits[t];
    bool match=false;
    for(int r = 0; r < vFullRecoHits.size(); r++) {
      CubeE reco = vFullRecoHits[r];
      if (truth.cubeN == reco.cubeN and truth.z == reco.z) {
        match=true;
        break;
      }
    }
    if (not match) dodraw = true;
  }
  if(dodraw)*/
  illustr->Draw(trueHits, /*leftRecoHits, rightRecoHits, fullRecoHits,*/ vFullRecoHits, evtID);
  return 0;
}
