#include "EnergyReader.h"

EnergyReader::EnergyReader() {};

void EnergyReader::ReadEnergy(int evID, vector<int>* vADCData, vector<int>* vTrueHitCubeN, vector<int>* vTrueHitZ, vector<double>* vTrueHitEnergy) {
   
  double responseMatrix[N_CUBES][N_SIDES][N_PMTS]; // response matrices for one layer for each side
  double energyArray[N_CUBES][N_CUBE_LAYERS]; // vector with deposited energies
  double peArray[N_SIDES][N_PMTS][N_CUBE_LAYERS]; // vector with PMTs
  double summedPE[N_PMT_LAYERS][N_SIDES][N_PMTS];

  for (int s = 0; s < N_SIDES; s++) {
    for (int j = 0; j < N_PMTS; j++) {
      for (int k = 0; k < N_CUBE_LAYERS; k++) {
        peArray[s][j][k] = 0.0;
      }
    }
  }

  for (int s = 0; s < N_SIDES; s++) {
    for (int m = 0; m < N_PMT_LAYERS; m++) {
      for (int j = 0; j < N_PMTS; j++) {
        summedPE[s][m][j] = 0.0;
      }
    }
  }

  for (int i = 0; i < N_CUBE_LAYERS; i++) {
    for (int c = 0; c < N_CUBES; c++) {
      energyArray[c][i] = 0.0;
    }
  }

  // Preparing and reading the matrices
  string path = "/home/Events Reconstruction/Matrix 64x256 scaled.dat";
  int nlines = 0;
  ifstream file;
  file.open(path);
  if (file.is_open()) cout << "File successfully opened." << endl;
  else if (!file.is_open()) cout << "File failed to open." << endl;
  while (nlines < N_SIDES*N_PMTS) {
    for (int i = 0; i < N_CUBES; i++) {
      double responseMatrixValue = 0.0;
      file >> responseMatrixValue;
      for (int s = 0; s < N_SIDES; s++) {
        if (s*N_PMTS <= nlines && nlines < (s+1)*N_PMTS) responseMatrix[i][s][nlines-s*N_PMTS] = responseMatrixValue; // x1, x2, y1, y2
      }
    }
    nlines++;
  }  
  file.close();
  // Matrices are prepared

  string path2 = "/home/Events Reconstruction/Result.dat";
  int nlines2 = 0;
  ifstream file2;
  file2.open(path2);
  if (file2.is_open()) cout << "File 2 successfully open" << endl;
  else if (!file2.is_open()) cout << "File 2 failed to open" << endl;
  while (nlines2 < 34165) {
    int x, y, z;
    int eventID;
    double energy = 0.0;
    file2 >> eventID >> energy >> x >> y >> z;
    if (eventID == evID) {
      energyArray[x + N_PMTS * y][z] = energy;
      vTrueHitEnergy->push_back(energy);
      vTrueHitCubeN->push_back(x+N_PMTS*y);
      vTrueHitZ->push_back(z);
    }
    nlines2++;
  }
  file2.close();
  
  // Multiplication of the matrices and the energy array to get pe arrays
  for (int k = 0; k < N_CUBE_LAYERS; k++) { // Looking at cube layer k
    for (int s = 0; s < N_SIDES; s++) { // Looking at particular side of the detector with PMTs
      for (int j = 0; j < N_PMTS; j++) { // Looking at particular channel PMT in that layer on that side
        for (int i = 0; i < N_CUBES; i++) {
          peArray[s][j][k] = peArray[s][j][k] + responseMatrix[i][s][j] * energyArray[i][k];
        }
      }
    }
  }
  // Matrices are multiplied
  
  
  // Summing PEs to get the values in each PMT along vertical axis
  for (int s = 0; s < N_SIDES; s++) {
    for (int j = 0; j < N_PMTS; j++) {
      for (int m = 0; m < N_PMT_LAYERS; m++) {
        if (s == 0 || s == 2) {
          if (m == 0) summedPE[m][s][j] = peArray[s][j][m];
          else summedPE[m][s][j] = peArray[s][j][2*m-1] + peArray[s][j][2*m];
        }
        if (s == 1 || s == 3) {
          if (m == N_PMT_LAYERS - 1) summedPE[m][s][j] = peArray[s][j][2*m];
          else summedPE[m][s][j] = peArray[s][j][2*m] + peArray[s][j][2*m+1];
        }
      }
    }
  }
  

  // Throwing a Poisson around each expectation value
  TRandom3 *randPoisson = new TRandom3();
  int peValue[N_PMT_LAYERS][N_SIDES][N_PMTS];
  for (int m = 0; m < N_PMT_LAYERS; m++) {
    for (int s = 0; s < N_SIDES; s++) {
      for (int j = 0; j < N_PMTS; j++) {
        double poissonMean = summedPE[m][s][j];
        if (summedPE[m][s][j] == 0) peValue[m][s][j] = 0;
        else if (summedPE[m][s][j] > 0) peValue[m][s][j] = round(randPoisson->Poisson(poissonMean)); // round(poissonMean);
      }
    }
  }
  // End of throwing Poisson

  // Throwing Gaussian around each poisson PE integer value
  TRandom3 *randGaus = new TRandom3();
  int adcValue[N_PMT_LAYERS][N_SIDES][N_PMTS];
  for (int m = 0; m < N_PMT_LAYERS; m++) {
    for (int s = 0; s < N_SIDES; s++) {
      for (int j = 0; j < N_PMTS; j++) {
        if (peValue[m][s][j] == 0) adcValue[m][s][j] = 0;
        else if (peValue[m][s][j] > 0) {
          double adcMean = PE_TO_ADC * peValue[m][s][j];
          double adcSigma = sqrt(adcMean);
          int adcPositive = round(randGaus->Gaus(adcMean, adcSigma)); // = round(adcMean);
          if (adcPositive >= 0) adcValue[m][s][j] = adcPositive;
          else if (adcPositive < 0) adcValue[m][s][j] = 0;
        }
        vADCData->push_back(adcValue[m][s][j]);
      }
    }
  }
  cout << "Energy has been read, and response has been simulated" << endl;
}
