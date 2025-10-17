#include "Illustration.h"

Illustration::Illustration() { }

void Illustration::DrawTable() {

  for (int i = 0; i < 19; i++) { // Drawing the table
    auto l = new TLine(0.05*(i+1), 0.05, 0.05*(i+1), 0.95);
    l->SetLineColor(1);
    l->SetLineWidth(2);
    l->Draw();

    l = new TLine(0.05, 0.05*(i+1), 0.95, 0.05*(i+1));
    l->SetLineColor(1);
    l->SetLineWidth(2);
    l->Draw();
  }

  for (int i = 0; i < 16; i++) { // Drawing numbers
    auto tex = new TLatex(0.06, 0.865-0.05*i, Form("%i",i));
    tex->SetTextColor(1);
    tex->SetTextSize(0.03);
    tex->Draw();

    tex = new TLatex(0.11+0.05*i, 0.065, Form("%i", i));
    tex->SetTextColor(1);
    tex->SetTextSize(0.03);
    tex->Draw();

    tex = new TLatex(0.91, 0.865-0.05*i, Form("%i", i));
    tex->SetTextColor(1);
    tex->SetTextSize(0.03);
    tex->Draw();

    tex = new TLatex(0.11+0.05*i, 0.915, Form("%i", i));
    tex->SetTextColor(1);
    tex->SetTextSize(0.03);
    tex->Draw();
  }
}

void Illustration::Draw(vector<CubeE> trueHits, /*vector<CubeE> leftRecoHits, vector<CubeE> rightRecoHits, vector<CubeE> fullRecoHits,*/ vector<CubeE> vFullRecoHits, int evID) {

  // Draw the detector one-layer grid
  int a = -1;
  TCanvas *c = new TCanvas("c", Form("Layer %i", a), 224, 330, 900, 900);
  TPad *p = new TPad("p", "p", 0.0, 0.0, 1, 1, 0, -1, -1);
  for (int i = 0; i < trueHits.size(); i++) {
    if (trueHits[i].z != a) {
      c = new TCanvas("c", Form("Layer %i", trueHits[i].z), 224, 330, 900, 900);
      p = new TPad("p", "p", 0.0, 0.0, 1, 1, 0, -1, -1);
      p->Draw();
      p->cd();
      TLatex *tex = new TLatex(0.39, 0.97, Form("Event %i, Layer %i", evID, trueHits[i].z));
      tex->SetTextSize(0.03);
      tex->Draw();
      Illustration::DrawTable();
    }

    // Red and blue boxes for true hits from adjacent layers
/*    for (int j = 0; j < trueHits.size(); j++) {
      if (trueHits[j].z == trueHits[i].z - 1) {
        auto box = new TBox(0.102+trueHits[j].x*0.05, 0.852-trueHits[j].y*0.05, 0.148+trueHits[j].x*0.05, 0.898-trueHits[j].y*0.05);
        box->SetFillColor(kRed-10);
        box->Draw();
      }
      if (trueHits[j].z == trueHits[i].z + 1) {
        auto box = new TBox(0.102+trueHits[j].x*0.05, 0.852-trueHits[j].y*0.05, 0.148+trueHits[j].x*0.05, 0.898-trueHits[j].y*0.05);
        box->SetFillColor(kBlue-10);
        box->Draw();
      }
    }
*/
    // Yellow box for true hits
    auto box = new TBox(0.102+trueHits[i].x*0.05, 0.852-trueHits[i].y*0.05, 0.148+trueHits[i].x*0.05, 0.898-trueHits[i].y*0.05);
    if (trueHits[i].e < 0.025) box->SetFillColor(kYellow-8);
    if (0.025 <= trueHits[i].e && trueHits[i].e < 0.1) box->SetFillColor(kYellow-10);
    if (0.1 <= trueHits[i].e && trueHits[i].e < 0.5) box->SetFillColor(kYellow-9);
    if (0.5 <= trueHits[i].e && trueHits[i].e < 1) box->SetFillColor(kYellow-7);
    if (1 <= trueHits[i].e) box->SetFillColor(kYellow);
    box->Draw();

    // Pattern box for reconstruction
    for (int j = 0; j < vFullRecoHits.size(); j++) {
      if (trueHits[i].z == vFullRecoHits[j].z) {
        auto box = new TBox(0.1+vFullRecoHits[j].x*0.05, 0.85-vFullRecoHits[j].y*0.05, 0.15+vFullRecoHits[j].x*0.05, 0.9-vFullRecoHits[j].y*0.05);
        box->SetFillColor(1);
        box->SetFillStyle(3021);
        box->Draw();
      }
    }
/*
    // Markers for Left Reco Hits, blue
    for (int j = 0; j < leftRecoHits.size(); j++) {
      if ((trueHits[i].z %2 == 0) && (leftRecoHits[j].z == 0.5*trueHits[i].z)) {
        auto marker = new TMarker(0.1125+0.05*leftRecoHits[j].x, 0.8875-0.05*leftRecoHits[j].y, 20);
        marker->SetMarkerColor(kBlue);
        marker->SetMarkerStyle(20);
        marker->SetMarkerSize(1.7);
        marker->Draw();
      }
      if ((trueHits[i].z %2 == 1) && (leftRecoHits[j].z == 0.5*(trueHits[i].z+1))) {
        auto marker = new TMarker(0.1125+0.05*leftRecoHits[j].x, 0.8875-0.05*leftRecoHits[j].y, 20);
        marker->SetMarkerColor(kBlue);
        marker->SetMarkerStyle(20);
        marker->SetMarkerSize(1.7);
        marker->Draw();
      }
    }

    // Markers for Right Reco Hits, red
    for (int j = 0; j < rightRecoHits.size(); j++) {
      if ((trueHits[i].z %2 == 0) && (rightRecoHits[j].z == 0.5*trueHits[i].z)) {
        auto marker = new TMarker(0.1375+0.05*rightRecoHits[j].x, 0.8875-0.05*rightRecoHits[j].y, 20);
        marker->SetMarkerColor(kRed);
        marker->SetMarkerStyle(45);
        marker->SetMarkerSize(1.8);
        marker->Draw();
      }
      if ((trueHits[i].z %2 == 1) && (rightRecoHits[j].z == 0.5*(trueHits[i].z-1))) {
        auto marker = new TMarker(0.1375+0.05*rightRecoHits[j].x, 0.8875-0.05*rightRecoHits[j].y, 20);
        marker->SetMarkerColor(kRed);
        marker->SetMarkerStyle(45);
        marker->SetMarkerSize(1.8);
        marker->Draw();
      }
    }

    // Markers for Full Reco Hits, green
    for (int j = 0; j < fullRecoHits.size(); j++) {
      if (trueHits[i].z == fullRecoHits[j].z) {
        auto marker = new TMarker(0.125+0.05*fullRecoHits[j].x, 0.8625-0.05*fullRecoHits[j].y, 20);        
        marker->SetMarkerColor(kGreen+2);
        marker->SetMarkerStyle(29);
        marker->SetMarkerSize(1.9);
        marker->Draw();
      }
    }
*/
    c->SaveAs(Form("/home/Events Reconstruction/vent %i, Layer %i.png", evID, trueHits[i].z));
    a = trueHits[i].z;
  }

  ofstream file;
  file.open("/home/Events Reconstruction/SavedResultsNew.dat", ios::app);
  if (file.is_open()) cout << "File successfully opened." << endl;
  else if (!file.is_open()) cout << "File failed to open." << endl;
  for (int i = 0; i < trueHits.size(); i++) {
    file << evID << "	" << trueHits[i].cubeN << "	" <<  trueHits[i].x << "	" << trueHits[i].y << "	" << trueHits[i].z << "	" << trueHits[i].e << "	" << 0 << endl;
  }
  for (int i = 0; i < vFullRecoHits.size(); i++) {
    file << evID << "	" << vFullRecoHits[i].cubeN << "	" << vFullRecoHits[i].x << "	" << vFullRecoHits[i].y << "	" << vFullRecoHits[i].z << "	" << vFullRecoHits[i].e << "	" << 1 << endl;
  }
  file.close();
}
