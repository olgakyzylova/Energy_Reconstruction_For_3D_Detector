#include "FullReco.h"

vector<double>* FullReco::TransferMatrixVec = new vector<double>();
vector<double>* FullReco::PrimeMatrixVec = new vector<double>();
vector<double>* FullReco::PMTVec = new vector<double>();
vector<int>* FullReco::hitCubeLayers = new vector<int>(); // Tracks the cube layer which has ADC deposits
vector<int>* FullReco::cubeIndex = new vector<int>();

const double CUTOFF = 0.02; // 20 keV or more
const double ADC_THRESHOLD = 5; 

//++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++

void FullReco::fcn(int& npar, double* deriv, double& f, double par[], int flag){

  double sum = 0;
  vector<double> pos_map = *TransferMatrixVec;
  vector<double> prime_map = *PrimeMatrixVec;
  vector<double> PMT_Value = *PMTVec;
  vector<int> layerN = *hitCubeLayers;
  vector<int> cubeN = *cubeIndex;

  Double_t model[2*N_PMT_LAYERS][N_TOTAL_PMTS];
  memset(model, 0, sizeof(model)); // initialize model entries to 0

  for (int k = 0; k < npar; k++){
    for (int j = 0; j < N_TOTAL_PMTS; j++) {
      if((layerN[k] & 1) == 0){ // even layers
        model[layerN[k]][j] += pos_map[cubeN[k] + N_CUBES * j] * par[k];
        model[layerN[k]+1][j] += prime_map[cubeN[k] + N_CUBES * j] * par[k];
      }
      if((layerN[k] & 1) == 1){ // odd layers
        model[layerN[k]][j] += prime_map[cubeN[k] + N_CUBES * j] * par[k];
        model[layerN[k]+1][j] += pos_map[cubeN[k] + N_CUBES * j] * par[k];
      }
    }
  }

  for (int Ly = 0; Ly < 2*N_PMT_LAYERS; Ly++){
    for (int j = 0; j < N_TOTAL_PMTS; j++) {
      if(PMT_Value[Ly*N_TOTAL_PMTS+j] < ADC_THRESHOLD and PE_TO_ADC*model[Ly][j] < ADC_THRESHOLD) continue; // ADC threshold
      double model_ADC = PE_TO_ADC * model[Ly][j];
      double true_ADC = PMT_Value[Ly * N_TOTAL_PMTS + j]; // actually an int
      if (true_ADC < ADC_THRESHOLD and model_ADC < ADC_THRESHOLD) continue; // ADC threshold cut
      sum += model_ADC - true_ADC;
      if (true_ADC > 0) {
        if (model_ADC == 0) sum += pow(true_ADC, 5);
        else sum += true_ADC * log(true_ADC/model_ADC);
      }
    }
  }
  f = sum;
}

//++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++

FullReco::FullReco() { }

int FullReco::initialize(string matrix_name, string prime_matrix_name) {

  ifstream f_matrix(matrix_name);
  TransferMatrixVec->clear();
  if (!f_matrix.good()) cout << "Transfer Matrix File failed! Quitting..." << endl;
  double M_ij;
  for (int j = 0; j < N_TOTAL_PMTS; j++) {
    for (int i = 0; i < N_CUBES; i++) {
      f_matrix >> M_ij;
      TransferMatrixVec->push_back(M_ij); 
    }
  }
  f_matrix.close();

  // Also opens the prime matrix
  ifstream f_prime_matrix(prime_matrix_name);
  PrimeMatrixVec->clear();
  if (!f_prime_matrix.good()) cout << "Prime Matrix File failed! Quitting..." << endl;
  double pM_ij;
  for (int j = 0; j < N_TOTAL_PMTS; j++) {
    for (int i = 0; i < N_CUBES; i++) {
      f_prime_matrix >> pM_ij;
      PrimeMatrixVec->push_back(pM_ij); 
    }
  }
  f_prime_matrix.close();
  fullreco_total_error = 0;
  return 0;
}

//++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++
// GenerateResidue
void FullReco::getSelection(vector<selection>& sel, vector<CubeE> hits) {

  vector<double> pos_map = *TransferMatrixVec;
  vector<double> prime_map = *PrimeMatrixVec;
  vector<double> PMT_Value = *PMTVec;
  Double_t res[2*N_PMT_LAYERS][N_TOTAL_PMTS]; // 42 x 32; residue in PE
  for (int i = 0; i < 2*N_PMT_LAYERS; i++) { // 42
    for (int j = 0; j < N_TOTAL_PMTS; j++) { // 32
      res[i][j]=PMT_Value[i*N_TOTAL_PMTS+j]/PE_TO_ADC;
    }
  }
  for (int h = 0; h < hits.size(); h++) { // Looping over hits vector
    CubeE hit = hits.at(h);
    int z = hit.z;
    for (int j = 0; j < N_TOTAL_PMTS; j++) {
      if ((z & 1) == 0) { // even layers
        res[z][j] -= pos_map[hit.cubeN + N_CUBES * j] * hit.e;
        res[z+1][j] -= prime_map[hit.cubeN + N_CUBES * j] * hit.e;
      }
      if ((z & 1) == 1) { // odd layers
        res[z][j] -= prime_map[hit.cubeN + N_CUBES * j] * hit.e;
        res[z+1][j] -= pos_map[hit.cubeN + N_CUBES * j] * hit.e;
      }
    }
  }

  // loop over plausible hit locations
  bool possibleX[N_CUBE_LAYERS][N_PMTS];
  memset(possibleX, 0, sizeof(possibleX));
  bool possibleY[N_CUBE_LAYERS][N_PMTS];
  memset(possibleY, 0, sizeof(possibleY));
  int XinLayer[N_CUBE_LAYERS];
  memset(XinLayer, 0, sizeof(XinLayer));
  int YinLayer[N_CUBE_LAYERS];
  memset(YinLayer, 0, sizeof(YinLayer));
  for (int Ly = 0; Ly < N_CUBE_LAYERS; Ly++) {
    for (int p = 0; p < N_TOTAL_PMTS; p++) {
      // Check above threshold as standard
      if (PMT_Value[Ly*N_TOTAL_PMTS+p] > ADC_THRESHOLD and PMT_Value[(Ly+1)*N_TOTAL_PMTS+p] > ADC_THRESHOLD) {
        if (p < N_PMTS) {
          possibleX[Ly][p] = true;
          XinLayer[Ly]+=1;
        }
        if (p >= N_PMTS){
          possibleY[Ly][p-N_PMTS] = true;
          YinLayer[Ly]+=1;
        }
      }
    }
  }
  //Generate Selection Matrices
  //ADD MATRIX
  for (int Ly = 0; Ly < N_CUBE_LAYERS; Ly++){
    int search_size = XinLayer[Ly]*YinLayer[Ly];
    if(search_size==0) continue;
    for (int x = 0; x < N_PMTS; x++){
      if (not possibleX[Ly][x]) continue; // skip small
      for (int y = 0; y < N_PMTS; y++) {
        if (not possibleY[Ly][y]) continue; // skip small
                int cube_index = y * N_PMTS + x;
                double norm=0;
                double power=0;
                for (int j = 0; j < N_TOTAL_PMTS; j++) {
                    norm += pos_map[cube_index + N_CUBES * j]*pos_map[cube_index + N_CUBES * j];
                    norm += prime_map[cube_index + N_CUBES * j]*prime_map[cube_index + N_CUBES * j];
                    if((Ly & 1) == 0){ // even layers
                      power += res[Ly][j] * pos_map[cube_index + N_CUBES * j];
                        power += res[Ly+1][j] * prime_map[cube_index + N_CUBES * j];
                    }
                    if((Ly & 1) == 1){ // odd layers
                      power += res[Ly][j] * prime_map[cube_index + N_CUBES * j];
                        power += res[Ly+1][j] * pos_map[cube_index + N_CUBES * j];
                    }
                }
                power/=norm;
                selection new_selection{Ly, "ADD", power, cube_index, 0};
                sel.push_back(new_selection);
            }
    }
  }
  //SPLIT MATRIX
  for (int h = 0; h < hits.size(); h++) {
    CubeE hit = hits.at(h);
    if (hit.e < .1) continue; //skip splits on hits <100keV
    int z = hit.z;
    int offset[] = {0, 1, 1, 1, 0, -1, -1, -1};
    for(int s = 0; s < 8; s++) {
      //if(s>3){continue;}//adjacent
      //int x = hit.x + (s&1)*(1-(s&2));//adjacent
      //int y = hit.y + (1-(s&1))*(1-(s&2));//adjacent
      int x = hit.x + offset[s];//diagonal
      int y = hit.y + offset[(s+2)%8];//diagonal
      int parent_index = hit.y * N_PMTS + hit.x;
      int child_index = y * N_PMTS + x;
      if (x <= N_PMTS and y <= N_PMTS and x >= 0 and y >= 0) {
        double norm = 0;
        double power = 0;
        for (int j = 0; j < N_TOTAL_PMTS; j++) {
          double split = pos_map[child_index + N_CUBES * j]-pos_map[parent_index + N_CUBES * j];
          norm += split*split;
          split = prime_map[child_index + N_CUBES * j]-prime_map[parent_index + N_CUBES * j];
          norm += split*split;
          if ((z & 1) == 0) { // even layers
            power -= res[z][j] * pos_map[parent_index + N_CUBES * j];
            power += res[z][j] * pos_map[child_index + N_CUBES * j];
            power -= res[z+1][j] * prime_map[parent_index + N_CUBES * j];
            power += res[z+1][j] * prime_map[child_index + N_CUBES * j];
          }
          if ((z & 1) == 1) { // odd layers
            power -= res[z][j] * prime_map[parent_index + N_CUBES * j];
            power += res[z][j] * prime_map[child_index + N_CUBES * j];
            power -= res[z+1][j] * pos_map[parent_index + N_CUBES * j];
            power += res[z+1][j] * pos_map[child_index + N_CUBES * j];
          }
        }
        power /= norm;
        selection new_selection{z, "SPLIT", power, parent_index, child_index};
        //sel.push_back(new_selection);//disable split selection for now
      }//end if
    }
  }//end split selection
  //SWAP MATRIX
  for (int a = 0; a < hits.size(); a++) {
    CubeE hitA = hits.at(a);
    int z = hitA.z;
    for (int b = a + 1; b < hits.size(); b++) {
      CubeE hitB = hits.at(b);
      int A_index = hitA.y*N_PMTS+hitA.x;
      int B_index = hitB.y*N_PMTS+hitB.x;
      int sA_index = hitB.y*N_PMTS+hitA.x; //swap index
      int sB_index = hitA.y*N_PMTS+hitB.x; //swap index
      if (hitA.x != hitB.x and hitA.y != hitB.y and hitA.z==hitB.z){
                double norm=0;
                double power=0;
                for (int j = 0; j < N_TOTAL_PMTS; j++) {
                    double swap = pos_map[sA_index + N_CUBES * j]+pos_map[sB_index + N_CUBES * j];
                    swap -= pos_map[A_index + N_CUBES * j]+pos_map[B_index + N_CUBES * j];
                    norm += swap*swap;
                    swap = prime_map[sA_index + N_CUBES * j]+prime_map[sB_index + N_CUBES * j];
                    swap -= prime_map[A_index + N_CUBES * j]+prime_map[B_index + N_CUBES * j];
                    norm += swap*swap;
                    if((z & 1) == 0){ // even layers
                        power += res[z][j] * pos_map[sA_index + N_CUBES * j];
                        power += res[z][j] * pos_map[sB_index + N_CUBES * j];
                        power -= res[z][j] * pos_map[A_index + N_CUBES * j];
                        power -= res[z][j] * pos_map[B_index + N_CUBES * j];
                        power += res[z+1][j] * prime_map[sA_index + N_CUBES * j];
                        power += res[z+1][j] * prime_map[sB_index + N_CUBES * j];
                        power -= res[z+1][j] * prime_map[A_index + N_CUBES * j];
                        power -= res[z+1][j] * prime_map[B_index + N_CUBES * j];
                    }
                    if((z & 1) == 1){ // odd layers
                        power += res[z][j] * prime_map[sA_index + N_CUBES * j];
                        power += res[z][j] * prime_map[sB_index + N_CUBES * j];
                        power -= res[z][j] * prime_map[A_index + N_CUBES * j];
                        power -= res[z][j] * prime_map[B_index + N_CUBES * j];
                        power += res[z+1][j] * pos_map[sA_index + N_CUBES * j];
                        power += res[z+1][j] * pos_map[sB_index + N_CUBES * j];
                        power -= res[z+1][j] * pos_map[A_index + N_CUBES * j];
                        power -= res[z+1][j] * pos_map[B_index + N_CUBES * j];
                    }
                }
                power/=norm;
                selection new_selection{z, "SWAP", power, A_index, B_index};
                sel.push_back(new_selection);
            }
        }
    }//end swap selection
  sort(sel.begin(), sel.end(), [](selection a,selection b) {return a.power > b.power;});
    cout<<"DEBUG SELECTION"<<endl;//Debug
    for(int si=0;si<sel.size();si++){//TODO
        if(si>4){continue;}
        selection s = sel[si];
        cout<<si<<":    "<<s.type<<"    "<<"("<<s.position%N_PMTS<<","<<s.position/N_PMTS<< ","<<s.layer<<","<<s.power<<")";
        if(s.type=="SPLIT"){
        cout<<"--->("<<s.alt_position%N_PMTS<<","<<s.alt_position/N_PMTS<< ","<<s.layer<<","<<s.power<<")";
        }
        if(s.type=="SWAP"){
        cout<<"<==>("<<s.alt_position%N_PMTS<<","<<s.alt_position/N_PMTS<< ","<<s.layer<<","<<s.power<<")";
        }
        cout<<endl;
    }
}


//++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++

void FullReco::processBlock(vector<CubeE>& hits, vector<int> layer, double ADC[2*N_PMT_LAYERS][N_TOTAL_PMTS]) {

  // Condition ADC so that error is measured in blocks
  int start = layer.front();
  int end = layer.back()+1;
  PMTVec->clear();
  for (int z = 0; z < (2 * N_PMT_LAYERS); z++) {
    for (int p = 0; p < N_TOTAL_PMTS; p++){
      if (z >= start and z <= end) PMTVec->push_back(ADC[z][p]);
      else  PMTVec->push_back(0); // ignore other blocks
    }
  }

  // Minuit initialization
  Double_t amin, edm, errdef;
  Int_t nvpar, nparx, icstat;
  double amin_min;

  const int MAX_SELECTIONS = 20; // move to DetGeo.h
  amin_min = 1e100;
  for (int k = 0; k < MAX_SELECTIONS; k++) {
    vector<CubeE> best_hits = hits;
		vector<selection> select;
		getSelection(select, best_hits);
		int selection_type=-1;//track which operation was done
		int search_limit = select.size();
		bool priority_operation=false;//TODO
		if(select.size()>N_TOTAL_PMTS){search_limit=N_TOTAL_PMTS;}//search most likely
		//if(select.size()>N_TOTAL_PMTS){search_limit=10;}//search most likely
		for(int s=0;s<search_limit;s++){
			vector<CubeE> test_hits = hits;
			selection sel = select[s];
      TMinuit* minuit = 0;
      minuit = new TMinuit(hits.size()+2);//is this safe?
      minuit->SetFCN(FullReco::fcn);
      minuit->SetPrintLevel(-1);
			//Add new test hit
			if(sel.type == "ADD" and not priority_operation){
				int x = sel.position%N_PMTS;
				int y = sel.position/N_PMTS;
				int z = sel.layer;
				double selEnergy = sel.power;
		    hitCubeLayers->clear();
		    cubeIndex->clear();
		    for (int h = 0; h < hits.size(); h++) {//add old list
		    	CubeE prev_hit = hits[h];
		      if ((prev_hit.x == x) and (prev_hit.y == y) and (prev_hit.z == z)){
						goto goto_delete_minuit;
					}
		      minuit->DefineParameter(h, Form("Var%d",h), prev_hit.e, .01, 0, 8.333);
		      hitCubeLayers->push_back(prev_hit.z);
		      cubeIndex->push_back(prev_hit.cubeN);
		    }
		    minuit->DefineParameter(hits.size(), "TestParameter", selEnergy, .01, 0, 8.333);
		    hitCubeLayers->push_back(z);
		    cubeIndex->push_back(sel.position);
				//evaluate outcome
		    minuit->Migrad();
		    minuit->mnstat(amin, edm, errdef, nvpar, nparx, icstat);
				//cout<<"ADD: "<<sel.position<<"	"<<amin_min-amin<<endl;
		    if (amin < amin_min) {
		    	amin_min = amin;
		    	double energy,error;
		    	for (int s = 0; s < hits.size(); s++) {
		    		minuit->GetParameter(s, energy, error);
		    		test_hits[s].e = energy;
		      }
		      minuit->GetParameter(hits.size(), energy, error);
		      CubeE new_hit(sel.position%N_PMTS,sel.position/N_PMTS,sel.layer,sel.position,energy);
		      test_hits.push_back(new_hit);
		      best_hits = test_hits;
					selection_type = 1;
		    }
			}//end ADD
			//Split test hits
			if(sel.type == "SPLIT"){
				int px = sel.position%N_PMTS;
				int py = sel.position/N_PMTS;
				int z = sel.layer;
				double splitE = sel.power;
		    hitCubeLayers->clear();
		    cubeIndex->clear();
		    for (int h = 0; h < hits.size(); h++) {//add old list
		    	CubeE prev_hit = hits[h];
		      if ((prev_hit.x == px) and (prev_hit.y == py) and (prev_hit.z == z)){
						minuit->DefineParameter(h,"Parent", prev_hit.e-splitE, .01, 0, 8.333);//take energy from parent
				    hitCubeLayers->push_back(prev_hit.z);
				    cubeIndex->push_back(prev_hit.cubeN);
					}else{
				    minuit->DefineParameter(h, Form("Var%d",h), prev_hit.e, .01, 0, 8.333);
				    hitCubeLayers->push_back(prev_hit.z);
				    cubeIndex->push_back(prev_hit.cubeN);
					}
		    }
				//Add new test hit
		    minuit->DefineParameter(hits.size(), "Child", splitE, .01, 0, 8.333);
		    hitCubeLayers->push_back(z);
		    cubeIndex->push_back(sel.alt_position);//use child position
				//evaluate outcome
		    minuit->Migrad();
		    minuit->mnstat(amin, edm, errdef, nvpar, nparx, icstat);
				//cout<<"SPLIT: "<<sel.position<<"->"<<sel.alt_position<<"	"<<amin_min-amin<<endl;
		    if (amin < amin_min and amin/amin_min<.99) {
		    	amin_min = amin;
		    	double energy,error;
		    	for (int s = 0; s < hits.size(); s++) {
		    		minuit->GetParameter(s, energy, error);
		    		test_hits[s].e = energy;
		      }
		      minuit->GetParameter(hits.size(), energy, error);
		      CubeE new_hit(sel.alt_position%N_PMTS,sel.alt_position/N_PMTS,sel.layer,sel.alt_position,energy);
		      test_hits.push_back(new_hit);
		      best_hits = test_hits;
					selection_type = 2;
		    }
			}//end SPLIT
			//Swap test hits
			if(sel.type == "SWAP"){
				int ax = sel.position%N_PMTS;
				int ay = sel.position/N_PMTS;
				int bx = sel.alt_position%N_PMTS;
				int by = sel.alt_position/N_PMTS;
				int z = sel.layer;
				bool pre_swapA=true;//is there no hit on the swap?
				bool pre_swapB=true;//is there no hit on the swap?
				double swapE = sel.power;
		    hitCubeLayers->clear();
		    cubeIndex->clear();
		    for (int h = 0; h < hits.size(); h++) {//add old list
		    	CubeE prev_hit = hits[h];
		      if ((prev_hit.x == ax) and (prev_hit.y == ay) and (prev_hit.z == z)){
						minuit->DefineParameter(h,"PivotA", prev_hit.e-swapE, .01, 0, 8.333);//take energy from pivot
				    hitCubeLayers->push_back(prev_hit.z);
				    cubeIndex->push_back(prev_hit.cubeN);
					}
		      else if ((prev_hit.x == bx) and (prev_hit.y == by) and (prev_hit.z == z)){
						minuit->DefineParameter(h,"PivotB", prev_hit.e-swapE, .01, 0, 8.333);//take energy from pivot
				    hitCubeLayers->push_back(prev_hit.z);
				    cubeIndex->push_back(prev_hit.cubeN);
					}
		      else if ((prev_hit.x == ax) and (prev_hit.y == by) and (prev_hit.z == z)){//swap with preexisting
						minuit->DefineParameter(h,"SwapA", prev_hit.e+swapE, .01, 0, 8.333);
				    hitCubeLayers->push_back(prev_hit.z);
				    cubeIndex->push_back(prev_hit.cubeN);
						pre_swapA=false;
					}
		      else if ((prev_hit.x == bx) and (prev_hit.y == ay) and (prev_hit.z == z)){//swap with preexisting
						minuit->DefineParameter(h,"SwapB", prev_hit.e+swapE, .01, 0, 8.333);
				    hitCubeLayers->push_back(prev_hit.z);
				    cubeIndex->push_back(prev_hit.cubeN);
						pre_swapB=false;
					}else{
				    minuit->DefineParameter(h, Form("Var%d",h), prev_hit.e, .01, 0, 8.333);
				    hitCubeLayers->push_back(prev_hit.z);
				    cubeIndex->push_back(prev_hit.cubeN);
					}
		    }
				if(pre_swapA){
					minuit->DefineParameter(hits.size(),"SwapA", swapE, .01, 0, 8.333);//take energy from pivot
					hitCubeLayers->push_back(z);
					cubeIndex->push_back(by*N_PMTS+ax);
				}
				if(pre_swapB){
					minuit->DefineParameter(hits.size()+1,"SwapB", swapE, .01, 0, 8.333);//take energy from pivot
					hitCubeLayers->push_back(z);
					cubeIndex->push_back(ay*N_PMTS+bx);
				}
				//evaluate outcome
		    minuit->Migrad();
		    minuit->mnstat(amin, edm, errdef, nvpar, nparx, icstat);
		    if (amin < amin_min and amin/amin_min<.99) {
		    	amin_min = amin;
		    	double energy,error;
		    	for (int s = 0; s < hits.size(); s++) {
		    		minuit->GetParameter(s, energy, error);
		    		test_hits[s].e = energy;
		      }
					if(pre_swapA){
				    minuit->GetParameter(hits.size(), energy, error);
				    CubeE A_hit(ax,by,z,by*N_PMTS+ax,energy);
				    test_hits.push_back(A_hit);
					}
					if(pre_swapB){
				    minuit->GetParameter(hits.size()+1, energy, error);
				    CubeE B_hit(bx,ay,z,ay*N_PMTS+bx,energy);
				    test_hits.push_back(B_hit);
					}
		      best_hits = test_hits;
					selection_type = 3;
					//priority_operation = true;
		    }
			}//end SWAP
      goto_delete_minuit: //GOTO
      delete minuit;
			//if(priority_operation){break;}//kill selection on priority operations
		}//end selection loop

    // HALTING
		if(selection_type==-1){//no selection
		    cout << "NO SELECTION" << endl;
		    goto goto_end_selection;
		}
		if(selection_type==1){//ADD
		  if (best_hits.size() == hits.size()) {
		    cout << "NO NEW LOCATIONS" << endl;
		    goto goto_end_selection;
		  }
		  CubeE add_hit=best_hits[hits.size()];
		  if (add_hit.e < CUTOFF){ // cutoff at 10 keV
		    cout << "LOW ENERGY! QUITTING " << add_hit.e << endl;
		    goto goto_end_selection;
		  }
		  else {
		    cout << "HIT ADDED (x,y,z,E): " << "(" << add_hit.x << "," << add_hit.y << "," << add_hit.z << "," << add_hit.e << ")" << endl;
		    hits = best_hits;
		  }
		}
		if(selection_type==2){//SPLIT
		  if (best_hits.size() == hits.size()) {
		    cout << "NO NEW LOCATIONS" << endl;
		    goto goto_end_selection;
		  }
		  CubeE add_hit=best_hits[hits.size()];
		  if (add_hit.e < CUTOFF){ // cutoff at 10 keV
		    cout << "LOW ENERGY! QUITTING " << add_hit.e << endl;
		    goto goto_end_selection;
		  }
		  else {
		    cout << "HIT SPLIT (x,y,z,E): " << "(" << add_hit.x << "," << add_hit.y << "," << add_hit.z << "," << add_hit.e << ")" << endl;
		    hits = best_hits;
		  }
		}
		if(selection_type==3){//SWAP
		  CubeE A_hit=best_hits[hits.size()];
			CubeE B_hit=best_hits[hits.size()+1];
		  cout << "HIT SWAP (x,y,z,E): " << "(" << A_hit.x << "," << A_hit.y << "," << A_hit.z << "," << A_hit.e << ")";
		  cout << "	(" << B_hit.x << "," << B_hit.y << "," << B_hit.z << "," << B_hit.e << ")"<<endl;
		  hits = best_hits;
		}
    cout << "size " << hits.size() << "    error " << amin_min << endl;
		//CULL small hits
  	sort(hits.begin(), hits.end(), [](CubeE a,CubeE b) {return a.e > b.e;});
		while(hits.size()>0 and hits.back().e < CUTOFF){
			hits.pop_back();
		}
  } // end of selection loop
	goto_end_selection:
	cout << "BLOCK size " << hits.size() << "    error " << amin_min << endl;
	fullreco_total_error+=amin_min;
}

//++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++

int FullReco::processFull(double ADC[2*N_PMT_LAYERS][N_TOTAL_PMTS], vector<CubeE>& vFullRecoHits) {

  // Timing
  TStopwatch timer;
  timer.Start();
  // initialize PMT vector and find layers with hits
  PMTVec->clear();
  vector<int> relevant_layer;// only loop over hit layers
  bool has_hit_in_prev_layer = false;
  for (int z = 0; z < (2 * N_PMT_LAYERS); z++) {
    bool has_hit_in_layer = false;
    for (int p = 0; p < N_TOTAL_PMTS; p++) {
      PMTVec->push_back(ADC[z][p]);
      if (ADC[z][p] > 0) has_hit_in_layer = true;
    }
    if (has_hit_in_layer and has_hit_in_prev_layer) relevant_layer.push_back(z-1);
    has_hit_in_prev_layer = has_hit_in_layer;
  }

  vector<CubeE> hits;
  // This solves all at once
  // processBlock(hits,relevant_layer,ADC);
  // This solves in layer blocks
  if (relevant_layer.size() == 0) cout << "ERROR NO LAYERS FOUND!!!" << endl;
  vector<int> block_of_layers;
  vector<CubeE> hits_in_block;
  int prev_layer = relevant_layer[0];
  block_of_layers.push_back(prev_layer);
  for (int i = 1; i < relevant_layer.size(); i++){ // possibly merge with previous loop?
    int z = relevant_layer.at(i);
    if (z != (prev_layer + 1)) {
      processBlock(hits_in_block, block_of_layers, ADC);
      for (int h = 0; h < hits_in_block.size(); h++) {
        hits.push_back(hits_in_block.at(h));
      }
      hits_in_block.clear();
      block_of_layers.clear();
    }
    block_of_layers.push_back(z);
    prev_layer = z;
  }

  //Process final block
  processBlock(hits_in_block, block_of_layers, ADC);
  for (int h = 0; h < hits_in_block.size(); h++) {
    hits.push_back(hits_in_block.at(h));
  }

  //Print results
  timer.Stop();
  cout << "END MINIMIZATION at " << hits.size() << " hits in " << timer.RealTime() << "s" << endl;
  sort(hits.begin(), hits.end(), [](CubeE a,CubeE b) {return a.e > b.e;} );
  for (int s = 0; s < hits.size(); s++) {
    CubeE h = hits[s];
    if (h.e > CUTOFF) {
      vFullRecoHits.push_back(h);
      cout << "Result: z = " << h.z << ",    Cube # = " << h.y*N_PMTS+h.x << ",    x = " << h.x << ",    y = " << h.y << ",    E = (" << h.e*1000 <<" +/- "<<.06*sqrt(h.e)*1000 << ") keV" << endl;
    }
  }
  return 1;
}

//++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++

//minimizes energies in locations provided from the full definition
int FullReco::processPredefined(double ADC[2*N_PMT_LAYERS][N_TOTAL_PMTS], int nhits, vector<int> cubelayer, vector<int> cubePosition, vector<double> initEnergy){

//inputs:
//ADC: every pmt in the detector adc value in an array with pmt layer as the first index and pmt number as the second
//nhits: number of hits to minimize
//cubelayer: array with the cube layer of each test hit
//cubePosition: array with the cube position of each test hit
//initEnergy: array with starting energy values of each test hit
  vector<CubeE> finalHits;
  PMTVec->clear();
  hitCubeLayers->clear();
  cubeIndex->clear();
  for (int z = 0; z < (2 * N_PMT_LAYERS); z++) {
    for (int p = 0; p < N_TOTAL_PMTS; p++) {
      PMTVec->push_back(ADC[z][p]);
    }
  }
  //memcpy(PMT_ADC,ADC,sizeof(ADC));//copy to local array
  vector<double> Energy_result;
  vector<double> Energy_error;
  //Minuit
  Double_t amin, edm, errdef;
  Int_t nvpar, nparx, icstat;
  TMinuit* minuit = 0;
  minuit = new TMinuit(nhits);
  minuit->SetFCN(FullReco::fcn);
  minuit->SetPrintLevel(-1);
  double init_error = 0;
  int init_nhits = nhits;
  double init_par[100];
  for (int s = 0; s < nhits; s++) {
    minuit->DefineParameter(s, Form("Var%d",s), initEnergy[s], 0.01, 0, 8.333); //TODO bring over error instead of .01
    hitCubeLayers->push_back(cubelayer[s]);
    cubeIndex->push_back(cubePosition[s]);
    init_par[s] = initEnergy[s];
  }
  fcn(init_nhits, 0, init_error, init_par, 0);
  minuit->Migrad();
  minuit->mnstat(amin, edm, errdef, nvpar, nparx, icstat);
  cout << "AMIN    " << amin << endl;
  cout << "PREDEFINED Size " << nhits << "	error " << amin << endl;
  for (int s = 0; s < nhits; s++) {
    double energy, energy_err; 
    minuit->GetParameter(s, energy, energy_err); //get energy+error from nth parameter
    Energy_result.push_back(energy);
    Energy_error.push_back(energy_err);    
    cout << "Result: z = " << cubelayer[s] << ",    Cube # = " << cubePosition[s] << ",    x = " << cubePosition[s]%16 << ",    y = " << cubePosition[s]/16 << ",    E = (" << energy*1000 <<" +/- "<<.06*sqrt(energy)*1000 << ") keV" << endl;
  }
  return 1;
}
