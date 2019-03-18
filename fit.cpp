/*  fit.cpp
 *  Class for correcting focal plane aberrations using interpolation of several polynomials 
 *  Corrects for tilt in the focal plane along with the characteristic x|theta aberrations
 *  The base constructor gives 5 polynomials, can override to give as many as needed
 *  G.M. Feb 2019
 *  Revised March 2019 to run without reopening and closing files as shown by KGH -- G.M.
 */

#include "fit.h"
#include "TCanvas.h"
#include <iostream>
#include <string>
#include "TMath.h"
#include <string>

using namespace std;
//Constructors
fit::fit() :
  tilt_space(new TCutG("tilt_space", 0)),
  tilt(new TF1("tilt", "pol1")),
  nfuncs(5) //default is 5 polynomials
{
  f = new TF1*[nfuncs];
  f_space = new TCutG*[nfuncs];
  for (int i = 0; i<nfuncs; i++) {
    char fname[20];
    char f_spacename[20];
    strcpy(fname, Form("f%d", i));
    strcpy(f_spacename, Form("f_space%d", i));
    f[i] = new TF1(fname, "pol3");
    f_space[i] = new TCutG(f_spacename, 0);
    c.push_back(0.0);
  }
}

fit::fit(int n) : //n is number of fits
  tilt_space(new TCutG("tilt_space", 0)),
  tilt(new TF1("tilt", "pol1")),
  nfuncs(n)
{
  f = new TF1*[nfuncs];
  f_space = new TCutG*[nfuncs];
  for (int i = 0; i<nfuncs; i++) {
    char fname[20];
    char f_spacename[20];
    strcpy(fname, Form("f%d", i));
    strcpy(f_spacename, Form("f_space%d", i));
    f[i] = new TF1(fname, "pol3");
    f_space[i] = new TCutG(f_spacename, 0);
    c.push_back(0.0);
  }
}

/*cut
 *Method for making cuts on the theta_notilt histogram
 *Makes one for each polynomial
 *Also names the Cuts appropriately
 */
void fit::cut() {
 
  cout<<"Draw "<< nfuncs << " fit cuts"<< endl;
  for (int i=0; i<nfuncs; i++) {
    cout<<i+1<<endl;
    x1_theta_notilt->Draw("colz");
    char f_spacename[20];
    strcpy(f_spacename, Form("f_space%d", i));
    while(c1->WaitPrimitive()) {}
    f_space[i] = (TCutG*) c1->GetPrimitive("CUTG");
    f_space[i]->SetName(f_spacename);
  }

}

/*interpolation method
 *Currently uses inverse distance weights to make a weighted average
 *of the polynomial values at a given theta, where the polynomial was recentered
 *at 0. Maybe other methods better? More polynomials == better fit
 *Current method better fit for peaks in between polynomials, not on them
 */
Float_t fit::interp(Float_t x, Float_t theta) {
  Float_t d[nfuncs];
  Float_t weight[nfuncs];
  Float_t denom = 0.0;
  Float_t value = 0.0; 
  for (int i=0; i<nfuncs; i++) {
    Float_t distance = TMath::Abs(f[i]->Eval(theta)-x);
    if (distance == 0.0) {
      return f[i]->Eval(theta)-c[i];
    } else {
      d[i] = 1.0/distance ;
      denom = denom+d[i];
    }
  }
  for (int i = 0; i<nfuncs; i++) {
    weight[i] = d[i]/denom;
    value = value+weight[i]*(f[i]->Eval(theta)-c[i]);
  }
  return value; 
}

/*Untilt
 *Takes data from original raw data file, sorts it based on cuts from histo file
 *makes a linear fit through x|theta and then sends data to flat line (untilted)
 *centered about 0. The new data is filled into a histogram and a tree written to
 *the fitting file
 */
void fit::untilt() {

  cout<<"Draw tilt cut"<<endl;
  h->Draw("colz");
  while(c1->WaitPrimitive()) {}
  tilt_space = (TCutG*) c1->GetPrimitive("CUTG");
  tilt_space->SetName("tilt_space");
 
  vector<Float_t> ft_set;
  vector<Float_t> thetat_set;

  for (int entry = 0; entry < nentries; entry++) {
    if (cutFlag_v[entry]&&tilt_space->IsInside(x1_v[entry], theta_v[entry])){  
      ft_set.push_back(x1_v[entry]);
      thetat_set.push_back(theta_v[entry]);
    }
  }

  TGraph *x1_theta_fit_t = new TGraph(ft_set.size(), &(ft_set[0]), &(thetat_set[0]));
  x1_theta_fit_t->Fit(tilt);

  for (int entry = 0; entry < nentries; entry++) {
    theta_v[entry] = theta_v[entry]-tilt->Eval(x1_v[entry]);
    if (cutFlag_v[entry]) {
      x1_theta_notilt->Fill(x1_v[entry], theta_v[entry]);
      x1_notilt->Fill(x1_v[entry]);
    }
  }
}


/*correct
 *Takes the untilted data and regions from cut to make TGraphs of fit regions
 *Each TGraph is individually fitted with 3rd order poly. Then subtracts 
 *the inerpolated position from the polynomials
 *from the actual position of each datum
 *Corrected data is then filled into histograms and the corrected tree
 */
void fit::correct() {

  vector<Float_t> f_set[nfuncs];//Use vectors since of unknown size
  vector<Float_t> theta_set[nfuncs];
 
  for (int entry = 0; entry < nentries; entry++) {
    for (int j = 0; j<nfuncs; j++) {
      if(cutFlag_v[entry] && f_space[j]->IsInside(x1_v[entry], theta_v[entry])) {
        f_set[j].push_back(x1_v[entry]);
        theta_set[j].push_back(theta_v[entry]);
        break;
      }
    }
  }

  TGraph **x1_theta_fit = new TGraph*[nfuncs];
  for (int i=0; i<nfuncs; i++) {
    x1_theta_fit[i] = new TGraph(f_set[i].size(), &(theta_set[i][0]), &(f_set[i][0]));
    x1_theta_fit[i]->Fit(f[i]);
    c[i] = f[i]->Eval(0.0);
    histoArray->Add(x1_theta_fit[i]);
  }

  for (int entry = 0; entry < nentries; entry++) {
    if (cutFlag_v[entry]) {
      x1_c = x1_v[entry] - interp(x1_v[entry], theta_v[entry]); 
      theta_c = theta_v[entry];
      x1_theta_c->Fill(x1_c, theta_c);
      x1_corrected->Fill(x1_c);
      correctTree->Fill();
    }
  }   
}

/*run
 *Pulls all of the data from the original file and stores in vectors for use in cut,untilt,correct
 *Is where all histos should be created and appended to histoArray
 *Writes all of the corrected data and histograms to a file
 */
void fit::run(char* dataName, char* storageName) {
  TFile *data = new TFile(dataName, "READ");
  TFile *storage = new TFile(storageName, "RECREATE");
  TTree *dataTree = (TTree*) data->Get("SortTree");
  c1 = new TCanvas();

  correctTree = new TTree("correctTree", "correctTree");
  histoArray = new TObjArray();
  h = (TH2F*) data->Get("x1_theta");
  x1_notilt = new TH1F("x1_notilt", "fp1 pos untilted theta", 1000, -300, 300);
  x1_theta_notilt = new TH2F("x1_theta_notilt", "fp1 pos vs theta notilt", 600,-300, 300, 600, -3, 3);
  x1_theta_c = new TH2F("x1_theta_corr", "fp1 pos vs theta corr", 600, -300, 300, 300, -3, 3);
  x1_corrected = new TH1F("x1_corr", "fp1 pos corr", 1000, -300, 300);
 
  histoArray->Add(x1_notilt);
  histoArray->Add(x1_theta_notilt);
  histoArray->Add(x1_theta_c);
  histoArray->Add(x1_corrected);
 
  dataTree->SetBranchAddress("x1", &x1_d);
  dataTree->SetBranchAddress("theta", &theta_d);
  dataTree->SetBranchAddress("cutFlag", &cutFlag_d);

  correctTree->Branch("x1_c", &x1_c, "x1_c/F");
  correctTree->Branch("theta_c", &theta_c, "theta_c/F");

  nentries = dataTree->GetEntries();
  for(int event = 0; event<nentries; event++) {
    dataTree->GetEntry(event);
    x1_v.push_back(x1_d);
    theta_v.push_back(theta_d);
    cutFlag_v.push_back(cutFlag_d);
  }
  
  untilt();
  cut();
  correct();
  c1->Close();
 
  x1_v.clear();
  theta_v.clear();

  histoArray->Write();
  data->Close();
  storage->Close();
}
