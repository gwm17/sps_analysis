/*  fit.cpp
 *  Class for correcting focal plane aberrations using interpolation of several polynomials 
 *  Corrects for tilt in the focal plane along with the characteristic x|theta aberrations
 *  The base constructor gives 5 polynomials, can override to give as many as needed
 *  G.M. Feb 2019
 */

#include "fit.h"
#include "TCanvas.h"
#include <iostream>
#include "TMath.h"
#include "TApplication.h"

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
    sprintf(fname, "f%d", i);
    sprintf(f_spacename, "f_space%d", i);
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
    sprintf(fname, "f%d", i);
    sprintf(f_spacename, "f_space%d", i);
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
void fit::cut(char* storageName) {
  TFile *storage = new TFile(storageName, "READ");
  TCanvas *c1 = new TCanvas();
  TH2F *h = (TH2F*) storage->Get("x1_theta_notilt");
 
  cout<<"Draw "<< nfuncs << " fit cuts"<< endl;
  for (int i=0; i<nfuncs; i++) {
    cout<<i+1<<endl;
    h->Draw("colz");
    char f_spacename[20];
    sprintf(f_spacename, "f_space%d", i);
    while(c1->WaitPrimitive()) {}
    f_space[i] = (TCutG*) c1->GetPrimitive("CUTG");
    f_space[i]->SetName(f_spacename);
  }

  storage->Close();
  c1->Close();
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
void fit::untilt(char* dataName,char* storageName, char* histoName) {
  TFile *data = new TFile(dataName, "READ");
  TFile *histo = new TFile(histoName, "READ");
  TFile *storage = new TFile(storageName, "RECREATE");
  TTree *tree = (TTree*) data->Get("DataTree");
  TCutG *s1a1_cut = (TCutG*) histo->Get("s1a1_cut");
  TCutG *x1x2_cut = (TCutG*) histo->Get("x1x2_cut");
  TCutG *fp1anode1_cut = (TCutG*) histo->Get("fp1anode1_cut");
  TH2F *h = (TH2F*) histo->Get("x1_theta");
  TH1F *x1_notilt = new TH1F("x1_notilt", "fp1 pos untilted theta", 1000, -300, 300);

  cout<<"Draw tilt cut"<<endl;
  TCanvas *c1 = new TCanvas();
  h->Draw("colz");
  while(c1->WaitPrimitive()) {}
  tilt_space = (TCutG*) c1->GetPrimitive("CUTG");
  tilt_space->SetName("tilt_space");
  c1->Close();
 
  TH2F *x1_theta_notilt = new TH2F("x1_theta_notilt", "fp1 pos vs theta notilt", 600,-300, 300, 600, -300, 300);
  Float_t tdiff1;
  Float_t tdiff2;
  Int_t Anode1;
  Int_t Scint1;
  Float_t theta_notilt;

  tree->SetBranchAddress("fp_plane1_tdiff", &tdiff1);
  tree->SetBranchAddress("fp_plane2_tdiff", &tdiff2);
  tree->SetBranchAddress("FP1", &Anode1);
  tree->SetBranchAddress("Scint1", &Scint1);
  vector<Float_t> ft_set;
  vector<Float_t> thetat_set;

  for (int i = 0; i < tree->GetEntries(); i++) {
    tree->GetEntry(i);
    tdiff1 = tdiff1*1/1.86;
    tdiff2 = tdiff2*1/1.86;
    Float_t theta = tdiff1-tdiff2;
    if ( fp1anode1_cut->IsInside(tdiff1, Anode1) && x1x2_cut->IsInside(tdiff1, tdiff2) && 
         s1a1_cut->IsInside(Scint1, Anode1) && tilt_space->IsInside(tdiff1, theta)) {
      ft_set.push_back(tdiff1);
      thetat_set.push_back(theta);
    }
  }
  TGraph *x1_theta_fit_t = new TGraph(ft_set.size(), &(ft_set[0]), &(thetat_set[0]));
  x1_theta_fit_t->Fit(tilt);

  TTree *new_tree = new TTree("corr_data", "corr_data");
  new_tree->Branch("theta", &theta_notilt, "theta/F");
  new_tree->Branch("x1", &tdiff1, "x1/F");
  for (int i = 0; i < tree->GetEntries(); i++) {
    tree->GetEntry(i);
    tdiff1 = tdiff1*1/1.86;
    tdiff2 = tdiff2*1/1.86;
    Float_t theta = tdiff1-tdiff2;
    theta_notilt = theta - tilt->Eval(tdiff1);
    if ( fp1anode1_cut->IsInside(tdiff1, Anode1) && x1x2_cut->IsInside(tdiff1, tdiff2) && 
         s1a1_cut->IsInside(Scint1, Anode1)) {
      x1_theta_notilt->Fill(tdiff1, theta_notilt);
      x1_notilt->Fill(tdiff1);
      new_tree->Fill();
    }
  }
  histo->Close();
  data->Close();
  x1_theta_fit_t->Write();
  x1_theta_notilt->Write();
  x1_notilt->Write();
  new_tree->Write();
  storage->Close();
}

/*sort
 *Takes untilted data set and takes slices from cut to 
 *make TGraphs of all of the regions for fitting
 *TGraphs are then individually fitted with a 3rd order polynomial
 *The number of polynomials is as usual the number specified in the 
 *constructor
 */
void fit::sort(char* storageName) {
  TFile *storage = new TFile(storageName, "UPDATE");
  TTree *tree = (TTree*) storage->Get("corr_data");

  Float_t x1;
  Float_t theta;
 
  tree->SetBranchAddress("x1", &x1);
  tree->SetBranchAddress("theta", &theta);
  vector<Float_t> f_set[nfuncs];//Use vectors since of unknown size
  vector<Float_t> theta_set[nfuncs];
 
  for (int i = 0; i < tree->GetEntries(); i++) {
    tree->GetEntry(i);
    for (int j = 0; j<nfuncs; j++) {
      if(f_space[j]->IsInside(x1, theta)) {
        f_set[j].push_back(x1);
        theta_set[j].push_back(theta);
      }
    }
  }
  TGraph **x1_theta_fit = new TGraph*[nfuncs];
  for (int i=0; i<nfuncs; i++) {
    x1_theta_fit[i] = new TGraph(f_set[i].size(), &(theta_set[i][0]), &(f_set[i][0]));
    x1_theta_fit[i]->Fit(f[i]);
    c[i] = f[i]->Eval(0.0);
    x1_theta_fit[i]->Write();
  }
  storage->Close();
}

/*correct
 *Takes the untilted data and now subtracts 
 *the inerpolated position from the polynomials
 *from the actual position of each datum
 *Corrected data is then written to file with a few histograms
 */
void fit::correct(char* storageName) {
  TFile *storage = new TFile(storageName, "UPDATE");
  TTree *tree = (TTree*) storage->Get("corr_data");

  TH2F *x1_theta_corr = new TH2F("x1_theta_corr", "fp1 pos vs theta corr", 600, -300, 300, 300, -150, 150);
  TH1F *x1_corr = new TH1F("x1_corr", "fp1 pos corr", 1000, -300, 300);
  Float_t x1;
  Float_t theta;

  tree->SetBranchAddress("x1", &x1);
  tree->SetBranchAddress("theta", &theta);

  for (int i = 0; i < tree->GetEntries(); i++) {
    tree->GetEntry(i);

    Float_t corr_x1 = x1 - interp(x1, theta); 
    x1_theta_corr->Fill(corr_x1, theta);
    x1_corr->Fill(corr_x1);
    
  }   
  x1_theta_corr->Write();
  x1_corr->Write();
  storage->Close();
}

void fit::run(char* dataName, char* storageName, char* histoName) {
  untilt(dataName, storageName, histoName);
  cut(storageName);
  sort(storageName);
  correct(storageName);
}
