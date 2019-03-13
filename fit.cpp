/*  fit.cpp
 *  Class for correcting focal plane aberrations using interpolation of several polynomials 
 *  Corrects for tilt in the focal plane along with the characteristic x|theta aberrations
 *  The base constructor gives 5 polynomials, can override to give as many as needed
 *  G.M. Feb 2019
 *  Modified by kgh Mar 2019 to adapt to own analysis, generalize order and phase space variable
 */

#include <iostream>
#include <cstring>

#include "fit.h"
#include "TCanvas.h"
#include "TMath.h"
#include "TApplication.h"

using namespace std;
//Constructors

fit::fit(int ipsvar, int iord, int n) { //ord is order, n is number of fits
                                        //psvar = 1 for theta, 2 for phi, 3 for y -- kgh/20190307

  psvar = (ipsvar > 0 && ipsvar < 4) ? ipsvar : 1; //defaults to theta
  ord = (iord > 0) ? iord : 3; //defaults to third order (strongest aberration)
  nfuncs = (n > 0) ? n : 3; //defaults to three polynomials

  tilt_space = new TCutG("tilt_space",0);
  tilt = new TF1("tilt","pol1");

  c = new Float_t[nfuncs];

  f = new TF1*[nfuncs];
  f_space = new TCutG*[nfuncs];

  for (int i = 0; i<nfuncs; i++) {
    char fname[20];
    char f_spacename[20];
    strcpy(fname, Form("f%d", i));
    strcpy(f_spacename, Form("f_space%d", i));
    f[i] = new TF1(fname, Form("pol%i",ord));
    f_space[i] = new TCutG(f_spacename, 0);
    c[i] = 0;
  }

  fix_tilt = true;

}

fit::~fit() {

  delete tilt_space; tilt_space=0;
  delete tilt; tilt=0;

  for (int i=0; i<nfuncs; i++) {
    delete f[i];
    delete f_space[i];
    f[i]=0;
    f_space[i]=0;
  }

  delete [] f; f=0;
  delete [] f_space; f_space=0;

  delete [] c;

}

/*cut
 *Method for making cuts on the theta_notilt histogram
 *Makes one for each polynomial
 *Also names the Cuts appropriately
 */
void fit::cut() {

  TCanvas *c1 = new TCanvas();
  c1->ToggleToolBar();
 
  cout<<"Draw "<< nfuncs << " fit cuts"<< endl;

  for (int i=0; i<nfuncs; i++) {
    cout<<i+1<<endl;
    switch(psvar) {
    case 1: thvx->Draw("colz"); break;
    case 2: phvx->Draw("colz"); break;
    case 3: yvx->Draw("colz"); break;
    }
    char f_spacename[20];
    sprintf(f_spacename, "f_space%d", i);
    f_space[i] = (TCutG*)c1->WaitPrimitive("CUTG");
    f_space[i]->SetName(f_spacename);
    f_space[i]->Write();
  }

  c1->Close();

}

/*interpolation method
 *Currently uses inverse distance weights to make a weighted average
 *of the polynomial values at a given theta, where the polynomial was recentered
 *at 0. Maybe other methods better? More polynomials == better fit
 *Current method better fit for peaks in between polynomials, not on them
 */
Float_t fit::interp(Float_t x, Float_t psvar_val) {

  Float_t d[nfuncs];
  Float_t weight[nfuncs];
  Float_t denom = 0.0;
  Float_t value = 0.0;

  for (int i=0; i<nfuncs; i++) {
    Float_t distance = TMath::Abs(f[i]->Eval(psvar_val)-x);
    if (distance == 0.0) {
      return f[i]->Eval(psvar_val)-c[i];
    } else {
      d[i] = 1.0/distance;
      denom = denom+d[i];
    }
  }

  for (int i = 0; i<nfuncs; i++) {
    weight[i] = d[i]/denom;
    value = value+weight[i]*(f[i]->Eval(psvar_val)-c[i]);
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

  TH2F *h; //main histogram from which to pull

  switch (psvar) {
  case 1: h = (TH2F*) infile->Get("thvx"); break;
  case 2: h = (TH2F*) infile->Get("phvx"); break;
  case 3: h = (TH2F*) infile->Get("yvx"); break;
  }

  cout<<"Draw tilt cut"<<endl;
  TCanvas *c1 = new TCanvas();
  c1->ToggleToolBar();
  h->Draw("colz");
  tilt_space = (TCutG*)c1->WaitPrimitive("CUTG");
  tilt_space->SetName("tilt_space");
  c1->Close();

  vector<Float_t> ft_set;
  vector<Float_t> psvar_set;

  for (int i=0; i<numevents; i++) {
    if (GenCuts(i) &&
	tilt_space->IsInside(x_ave_v[i], psvar_v[i])) {

      ft_set.push_back(x_ave_v[i]);
      psvar_set.push_back(psvar_v[i]);

    }
  }

  TGraph *x1_psvar_fit_t = new TGraph(ft_set.size(), &(ft_set[0]), &(psvar_set[0]));
  x1_psvar_fit_t->Fit(tilt);

  thvx->Reset();
  phvx->Reset();
  yvx->Reset();
  x_ave_hist->Reset();

  for (int i=0; i<numevents; i++) {

    psvar_v[i] -= tilt->Eval(x_ave_v[i]);

    switch(psvar) {
    case 1: theta_v[i] = psvar_v[i]; break;
    case 2: phi_v[i] = psvar_v[i]; break;
    case 3: y_ave_v[i] = psvar_v[i]; break;
    }

    if (GenCuts(i)) {

      x_ave_hist->Fill(x_ave_v[i]);

      thvx->Fill(x_ave_v[i], theta_v[i]);
      phvx->Fill(x_ave_v[i], phi_v[i]);
      yvx->Fill(x_ave_v[i], y_ave_v[i]);

    }

  }

}

/*sort
 *Takes untilted data set and takes slices from cut to 
 *make TGraphs of all of the regions for fitting
 *TGraphs are then individually fitted with a 3rd order polynomial
 *The number of polynomials is as usual the number specified in the 
 *constructor
 */
void fit::sort() {

  vector<Float_t> f_set[nfuncs];
  vector<Float_t> psvar_set[nfuncs];

  for (int i=0; i<numevents; i++) {
    for (int j = 0; j<nfuncs; j++) {
      if (f_space[j]->IsInside(x_ave_v[i], psvar_v[i]) && GenCuts(i)) {
        f_set[j].push_back(x_ave_v[i]);
        psvar_set[j].push_back(psvar_v[i]);
      }
    }
  }

  TGraph **x_psvar_fit = new TGraph*[nfuncs];

  for (int i=0; i<nfuncs; i++) {
    x_psvar_fit[i] = new TGraph(f_set[i].size(), &(psvar_set[i][0]), &(f_set[i][0]));
    x_psvar_fit[i]->Fit(f[i]);
    c[i] = f[i]->Eval(0.0);
    x_psvar_fit[i]->Write();
  }
  
}

/*correct
 *Takes the untilted data and now subtracts 
 *the inerpolated position from the polynomials
 *from the actual position of each datum
 *Corrected data is then written to file with a few histograms
 */
void fit::correct() {

  thvx->Reset();
  phvx->Reset();
  yvx->Reset();
  x_ave_hist->Reset();

  for (int i=0; i<numevents; i++) {

    if (GenCuts(i)) {

      x_ave_v[i] -= interp(x_ave_v[i], psvar_v[i]); 

      x_ave_hist->Fill(x_ave_v[i]);

      thvx->Fill(x_ave_v[i], theta_v[i]);
      phvx->Fill(x_ave_v[i], phi_v[i]);
      yvx->Fill(x_ave_v[i], y_ave_v[i]);

      x_ave = x_ave_v[i];
      theta = theta_v[i];
      y_ave = y_ave_v[i];
      phi = phi_v[i];
      fp1_tdiff = fp1_tdiff_v[i];
      fp2_tdiff = fp2_tdiff_v[i];
      fp1_tsum = fp1_tsum_v[i];
      fp2_tsum = fp2_tsum_v[i];
      Anode1 = Anode1_v[i];
      Scint1 = Scint1_v[i];

      outtree->Fill();

    }    

  }

  ps_canv->cd(1);
  thvx->Draw("colz");
  ps_canv->cd(2);
  phvx->Draw("colz");
  ps_canv->cd(3);
  yvx->Draw("colz");

}

void fit::run(TString infileName, TString intreeName, TString outfileName, TString outtreeName) {

  infile = new TFile(infileName, "READ");
  intree = (TTree*) infile->Get(intreeName);
  outfile = new TFile(outfileName, "RECREATE");
  outtree = new TTree(outtreeName, outtreeName);

  x_ave_hist = new TH1F("x_ave_hist","x_ave_hist",1024,-300,300);

  thvx = new TH2F("thvx","thvx",512,-350,350,512,-90,90);
  phvx = new TH2F("phvx","phvx",512,-350,350,512,-5,5);
  yvx = new TH2F("yvx","yvx",512,-350,350,512,-20,20);

  ps_canv = new TCanvas("ps_canv","ps_canv");
  ps_canv->Divide(1,3);

  intree->SetBranchAddress("x_ave", &ix_ave);
  intree->SetBranchAddress("theta", &itheta);
  intree->SetBranchAddress("phi", &iphi);
  intree->SetBranchAddress("y_ave", &iy_ave);
  intree->SetBranchAddress("fp1_tdiff", &ifp1_tdiff);
  intree->SetBranchAddress("fp2_tdiff", &ifp2_tdiff);
  intree->SetBranchAddress("fp1_tsum", &ifp1_tsum);
  intree->SetBranchAddress("fp2_tsum", &ifp2_tsum);
  intree->SetBranchAddress("Anode1", &iAnode1);
  intree->SetBranchAddress("Scint1", &iScint1);

  outtree->Branch("x_ave", &x_ave);
  outtree->Branch("theta", &theta);
  outtree->Branch("phi", &phi);
  outtree->Branch("y_ave", &y_ave);
  outtree->Branch("fp1_tdiff", &fp1_tdiff);
  outtree->Branch("fp2_tdiff", &fp2_tdiff);
  outtree->Branch("fp1_tsum", &fp1_tsum);
  outtree->Branch("fp2_tsum", &fp2_tsum);
  outtree->Branch("Anode1", &Anode1);
  outtree->Branch("Scint1", &Scint1);

  numevents = intree->GetEntries();

  x_ave_v = new Float_t[numevents];
  theta_v = new Float_t[numevents];
  phi_v = new Float_t[numevents];
  y_ave_v = new Float_t[numevents];
  fp1_tdiff_v = new Float_t[numevents];
  fp2_tdiff_v = new Float_t[numevents];
  fp1_tsum_v = new Float_t[numevents];
  fp2_tsum_v = new Float_t[numevents];
  Anode1_v = new Float_t[numevents];
  Scint1_v = new Float_t[numevents];

  psvar_v = new Float_t[numevents];

  for (int i=0; i<numevents; i++) {

    intree->GetEntry(i);

    x_ave_v[i] = ix_ave;
    theta_v[i] = itheta;
    phi_v[i] = iphi;
    y_ave_v[i] = iy_ave;
    fp1_tdiff_v[i] = ifp1_tdiff;
    fp2_tdiff_v[i] = ifp2_tdiff;
    fp1_tsum_v[i] = ifp1_tsum;
    fp2_tsum_v[i] = ifp2_tsum;
    Anode1_v[i] = iAnode1;
    Scint1_v[i] = iScint1;

    switch(psvar) {
    case 1: psvar_v[i] = itheta; break;
    case 2: psvar_v[i] = iphi; break;
    case 3: psvar_v[i] = iy_ave; break;
    }

  }

  particleID_cut = (TCutG*) infile->Get("particleID_cut");
  tdiff1_v_tdiff2_cut = (TCutG*) infile->Get("tdiff1_v_tdiff2_cut");
  tsum1_v_tsum2_cut = (TCutG*) infile->Get("tsum1_v_tsum2_cut");

  if (fix_tilt) {
    cout << "Untilting...\n";
    untilt();
  }
  cout << "Making cuts...\n";
  cut();
  cout << "Sorting...\n";
  sort();
  cout << "Correcting...\n";
  correct();

  outfile->cd();

  cout << "Writing objects...\n";

  particleID_cut->Write();
  tdiff1_v_tdiff2_cut->Write();
  tsum1_v_tsum2_cut->Write();
  thvx->Write();
  phvx->Write();
  yvx->Write();
  x_ave_hist->Write();
  ps_canv->Write();
  outtree->Write();

  infile->Close();
  outfile->Close();

  delete [] x_ave_v;
  delete [] theta_v;
  delete [] y_ave_v;
  delete [] phi_v;
  delete [] fp1_tdiff_v;
  delete [] fp2_tdiff_v;
  delete [] fp1_tsum_v;
  delete [] fp2_tsum_v;
  delete [] Anode1_v;
  delete [] Scint1_v;

  delete [] psvar_v;

  cout << "Done.\n";

}

void fit::ToggleFixTilt() {fix_tilt = !fix_tilt;}

Bool_t fit::GenCuts(Int_t elem) {

  if (tsum1_v_tsum2_cut->IsInside(fp2_tsum_v[elem], fp1_tsum_v[elem]) &&
      tdiff1_v_tdiff2_cut->IsInside(fp2_tdiff_v[elem], fp1_tdiff_v[elem]) && 
      particleID_cut->IsInside(Scint1_v[elem], Anode1_v[elem]))
    return true;

  return false;

}
