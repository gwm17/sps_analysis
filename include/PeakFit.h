/* PeakFit()
 *
 * Class to make complex fits of peaks in spectra
 * Currently can provide a fit over a determined range with a preassigned number of peaks
 * in that range. Will return a reduced chi-square value as an initial test of goodness of fit
 * Current peak shapes allowed are gaussians and breit-wigner distributions
 *
 * Gordon M. -- July 2019
 */

#ifndef PEAKFIT_H
#define PEAKFIT_H

#include <TROOT.h>
#include <TH1.h>
#include <TF1.h>
#include <TFile.h>
#include <TTree.h>
#include <TCanvas.h>
#include <TSpectrum.h>
#include <iostream>
#include <string>
#include <vector>
#include <fstream>

using namespace std;

//Class for the full function. Use a class to avoid making program global variables for 
//number of each type of peak
class MyFunc {
  public:
    Int_t nGaussians, nBW;
    MyFunc() {}
    Double_t operator() (Double_t *x, Double_t *p) {
      Double_t value=0;
      for(int i=0; i<nGaussians; i++) {
        value += gausFunc(x, &p[i*3]);
      }
      for(int i=0; i<nBW; i++) {
        value += bwFunc(x, &p[nGaussians*3+i*2]);
      }
      return value;
    }
    Double_t gausFunc(Double_t *x, Double_t *p) {
      Double_t arg = 0;
      if(p[2] != 0) arg = (x[0]-p[1])/p[2];
      return p[0]*TMath::Exp(-0.5*arg*arg);
    }
    Double_t bwFunc(Double_t *x, Double_t *p) {
      Double_t denom = TMath::Pi()*2.0*((x[0]-p[0])*(x[0]-p[0])+p[1]*p[1]/4.0);
      return p[1]/denom;
    }
};

class PeakFit {

  public:
    PeakFit();
    ~PeakFit();
    void createFunctions();
    void getRange();
    void getRanges();
    void bckgndRemoval();
    void findPeaks();
    bool getHisto(char* filename, char* histoname);
    void saveResults(char* filename);
    void fitHisto();
    bool drawFit();
    
  private:
    vector<TF1*> gaussians, breitwigners;
    TF1* multigaus;
    MyFunc func;
    TSpectrum *spec;
    Double_t BIN_WIDTH;
    vector<Float_t> bw_min, bw_max;
    vector<Float_t> g_min, g_max;
    vector<Int_t> gaus_id, bw_id;
    Int_t nPeaks, nBW, nGaussians, totalParams;
    Float_t fullMax, fullMin;
    Double_t *params;
    Double_t *new_params;
    vector<Double_t> bw_params;
    Double_t chisq, r_chisq;
    Int_t ndf;
    TH1F *raw_histo;
    TH1F *histo;
    TH1 *bckgnd;
};


#endif
