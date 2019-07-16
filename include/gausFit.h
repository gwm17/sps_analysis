
/*gausFit()
 * Class to make multifitted gaussian peaks in position spectra
 * from the SESPS
 * 
 *
 */

#ifndef GAUSFIT_H
#define GAUSFIT_H

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

class gausFit {

  public:
    gausFit(int numGaus);
    ~gausFit();
    void createGaus();
    void getRanges();
    void getHisto(char* filename);
    void saveResults(string filename);
    void fitIndividuals();
    bool drawFit();
    
  private:
    vector<TF1*> gaussians;
    TF1* multigaus;
    TF1* quad_bckgnd;
    TSpectrum *spec;
    const float BIN_WIDTH = 0.5;
    vector<Float_t> min;
    vector<Float_t> max;
    Int_t nGaussians;
    Float_t fullMax, fullMin;
    Double_t *params;
    Double_t *new_params;
    Double_t chisq, r_chisq;
    Int_t ndf;
    TH1F *histo;
};

#endif
