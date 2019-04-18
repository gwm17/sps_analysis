#ifndef BACKGROUND_H
#define BACKGROUND_H

#include "TROOT.h"
#include "TFile.h"
#include "TH1.h"
#include "TF1.h"
#include "TCanvas.h"
#include "TSpectrum.h"

using namespace std;

class Backgnd {

  public: 
    Backgnd();
    void run(char* inputname, char* outputname);

  private:
    TH1 *backgndHisto;
    TH1F *cleanHisto;
    TH1F *rawHisto;

};

#endif
