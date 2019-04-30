/*background.h
 *Class for use in SPS analysis to remove background from a spectrum.
 *Designed to be fed a corrected x-position spectrum and then wipe the 
 *background using the TSpectrum tool from ROOT. If the background isn't very smooth, adjust
 *the number of iterations used by the Background() function (default is 10). Returns the 
 *original histogram, the background histogram, and the original minus the background histogram.
 *
 *Gordon M. -- April 2019
 *
 */

#ifndef BACKGROUND_H
#define BACKGROUND_H

#include "TROOT.h"
#include "TFile.h"
#include "TH1.h"
#include "TF1.h"
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
