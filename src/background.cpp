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

#include "background.h"

Backgnd::Backgnd() {

}

void Backgnd::run(char* inputname, char* outputname) {
  TFile *inputFile = new TFile(inputname, "READ");
  TFile *outputFile = new TFile(outputname, "RECREATE");
  rawHisto = (TH1F*) inputFile->Get("x1_corr"); //Pulls corrected position spectrum
  cleanHisto = new TH1F("x1_corr_nobckgnd","x1_corr_nobckgnd", 1000, -300, 300);

  TSpectrum *s = new TSpectrum(100);
  backgndHisto = s->Background(rawHisto); //If choppy: Background(rawHisto, numIterations) 
  cleanHisto->Add(rawHisto, backgndHisto, 1, -1);

  rawHisto->Write();
  backgndHisto->Write();
  cleanHisto->Write();
  inputFile->Close();
  outputFile->Close();
}
