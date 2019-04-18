#include "background.h"

Backgnd::Backgnd() {

}

void Backgnd::run(char* inputname, char* outputname) {
  TFile *inputFile = new TFile(inputname, "READ");
  TFile *outputFile = new TFile(outputname, "RECREATE");
  rawHisto = (TH1F*) inputFile->Get("x1_corr");
  cleanHisto = new TH1F("x1_corr_nobckgnd","x1_corr_nobckgnd", 1000, -300, 300);

  TSpectrum *s = new TSpectrum(100);
  backgndHisto = s->Background(rawHisto); 
  cleanHisto->Add(rawHisto, backgndHisto, 1, -1);

  rawHisto->Write();
  backgndHisto->Write();
  cleanHisto->Write();
  inputFile->Close();
  outputFile->Close();
}
