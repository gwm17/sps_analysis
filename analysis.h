/*
analysis.h
ROOT analysis standalone for the splitpole spectrograph
3 steps of sorting, takes a raw data file and makes a file of histograms and cuts 
input arguments: data file name and storage file name

Dec. 2018
Gordon M.
*/

#ifndef ANALYSIS_H
#define ANALYSIS_H

#include "TROOT.h"
#include "TH1.h"
#include "TH2.h"
#include "TObjArray.h"
#include "TCutG.h"
#include "TTree.h"
#include "TFile.h"

using namespace std;

class analysis
{

  public:
    analysis();
    void run(char* dataName, char* storageName);
  
  private:
    void sort_raw(char* dataName, char* storageName);
    void sort_tclean(char* dataName, char* storageName);
    void sort_full(char* dataName, char* storageName);
    int Tsum1Check(Float_t value);
    int Tsum2Check(Float_t value);
    int SiTimeCheck(Int_t value);

    TCutG *s1a1_cut;
    TCutG *x1x2_cut;
    TCutG *fp1anode1_cut;
    TCutG *sicoinc_cut;
    TCutG *theta_cut;
    Int_t minSi;
    Int_t maxSi;
    Int_t max1;
    Int_t min1;
    Int_t max2;
    Int_t min2;

} ;

#endif      
