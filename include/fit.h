/*  fit.h
 *  Class for correcting focal plane aberrations using interpolation of several polynomials 
 *  Corrects for tilt in the focal plane along with the characteristic x|theta aberrations
 *  The base constructor gives 5 polynomials, can override to give as many as needed
 *  G.M. Feb 2019
 *  Revised March 2019 to run without reopen and closing files as shown by KGH -- G.M.
 */

#ifndef FIT_H
#define FIT_H

#include "TROOT.h"
#include "TH1.h"
#include "TH2.h"
#include "TF1.h"
#include "TCutG.h"
#include "TGraph.h"
#include "TFile.h"
#include "TTree.h"
#include "TVectorF.h"
#include <vector>

using namespace std;

class fit 
{

  public:
    fit(); //base constructor; 5 polynomials
    fit(int n); //override; n polynomials
    void run(char* dataName, char* fileName);

  private:
    void untilt();
    void cut();
    void correct();
    Float_t interp(Float_t x, Float_t theta);
    
    //sets of data for fitting
    TCutG *tilt_space; 
    TCutG **f_space;   //x|theta
    //fit functions  
    TF1 *tilt;
    TF1 **f;//x|theta
    
    //vectors to store original data
    vector<Float_t> x1_v;
    vector<Float_t> theta_v;
    vector<Int_t> cutFlag_v;
    
    //data branch variables
    Float_t x1_d,
    theta_d;
    Int_t cutFlag_d;

    //corrected branch variables
    Float_t x1_c,
    theta_c;
    
    //global histogram objects
    TObjArray *histoArray;
    TH2F *h;
    TTree *correctTree;
    int nentries;
    
    TH2F *x1_theta_c;
    TH1F *x1_corrected;
    TH2F *x1_theta_notilt;
    TH1F *x1_notilt;

    vector<Float_t> c; //"centers" of the polynomials
    int nfuncs; //number of polynomials
    TCanvas *c1;
};

#endif
