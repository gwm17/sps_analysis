/*  fit.h
 *  Class for correcting focal plane aberrations using interpolation of several polynomials 
 *  Corrects for tilt in the focal plane along with the characteristic x|theta aberrations
 *  The base constructor gives 5 polynomials, can override to give as many as needed
 *  G.M. Feb 2019
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
    void run(char* dataName, char* fileName, char* histoName);

  private:
    void untilt(char* dataName, char* storageName, char* histoName);
    void sort(char* storageName);
    void cut(char* storageName);
    void correct(char* storageName);
    Float_t interp(Float_t x, Float_t theta);
    
    //sets of data for fitting
    TCutG *tilt_space; 
    TCutG **f_space;   //x|theta
    //fit functions  
    TF1 *tilt;
    TF1 **f;//x|theta
    
    vector<Float_t> c; //"centers" of the polynomials
    int nfuncs; //number of polynomials
};

#endif
