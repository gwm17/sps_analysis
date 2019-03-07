/*  fit.h
 *  Class for correcting focal plane aberrations using interpolation of several polynomials 
 *  Corrects for tilt in the focal plane along with the characteristic x|theta aberrations
 *  The base constructor gives 5 polynomials, can override to give as many as needed
 *  G.M. Feb 2019
 *  Modified by kgh Mar 2019 to adapt to own analysis, generalize order and phase space variable
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
  fit(int psvar, int ord=2, int n=5); //phase space var, order, n polynomials
                                      //psvar = 1 for theta, 2 for phi, 3 for y -- kgh/20190307
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
