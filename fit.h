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
#include "TString.h"
#include <vector>

using namespace std;

class fit {

 public:

  fit(int ipsvar, int iord=2, int n=5); //phase space var, order, n polynomials
                                        //psvar = 1 for theta, 2 for phi, 3 for y -- kgh/20190307

  ~fit();

  void run(TString infileName, TString intreeName, TString outfileName, TString outtreeName);

  void ToggleFixTilt(); //kgh/20190307

  private:
    void untilt();
    void sort();
    void cut();
    void correct();
    Float_t interp(Float_t x, Float_t psvar_val);

    Bool_t GenCuts(Int_t elem);
    
    //sets of data for fitting
    TCutG *tilt_space; 
    TCutG **f_space;   //x|theta
    //fit functions  
    TF1 *tilt;
    TF1 **f;//x|theta

    //----------------- kgh/20190307 -----------------//
    Int_t psvar, ord;
    Int_t numevents;

    TFile *infile, *outfile;
    TTree *intree, *outtree;

    //phase space variables to update with corrections + necessary tree values
    Float_t
      *x_ave_v,
      *theta_v,
      *y_ave_v,
      *phi_v,
      *fp1_tdiff_v,
      *fp2_tdiff_v,
      *fp1_tsum_v,
      *fp2_tsum_v,
      *Anode1_v,
      *Scint1_v;

    //values from intree
    Float_t
      ix_ave,
      itheta,
      iy_ave,
      iphi,
      ifp1_tdiff,
      ifp2_tdiff,
      ifp1_tsum,
      ifp2_tsum,
      iAnode1,
      iScint1;

    //values from outtree
    Float_t
      x_ave,
      theta,
      y_ave,
      phi,
      fp1_tdiff,
      fp2_tdiff,
      fp1_tsum,
      fp2_tsum,
      Anode1,
      Scint1;

    //histograms, cuts, etc. we're going to keep track of with each iteration:
    TH1F
      *x_ave_hist;

    TH2F
      *thvx,
      *phvx,
      *yvx;

    TCanvas *ps_canv; //phase space canvas

    TCutG
      *particleID_cut,
      *tdiff1_v_tdiff2_cut,
      *tsum1_v_tsum2_cut;

    Float_t *psvar_v; //phase space variable we wanna work with

    Bool_t fix_tilt; //toggles whether to fix tilt or not (e.g. if it's already been done)

    //------------------------------------------------//

    Float_t *c; //"centers" of the polynomials
    int nfuncs; //number of polynomials
};

#endif
