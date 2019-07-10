/*
analysis.h
ROOT analysis standalone for the splitpole spectrograph
3 steps of sorting, takes a raw data file and makes a file of histograms and cuts 
input arguments: data file name and storage file name

Dec. 2018
Gordon M.

Revised March 2019 to run without reopening and reclosing files as shown by KGH -- Gordon M.
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
    ~analysis();
    void run(char* dataName, char* storageName);
  
  private:
    /*functions*/
    void sort_raw();
    void sort_tclean();
    void sort_full();
    int notEmpty(Int_t value);
    int TCheck1Check(Float_t value);
    int TCheck2Check(Float_t value);
    int SiTimeCheck(Int_t value);
    void Reset();
    void GetWeights();

    /*Tree for storing final paramters*/    
    TTree *sortTree; 

    /*storage vectors, so raw file only needs opened once*/
    vector<Int_t> anode1_v,
    anode2_v,
    scint1_v,
    scint2_v;
    
    vector<Float_t> tdiff1_v,
    tdiff2_v,
    tsum1_v,
    tsum2_v,
    scint1_time_v,
    anode1_time_v,
    anode2_time_v;

    vector<vector<Int_t>> mtdc_v;
    
    Float_t w1, w2;

    /*raw branch variables*/
    Int_t anode1_d,
    anode2_d,
    scint1_d,
    scint2_d;

    Float_t tsum1_d,
    tsum2_d,
    tdiff2_d,
    tdiff1_d,
    scint1_time_d,
    anode1_time_d,
    anode2_time_d;

    vector<Int_t> *mtdc_d;
    
    /*new tree variables*/
    Float_t tdiff1_n,
    tdiff2_n,
    tsum1_n,
    tsum2_n,
    tcheck2_n,
    tcheck1_n,
    theta_n,
    phi_n,
    y1_n,
    scint1_time_n,
    rf_scint_wrapped_n,
    x_avg_n,
    y2_n;
    
    Int_t anode1_n,
    anode2_n,
    cutFlag_n,
    coincFlag_n,
    scint1_n;

    /*keep track of number of entries for loops*/
    int nentries;
 
    /*histograms, cuts, etc.*/
    TH1F *fp1_tsum,
    *fp1_tdiff,
    *fp2_tsum,
    *fp2_tdiff,
    *si_time,
    *fp1_tdiff_ts1a1gate,
    *fp1_tdiff_all,
    *fp1_tdiff_all_closed,
    *xavg,
    *xdiff,
    *fp1_y,
    *fp2_y,
    *yavg,
    *fp1_tdiff_all_sitime,
    *fp1_tdiff_all_sitime_closed,
    *fp1_tcheck,
    *fp2_tcheck,
    *phi_hist;

    TH2F *scint1_anode1,
    *fp1_anode1,
    *fp2_anode2,
    *x1_x2,
    *x1_theta,
    *fp1_anode_ts1a1gate,
    *fp1_plastic_time,
    *fp1_rf_scint_wrapped,
    *fp1_tdiffsum;

    TObjArray *histoArray;

    TCutG *s1a1_cut;
    TCutG *x1x2_cut;
    TCutG *fp1anode1_cut;
    TCutG *sicoinc_cut;
    TCutG *theta_cut;
    TCutG *fp1plast_cut;
    TCutG *fp1rfwrap_cut;
    Int_t minSi;
    Int_t maxSi;
    Int_t max1;
    Int_t min1;
    Int_t max2;
    Int_t min2;

} ;

#endif      
