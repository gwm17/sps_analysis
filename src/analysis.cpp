/*
analysis.cpp
ROOT analysis standalone for the splitpole spectrograph
3 steps of sorting, takes a raw data file and makes a file of histograms and cuts 
input arguments: data file name and storage file name

Includes scattering chamber si detector for coincidence measurements (pre-SABRE)
Comment out relevant sections (labeled) if not used in run (search GLORP)

Dec. 2018
Gordon M.

Revised March 2019 to run without reopening and reclosing files as shown by KGH -- Gordon M.
*/

#include "analysis.h"
#include "TCanvas.h"
#include "FP_kinematics.h"
#include <iostream>
//#include "TApplication.h"
using namespace std;

//constructor
analysis::analysis() : 
  s1a1_cut(new TCutG("s1a1_cut", 0)),
  x1x2_cut(new TCutG("x1x2_cut", 0)),
  fp1anode1_cut(new TCutG("fp1anode_cut",0)),
  theta_cut(new TCutG("theta_cut",0)),
  fp1plast_cut(new TCutG("fp1plast_cut",0)),
  max1(100000), min1(-100000), max2(100000),  min2(-100000)
{
}
analysis::~analysis() {
  delete mtdc_d;
}
/*dump empties*/
int analysis::notEmpty(Int_t value) {
  if (value>1.0) {return 1;}
  else {return 0;}
}
/* functions for checking if values are in bounds for 1D histos */
int analysis::TCheck1Check(Float_t value) {

  if (value<max1 && value>min1) {return 1;}
  else {return 0;}

}

int analysis::TCheck2Check(Float_t value) {

  if (value<max2 && value>min2) {return 1;}
  else {return 0;}

}

int analysis::SiTimeCheck(Int_t value) {

  if (value<maxSi && value>minSi) {return 1;}
  else {return 0;}

}

void analysis::Reset() {
  tdiff1_n = -1e6;
  tdiff2_n = -1e6;
  tsum1_n = -1e6;
  tsum2_n = -1e6;
  tcheck2_n = -1e6;
  tcheck1_n = -1e6;
  theta_n = -1e6;
  phi_n = -1e6;
  y1_n = -1e6;
  y2_n = -1e6;
  
  anode1_n = -1,
  anode2_n = -1,
  cutFlag_n = -1,
  scint1_n = -1;
}
/*end of check functions*/

void analysis::GetWeights() {
  int Zt = 6, At = 12, Zp = 1, Ap = 2, Ze = 1, Ae = 1;
  double Ep = 16.0, angle = 20.0, B = 8840;
  w1 = (Wire_Dist()/2.0-Delta_Z(Zt,At,Zp,Ap,Ze,Ae,Ep,angle,B))/Wire_Dist();
  w2 = 1.0-w1;
}

/*sort_raw
 *First sort, takes the data and makes tsum plots
 *Gates are then applied on the sum data
 *Includes Si time cut if there is to be coincidence
 *in the scattering chamber
 */
void analysis::sort_raw() {

  TCanvas *c1 = new TCanvas();
  for (int entry = 0; entry < nentries; entry++) {
     vector<Int_t> mtdc = mtdc_v[entry];
     if(notEmpty(mtdc[1]) && notEmpty(mtdc[2])){
       Float_t tdiff1 = tdiff1_v[entry]*1/1.83;
       Float_t tcheck1 = tsum1_v[entry]/2.0-anode1_time_v[entry]*0.0625;
       fp1_tsum->Fill(tsum1_v[entry]);
       fp1_tdiff->Fill(tdiff1);
       fp1_tcheck->Fill(tcheck1);
     }
     if(notEmpty(mtdc[3]) && notEmpty(mtdc[4])){
       Float_t tdiff2 = tdiff2_v[entry]*1/1.969;
       Float_t tcheck2 = tsum2_v[entry]/2.0-anode2_time_v[entry]*0.0625;
       fp2_tsum->Fill(tsum2_v[entry]);
       fp2_tdiff->Fill(tdiff2);
       fp2_tcheck->Fill(tcheck2);
     }
     
     //Si scattering chamber coincidence GLORP
     /*for(int i=16; i<32; i++) {
       if (mtdc[i] != 0) si_time->Fill(mtdc[i]);
     }*/
     //////////////////////////////////
 
  }
 

//Where cuts are made; WaitPrimitive returns true until a double click on canvas
  fp1_tcheck->Draw();
  while(c1->WaitPrimitive()) {}
  cout<<"enter fp1_tcheck min: ";
  cin >> min1;
  cout<< "enter fp1_tcheck max: ";
  cin >> max1;

  fp2_tcheck->Draw();
  while(c1->WaitPrimitive()) {}
  cout<<"enter fp2_tcheck min: ";
  cin >> min2;
  cout<< "enter fp2_tcheck max: ";
  cin >> max2;
  
  //Si coincidence GLORP
  /*si_time->Draw();
  while(c1->WaitPrimitive()) {}
  cout<<"enter si time min: ";
  cin >> minSi;
  cout << "enter si time max: ";
  cin >> maxSi;*/
  ////////////////

  c1->Close();
}

/*sort_tclean
 *Takes tsum sorted data and now makes EdE x1_x2 and fp-anode
 *histograms for a final round of cuts
 *If there is Si EdE, include those cuts here
 */
void analysis::sort_tclean() {

  TCanvas *c1 = new TCanvas();
  for (int entry = 0; entry <nentries; entry++) {
    vector<Int_t> mtdc = mtdc_v[entry];
    if (notEmpty(mtdc[1]) && notEmpty(mtdc[2]) && notEmpty(mtdc[3]) && notEmpty(mtdc[4])) {
      Float_t tdiff1 = tdiff1_v[entry]*1/1.83;
      Float_t tdiff2 = tdiff2_v[entry]*1/1.969;
      Float_t tcheck1 = tsum1_v[entry]/2.0-anode1_time_v[entry]*0.0625;
      Float_t tcheck2 = tsum2_v[entry]/2.0-anode2_time_v[entry]*0.0625;
      Float_t theta = (tdiff2-tdiff1)/36.0; //36 mm separation between wires     

      if (TCheck1Check(tcheck1)){
        if(notEmpty(anode1_v[entry])) { 
          scint1_anode1->Fill(scint1_v[entry], anode1_v[entry]);
          fp1_anode1->Fill(tdiff1, anode1_v[entry]);
        }
        x1_x2->Fill(tdiff1, tdiff2);
        x1_theta->Fill(tdiff1, theta);
      }
      if (TCheck2Check(tcheck2)) {
        fp2_anode2->Fill(tdiff2, anode2_v[entry]);
      }
    }
  }


//again cuts, but now use GetPrimitive to retrieve obj CUTG (cast as TCutG)
  /*scint1_anode1->Draw("colz");
  TCutG *cut;
  while(c1->WaitPrimitive()) {}
  cut = (TCutG*)c1->GetPrimitive("CUTG");
  s1a1_cut = cut;
  s1a1_cut->SetName("s1a1_cut");
  s1a1_cut->SetVarX("scint");
  s1a1_cut->SetVarY("anode1");
  histoArray->Add(s1a1_cut);*/

  x1_x2->Draw("colz");
  while(c1->WaitPrimitive()) {}
  x1x2_cut = (TCutG*)c1->GetPrimitive("CUTG");
  x1x2_cut->SetName("x1x2_cut");
  x1x2_cut->SetVarX("x1");
  x1x2_cut->SetVarY("x2");
  histoArray->Add(x1x2_cut);

  fp1_anode1->Draw("colz");
  while(c1->WaitPrimitive()) {}
  fp1anode1_cut = (TCutG*)c1->GetPrimitive("CUTG");
  fp1anode1_cut->SetName("fp1anode1_cut");
  fp1anode1_cut->SetVarX("x1");
  fp1anode1_cut->SetVarY("anode1");
  histoArray->Add(fp1anode1_cut);


  /*x1_theta->Draw("colz");
  while(c1->WaitPrimitive()) {}
  theta_cut = (TCutG*)c1->GetPrimitive("CUTG");
  theta_cut->SetName("theta_cut");
  theta_cut->SetVarX("x1");
  theta_cut->SetVarY("theta");
  histoArray->Add(theta_cut);*/

  for(int i=0; i<nentries; i++) {
    vector<Int_t> mtdc = mtdc_v[i];
    Float_t tdiff1 = tdiff1_v[i]*1/1.86;
    Float_t anode1 = anode1_v[i];
    if(fp1anode1_cut->IsInside(tdiff1, anode1)) {
      Float_t scint1_time = scint1_time_v[i]*0.0625;
      fp1_plastic_time->Fill(tdiff1, scint1_time);
      Float_t rf_scint_time_wrapped = fmod(mtdc[9]*0.0625-scint1_time,164.95);
      fp1_rf_scint_wrapped->Fill(tdiff1, rf_scint_time_wrapped);
    }
  }
  
  fp1_plastic_time->Draw("colz");
  while(c1->WaitPrimitive()) {}
  fp1plast_cut = (TCutG*)c1->GetPrimitive("CUTG");
  fp1plast_cut->SetName("fp1plast_cut");
  fp1plast_cut->SetVarX("x1");
  fp1plast_cut->SetVarY("scint_time");
  histoArray->Add(fp1plast_cut);
  
  fp1_rf_scint_wrapped->Draw("colz");
  while(c1->WaitPrimitive()) {}
  fp1rfwrap_cut = (TCutG*)c1->GetPrimitive("CUTG");
  fp1rfwrap_cut->SetName("fp1rfwrap_cut");
  fp1rfwrap_cut->SetVarX("x1");
  fp1rfwrap_cut->SetVarY("rf_scint_wrapped");
  histoArray->Add(fp1rfwrap_cut);
  c1->Close();
}

/*sort_full
 *Takes data through the full range of cuts and produces
 *the majority of the histograms 
 */
void analysis::sort_full() {

  GetWeights();
  for (int entry = 0; entry < nentries; entry++) {
    vector<Int_t> mtdc = mtdc_v[entry];
    cutFlag_n = 0;
    coincFlag_n = 0;
    if (notEmpty(mtdc[1]) && notEmpty(mtdc[2]) && notEmpty(mtdc[3]) && notEmpty(mtdc[4])) {
      tdiff1_n = tdiff1_v[entry]*1/1.83;
      tdiff2_n = tdiff2_v[entry]*1/1.969;
      tcheck1_n = tsum1_v[entry]/2.0-anode1_time_v[entry]*0.0625;
      tcheck2_n = tsum2_v[entry]/2.0-anode2_time_v[entry]*0.0625;
      tsum1_n = tsum1_v[entry];
      tsum2_n = tsum2_v[entry];
      x_avg_n = tdiff1_n*w1+tdiff2_n*w2;
      theta_n = (tdiff2_n-tdiff1_n)/36.0;
      y1_n = anode1_time_v[entry]-scint1_time_v[entry];
      y2_n = anode2_time_v[entry]-scint1_time_v[entry];
      scint1_time_n = (Float_t)scint1_time_v[entry]*0.0625;
      rf_scint_wrapped_n = fmod(mtdc[9]*0.0625-scint1_time_n, 164.95);
      scint1_n = scint1_v[entry];
      anode1_n = anode1_v[entry];
      anode2_n = anode2_v[entry];
      phi_n = (y2_n-y1_n)/36.0;
      if (TCheck1Check(tcheck1_n) && TCheck2Check(tcheck2_n)){
     
        if (//s1a1_cut->IsInside(scint1_n, anode1_n) && 
            fp1plast_cut->IsInside(tdiff1_n, scint1_time_n)) {

          fp1_tdiff_ts1a1gate->Fill(tdiff1_n);
          fp1_anode_ts1a1gate->Fill(tdiff1_n, anode1_n);

          if(x1x2_cut->IsInside(tdiff1_n,tdiff2_n) && fp1anode1_cut->IsInside(tdiff1_n,anode1_n)){

            fp1_tdiff_all->Fill(tdiff1_n);
            //if (theta_cut->IsInside(tdiff1_n, theta_n)) fp1_tdiff_all_closed->Fill(tdiff1_n);
            fp1_tdiffsum->Fill(tdiff1_n, tsum1_n);
            xdiff->Fill(theta_n);
            xavg->Fill(x_avg_n);
            fp1_y->Fill(y1_n);
            phi_hist->Fill(phi_n);
            cutFlag_n = 1;         
          //Si scattering chamber coincidence GLORP
           /* for (int i = 16; i<32; i++) {
              if (SiTimeCheck(mtdc[i])){
                fp1_tdiff_all_sitime->Fill(tdiff1_n);
                coincFlag_n = 1;
                if(theta_cut->IsInside(tdiff1_n,theta_n)) 
                  fp1_tdiff_all_sitime_closed->Fill(tdiff1_n);
                break;
              }
            }*/
          ///////////////////////////////////////

          }
        }
      }
      sortTree->Fill();
    }
  }
}

/*run
 *runs all three sorts in proper order
 */
void analysis::run(char* dataName, char* storageName) {
  TFile *data = new TFile(dataName, "READ");
  TFile *storage = new TFile(storageName, "RECREATE");
  TTree *dataTree = (TTree*) data->Get("DataTree");
  sortTree = new TTree("SortTree", "SortTree");
  histoArray = new TObjArray();

  fp1_tsum = new TH1F("fp1_tsum", "fp1 tsum", 8192, 0, 8191);
  fp1_tdiff = new TH1F("fp1_tdiff", "fp1 position", 1200, -300, 300);
  fp2_tsum = new TH1F("fp2_tsum", "fp2 tsum", 8192, 0, 8191);
  fp2_tdiff = new TH1F("fp2_tdiff", "fp2 position", 1200, -300, 300);
  si_time = new TH1F("si_time", "si timing", 65535, 0, 65535);
  scint1_anode1 = new TH2F("scint1_anode1", "E_dE", 512, 0, 4095, 512, 0, 4095);
  fp1_anode1 = new TH2F("fp1_anode","fp1 pos vs anode",600, -300, 300, 512, 0, 4095);
  fp2_anode2 = new TH2F("fp2_anode","fp2 pos vs anode",600, -300, 300, 512, 0, 4095);
  x1_x2 = new TH2F("x1_x2", "fp1 pos vs fp2 pos", 600,-300,300,600,-300,300);
  x1_theta = new TH2F("x1_theta", "fp pos vs theta", 600,-300,300, 600,-3,3);
  fp1_tdiff_ts1a1gate = new TH1F("fp1_tdiff_ts1a1gate", "fp1 pos gated s1a1 and tsum", 1200, -300, 300);
  fp1_anode_ts1a1gate = new TH2F("fp1_anode_ts1a1gate", "fp1 pos vs anode gated s1a1 and tsum", 600, -300, 300, 512, 0, 4095);
  fp1_tdiff_all = new TH1F("fp1_tdiff_all", "fp1 pos all gates", 1200, -300, 300);
  fp1_tdiff_all_closed = new TH1F("fp1_tdiff_all_closed", "fp1 pos all gates & closed slits", 1200, -300, 300);
  xavg = new TH1F("xavg", "Avg position", 1200,-300,300);
  xdiff = new TH1F("xdiff", "Theta", 1200, -300, 300);
  fp1_y = new TH1F("fp1_y", "fp1 y pos", 8192, -4095, 4096); 
  fp2_y = new TH1F("fp2_y", "fp2 y pos", 8192, -4095, 4096);
  fp1_tdiffsum = new TH2F("fp1_tdiffsum", "fp1_tdiffsum", 600, -300,300,512,0,8191);
  fp1_tdiff_all_sitime = new TH1F("fp1_tdiff_all_sitime", "fp1 pos all w/coinc time", 1200, -300, 300);  
  fp1_tdiff_all_sitime_closed = new TH1F("fp1_tdiff_all_sitime_closed", "fp1 pos all w/coinc time & closed slits", 1200, -300, 300);  
  phi_hist = new TH1F("phi", "phi", 8192, -4095, 4096);
  fp1_plastic_time = new TH2F("fp1_plastic_time","fp1_plastic_time",600,-300,300,600,0,8191);
  fp1_rf_scint_wrapped = new TH2F("fp1_rf_scint_wrapped","fp1_rf_scint_wrapped",600,-300,300,600,0,8191);
  fp1_tcheck = new TH1F("fp1_tcheck", "fp1_tcheck",8192,0,8191);
  fp2_tcheck = new TH1F("fp2_tcheck", "fp2_tcheck",8192,0,8191);
 
  histoArray->Add(fp1_tdiff_ts1a1gate);
  histoArray->Add(fp1_anode_ts1a1gate);
  histoArray->Add(fp1_tdiff_all);
  histoArray->Add(fp1_tdiff_all_closed);
  histoArray->Add(fp1_tdiffsum);
  histoArray->Add(xavg);
  histoArray->Add(xdiff);
  histoArray->Add(fp1_y);
  histoArray->Add(fp2_y);
  histoArray->Add(fp1_tdiff_all_sitime);
  histoArray->Add(fp1_tdiff_all_sitime_closed);
  histoArray->Add(phi_hist);
  histoArray->Add(scint1_anode1);
  histoArray->Add(fp1_anode1);
  histoArray->Add(fp2_anode2);
  histoArray->Add(x1_x2);
  histoArray->Add(x1_theta);
  histoArray->Add(fp1_tsum);
  histoArray->Add(fp1_tdiff); 
  histoArray->Add(fp2_tsum);
  histoArray->Add(fp2_tdiff);
  histoArray->Add(si_time);
  histoArray->Add(fp1_plastic_time);
  histoArray->Add(fp1_rf_scint_wrapped);

  dataTree->SetBranchAddress("anode1", &anode1_d);
  dataTree->SetBranchAddress("anode2", &anode2_d);
  dataTree->SetBranchAddress("scint1", &scint1_d);
  dataTree->SetBranchAddress("scint2", &scint2_d);
  dataTree->SetBranchAddress("fp_plane1_tdiff", &tdiff1_d);
  dataTree->SetBranchAddress("fp_plane2_tdiff", &tdiff2_d);
  dataTree->SetBranchAddress("fp_plane1_tsum", &tsum1_d);
  dataTree->SetBranchAddress("fp_plane2_tsum", &tsum2_d);
  dataTree->SetBranchAddress("mtdc1", &mtdc_d);
  dataTree->SetBranchAddress("anode1_time", &anode1_time_d);
  dataTree->SetBranchAddress("anode2_time", &anode1_time_d);
  dataTree->SetBranchAddress("plastic_time", &scint1_time_d);

  sortTree->Branch("x1", &tdiff1_n, "x1/F");
  sortTree->Branch("x2", &tdiff2_n, "x2/F");
  sortTree->Branch("tsum1", &tsum1_n, "tsum1/F");
  sortTree->Branch("tsum2", &tsum2_n, "tsum2/F");
  sortTree->Branch("tcheck2", &tcheck2_n, "tcheck2/F");
  sortTree->Branch("tcheck1", &tcheck1_n, "tcheck1/F");
  sortTree->Branch("theta", &theta_n, "theta/F");
  sortTree->Branch("phi", &phi_n, "phi/F");
  sortTree->Branch("y1", &y1_n, "y1/F");
  sortTree->Branch("y2", &y2_n, "y2/F");
  sortTree->Branch("anode1", &anode1_n, "anode1/I");
  sortTree->Branch("anode2", &anode2_n, "anode2/I");
  sortTree->Branch("scint", &scint1_n, "scint/I");
  sortTree->Branch("rf_scint_wrapped", &rf_scint_wrapped_n,"rf_scint_wrapped/F");
  sortTree->Branch("scint_time", &scint1_time_n,"scint_time/F");
  sortTree->Branch("cutFlag", &cutFlag_n, "cutFlag/I");
  sortTree->Branch("coincFlag", &coincFlag_n, "coincFlag/I");

  nentries = dataTree->GetEntries();
  cout<<"entries: "<<nentries<<endl;
  mtdc_v.resize(nentries);
  for (int entry = 0; entry<nentries; entry++) {
    dataTree->GetEntry(entry);
    mtdc_v[entry].resize(32);
    anode1_v.push_back(anode1_d);
    anode2_v.push_back(anode2_d);
    scint2_v.push_back(scint2_d);
    scint1_v.push_back(scint1_d);
    tdiff1_v.push_back(tdiff1_d);
    tdiff2_v.push_back(tdiff2_d);
    tsum1_v.push_back(tsum1_d);
    tsum2_v.push_back(tsum2_d);
    scint1_time_v.push_back(scint1_time_d);
    anode1_time_v.push_back(anode1_time_d);
    anode2_time_v.push_back(anode2_time_d);
    for (int i=0; i<32; i++) {
      mtdc_v[entry][i] = (*mtdc_d)[i];
    }
  }
  
  sort_raw();
  sort_tclean();
  sort_full();

  anode1_v.clear();
  anode2_v.clear();
  scint2_v.clear();
  scint1_v.clear();
  tdiff1_v.clear();
  tdiff2_v.clear();
  tsum1_v.clear();
  tsum2_v.clear();
  mtdc_v.clear();
  scint1_time_v.clear();
  anode1_time_v.clear();
  anode2_time_v.clear();

  sortTree->Write(sortTree->GetName(), TObject::kOverwrite);
  histoArray->Write();
  data->Close();
  storage->Close();
}
