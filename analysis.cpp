/*
analysis.cpp
ROOT analysis standalone for the splitpole spectrograph
3 steps of sorting, takes a raw data file and makes a file of histograms and cuts 
input arguments: data file name and storage file name

Includes scattering chamber si detector for coincidence measurements (pre-SABRE)
Comment out relevant sections (labeled) if not used in run (search GLORP)

Dec. 2018
Gordon M.
*/

#include "analysis.h"
#include "TCanvas.h"
#include <iostream>
//#include "TApplication.h"
using namespace std;

//constructor
analysis::analysis() : 
  s1a1_cut(new TCutG("s1a1_cut", 0)),
  x1x2_cut(new TCutG("x1x2_cut", 0)),
  fp1anode1_cut(new TCutG("fp1anode_cut",0)),
  theta_cut(new TCutG("theta_cut",0)),
  max1(100000), min1(-100000), max2(100000),  min2(-100000)
{
}

/* functions for checking if values are in bounds for 1D histos */
int analysis::Tsum1Check(Float_t value) {

  if (value<max1 && value>min1) {return 1;}
  else {return 0;}

}

int analysis::Tsum2Check(Float_t value) {

  if (value<max2 && value>min2) {return 1;}
  else {return 0;}

}

int analysis::SiTimeCheck(Int_t value) {

  if (value<maxSi && value>minSi) {return 1;}
  else {return 0;}

}
/*end of check functions*/

/*sort_raw
 *First sort, takes the data and makes tsum plots
 *Gates are then applied on the sum data
 *Includes Si time cut if there is to be coincidence
 *in the scattering chamber
 */
void analysis::sort_raw(char* dataName, char* storageName) {

  TCanvas *c1 = new TCanvas();
  TFile *data_file = new TFile(dataName, "READ");
  TFile *histo_file = new TFile(storageName, "RECREATE");
  TTree *data_tree = (TTree*) data_file->Get("DataTree");

  TObjArray *histo_array = new TObjArray();
  TH1F *fp1_tsum = new TH1F("fp1_tsum", "fp1 tsum", 8192, 0, 8191);
  TH1F *fp1_tdiff = new TH1F("fp1_tdiff", "fp1 position", 1000, -300, 300);
  TH1F *fp2_tsum = new TH1F("fp2_tsum", "fp2 tsum", 8192, 0, 8191);
  TH1F *fp2_tdiff = new TH1F("fp2_tdiff", "fp2 position", 1000, -300, 300);
  TH2F *fp_tsum = new TH2F("fp_tsum", "focal plane sum times", 500, 0, 8191, 500, 0, 8191);
  TH1F *si_time = new TH1F("si_time", "si timing", 65535, 0, 65535);

  histo_array->Add(fp1_tsum);
  histo_array->Add(fp1_tdiff); 
  histo_array->Add(fp2_tsum);
  histo_array->Add(fp2_tdiff);
  histo_array->Add(si_time);

  Float_t tsum1;
  Float_t tsum2;
  Float_t tdiff1;
  Float_t tdiff2;
  Int_t mtdc_data[32];

  data_tree->SetBranchAddress("fp_plane1_tsum", &tsum1);
  data_tree->SetBranchAddress("fp_plane2_tsum", &tsum2);
  data_tree->SetBranchAddress("fp_plane1_tdiff", &tdiff1);
  data_tree->SetBranchAddress("fp_plane2_tdiff", &tdiff2);
  data_tree->SetBranchAddress("mTDC.Data", &mtdc_data);

  for (int i = 0; i < data_tree->GetEntries(); i++) {
   
     data_tree->GetEntry(i);
     
     tdiff1 = tdiff1*1/1.86;
     tdiff2 = tdiff2*1/1.86;
     
     fp1_tsum->Fill(tsum1);
     fp2_tsum->Fill(tsum2);
     fp_tsum->Fill(tsum1, tsum2);

     fp1_tdiff->Fill(tdiff1);
     fp2_tdiff->Fill(tdiff2);
     
     //Si scattering chamber coincidence GLORP
     for(int i=16; i<32; i++) {
       if (mtdc_data[i] != 0) si_time->Fill(mtdc_data[i]);
     }
     //////////////////////////////////
 
  }
 
  data_file->Close();//curently never adjusts raw data file

//Where cuts are made; WaitPrimitive returns true until a double click on canvas
  fp1_tsum->Draw();
  while(c1->WaitPrimitive()) {}
  cout<<"enter fp1_tsum min: ";
  cin >> min1;
  cout<< "enter fp1_tsum max: ";
  cin >> max1;

  fp2_tsum->Draw();
  while(c1->WaitPrimitive()) {}
  cout<<"enter fp2_tsum min: ";
  cin >> min2;
  cout<< "enter fp2_tsum max: ";
  cin >> max2;
  
  //Si coincidence GLORP
  si_time->Draw();
  while(c1->WaitPrimitive()) {}
  cout<<"enter si time min: ";
  cin >> minSi;
  cout << "enter si time max: ";
  cin >> maxSi;
  ////////////////

  c1->Close();
  histo_array->Write();
  histo_file->Close();
}

/*sort_tclean
 *Takes tsum sorted data and now makes EdE x1_x2 and fp-anode
 *histograms for a final round of cuts
 *If there is Si EdE, include those cuts here
 */
void analysis::sort_tclean(char* dataName, char* storageName) {

  TCanvas *c1 = new TCanvas();
  TFile *data_file = new TFile(dataName, "READ");
  TFile *histo_file = new TFile(storageName, "UPDATE");
  TTree *data_tree = (TTree*) data_file->Get("DataTree");

  TObjArray *histo_array = new TObjArray();
  TH2F *scint1_anode1 = new TH2F("scint1_anode1", "E_dE", 512, 0, 4095, 512, 0, 4095);
  TH2F *fp1_anode1 = new TH2F("fp1_anode","fp1 pos vs anode",500, -300, 300, 512, 0, 4095);
  TH2F *fp2_anode2 = new TH2F("fp2_anode","fp2 pos vs anode",500, -300, 300, 512, 0, 4095);
  TH2F *x1_x2 = new TH2F("x1_x2", "fp1 pos vs fp2 pos", 500,-300,300,500,-300,300);
  TH2F *x1_theta = new TH2F("x1_theta", "fp pos vs theta", 600,-300,300, 600,-300,300);
  
  histo_array->Add(scint1_anode1);
  histo_array->Add(fp1_anode1);
  histo_array->Add(fp2_anode2);
  histo_array->Add(x1_x2);
  histo_array->Add(x1_theta);
  
  Float_t tdiff1;
  Float_t tdiff2;
  Float_t tsum2;
  Float_t tsum1;
  Int_t Scint;
  Int_t Anode1;
  Int_t Anode2;
  Int_t Cathode;
  Int_t mtdc_data[32];

  data_tree->SetBranchAddress("fp_plane1_tsum", &tsum1);
  data_tree->SetBranchAddress("fp_plane2_tsum", &tsum2);
  data_tree->SetBranchAddress("fp_plane1_tdiff", &tdiff1);
  data_tree->SetBranchAddress("fp_plane2_tdiff", &tdiff2);
  data_tree->SetBranchAddress("FP1", &Anode1);
  data_tree->SetBranchAddress("FP2", &Anode2);
  data_tree->SetBranchAddress("Scint1", &Scint);
  data_tree->SetBranchAddress("Cath", &Cathode);
  data_tree->SetBranchAddress("mTDC.Data", &mtdc_data);

  for (int i = 0; i < data_tree->GetEntries(); i++) {
 
     data_tree->GetEntry(i);

     tdiff1 = tdiff1*1/1.86;
     tdiff2 = tdiff2*1/1.86;
     Float_t theta = tdiff1-tdiff2;     

     if (Tsum1Check(tsum1)){
       scint1_anode1->Fill(Scint, Anode1);
       fp1_anode1->Fill(tdiff1, Anode1);
       fp2_anode2->Fill(tdiff2, Anode2);
       x1_x2->Fill(tdiff1, tdiff2);
       x1_theta->Fill(tdiff1, theta);

     }
     if (Tsum2Check(tsum2)) {}
  }

  data_file->Close();

//again cuts, but now use GetPrimitive to retrieve obj CUTG (cast as TCutG)
  scint1_anode1->Draw("colz");
  TCutG *cut;
  while(c1->WaitPrimitive()) {}
  cut = (TCutG*)c1->GetPrimitive("CUTG");
  s1a1_cut = cut;
  s1a1_cut->SetName("s1a1_cut");

  x1_x2->Draw("colz");
  while(c1->WaitPrimitive()) {}
  cut = (TCutG*)c1->GetPrimitive("CUTG");
  x1x2_cut = cut;
  x1x2_cut->SetName("x1x2_cut");

  fp1_anode1->Draw("colz");
  while(c1->WaitPrimitive()) {}
  cut = (TCutG*)c1->GetPrimitive("CUTG");
  fp1anode1_cut = cut;
  fp1anode1_cut->SetName("fp1anode1_cut");


  x1_theta->Draw("colz");
  while(c1->WaitPrimitive()) {}
  cut = (TCutG*)c1->GetPrimitive("CUTG");
  theta_cut = cut;
  theta_cut->SetName("theta_cut");
 
  c1->Close();
  histo_array->Write();
  histo_file->Close();
}

/*sort_full
 *Takes data through the full range of cuts and produces
 *the majority of the histograms 
 */
void analysis::sort_full(char* dataName, char* storageName) {

  TFile *data_file = new TFile(dataName, "READ");
  TFile *histo_file = new TFile(storageName, "UPDATE");
  TTree *data_tree = (TTree*) data_file->Get("DataTree");
  
  TH1F *fp1_tdiff_ts1a1gate = new TH1F("fp1_tdiff_ts1a1gate", "fp1 pos gated s1a1 and tsum", 1000, -300, 300);
  TH1F *fp2_tdiff_ts1a2gate = new TH1F("fp2_tdiff_ts1a2gate", "fp2 pos gated s1a2 and tsum", 1000, -300, 300);
  TH2F *fp1_anode_ts1a1gate = new TH2F("fp1_anode_ts1a1gate", "fp1 pos vs anode gated s1a1 and tsum", 500, -300, 300, 512, 0, 4095);
  TH1F *fp1_tdiff_all = new TH1F("fp1_tdiff_all", "fp1 pos all gates", 1000, -300, 300);
  TH1F *fp1_tdiff_all_closed = new TH1F("fp1_tdiff_all_closed", "fp1 pos all gates & closed slits", 1000, -300, 300);
  TH1F *xavg = new TH1F("xavg", "Avg position", 1000,-300,300);
  TH1F *xdiff = new TH1F("xdiff", "Theta", 1000, -300, 300);
  TH1F *fp1_y = new TH1F("fp1_y", "fp1 y pos", 4096, 0, 4095); 
  TH1F *fp2_y = new TH1F("fp2_y", "fp2 y pos", 4096, 0, 4095);
  TH2F *fp1_tdiffsum = new TH2F("fp1_tdiffsum", "fp1_tdiffsum", 500, -300,300,512,0,8191);
  TH1F *fp1_tdiff_all_sitime = new TH1F("fp1_tdiff_all_sitime", "fp1 pos all w/coinc time", 1000, -300, 300);  
  TH1F *fp1_tdiff_all_sitime_closed = new TH1F("fp1_tdiff_all_sitime_closed", "fp1 pos all w/coinc time & closed slits", 1000, -300, 300);  
  TH1F *phi_hist = new TH1F("phi", "phi", 8192, -4095, 4095);
 
  TObjArray *histo_array = new TObjArray();
  histo_array->Add(fp1_tdiff_ts1a1gate);
  histo_array->Add(fp2_tdiff_ts1a2gate);
  histo_array->Add(fp1_anode_ts1a1gate);
  histo_array->Add(s1a1_cut);
  histo_array->Add(x1x2_cut);
  histo_array->Add(fp1anode1_cut);
  histo_array->Add(theta_cut);
  histo_array->Add(fp1_tdiff_all);
  histo_array->Add(fp1_tdiff_all_closed);
  histo_array->Add(fp1_tdiffsum);
  histo_array->Add(xavg);
  histo_array->Add(xdiff);
  histo_array->Add(fp1_y);
  histo_array->Add(fp2_y);
  histo_array->Add(fp1_tdiff_all_sitime);
  histo_array->Add(fp1_tdiff_all_sitime_closed);
  histo_array->Add(phi_hist);
 
  Float_t tdiff1;
  Float_t tdiff2;
  Float_t tsum2;
  Float_t tsum1;
  Int_t Scint;
  Int_t Anode1;
  Int_t Anode2;
  Int_t Cathode;
  Int_t mtdc_data[32];

  data_tree->SetBranchAddress("fp_plane1_tdiff", &tdiff1);
  data_tree->SetBranchAddress("fp_plane2_tdiff", &tdiff2);
  data_tree->SetBranchAddress("fp_plane1_tsum", &tsum1);
  data_tree->SetBranchAddress("fp_plane2_tsum", &tsum2);
  data_tree->SetBranchAddress("FP1", &Anode1);
  data_tree->SetBranchAddress("FP2", &Anode2);
  data_tree->SetBranchAddress("Scint1", &Scint);
  data_tree->SetBranchAddress("Cath", &Cathode);
  data_tree->SetBranchAddress("mTDC.Data", &mtdc_data);

  for (int i = 0; i < data_tree->GetEntries(); i++) {
    
    data_tree->GetEntry(i);
    
    tdiff1 = tdiff1*1/1.86;
    tdiff2 = tdiff2*1/1.86;
    Float_t x_avg = tdiff1*0.8+tdiff2*0.2;
    Float_t theta = tdiff1-tdiff2;
    Float_t y1 = mtdc_data[5]-mtdc_data[7];

    if (Tsum1Check(tsum1)){
     
      if (s1a1_cut->IsInside(Scint, Anode1)) {

        fp1_tdiff_ts1a1gate->Fill(tdiff1);
        fp1_anode_ts1a1gate->Fill(tdiff1, Anode1);

        if (x1x2_cut->IsInside(tdiff1, tdiff2) && fp1anode1_cut->IsInside(tdiff1, Anode1)){

          fp1_tdiff_all->Fill(tdiff1);
          if (theta_cut->IsInside(tdiff1, theta)) fp1_tdiff_all_closed->Fill(tdiff1);
          fp1_tdiffsum->Fill(tdiff1, tsum1);
          xdiff->Fill(theta);
          xavg->Fill(x_avg);
          fp1_y->Fill(y1);
          phi_hist->Fill(tdiff2-tdiff1);
          
          //Si scattering chamber coincidence GLORP
          for (int i = 16; i<32; i++) {
            if (SiTimeCheck(mtdc_data[i])){
              fp1_tdiff_all_sitime->Fill(tdiff1);
              if(theta_cut->IsInside(tdiff1,theta)) fp1_tdiff_all_sitime_closed->Fill(tdiff1);
              break;
            }
          }
          ///////////////////////////////////////

        }
      }
    }
  }

  data_file->Close();
  histo_array->Write();
  histo_file->Close();
}

/*run
 *runs all three sorts in proper order
 */
void analysis::run(char* dataName, char* storageName) {
  sort_raw(dataName, storageName);
  sort_tclean(dataName, storageName);
  sort_full(dataName, storageName);
}
