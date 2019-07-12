/*gausFit()
 * Class to make multifitted gaussian peaks in position spectra
 * from the SESPS
 * 
 *
 */

#include "gausFit.h"

using namespace std;

gausFit::gausFit(int numGaus) {
  nGaussians = numGaus;
  params = new Double_t[nGaussians*3];
}

gausFit::~gausFit() {
  delete[] params;
  delete[] new_params;
}

void gausFit::getHisto(char* filename) {
  TFile *file = new TFile(filename, "READ");
  histo = (TH1F*) file->Get("xavg");
}

void gausFit::getRanges() {
  cout<<"Enter in the range for each gaussian to be fitted"<<endl;
  float min_i, max_i;
  TCanvas *c1 = new TCanvas();
  for (int i=0; i<nGaussians; i++) {
    min_i = 0;
    max_i = 0;
    cout<<"For gaussian "<<i<<": "<<endl;
    histo->Draw();
    while(c1->WaitPrimitive()) {}
    cout<<" Min = ";
    cin>>min_i;
    cout<<" Max = ";
    cin>>max_i;
    min.push_back(min_i);
    max.push_back(max_i);
  }
  cout<<"Enter the range for the entire fit: "<<endl;
  histo->Draw();
  while(c1->WaitPrimitive()) {}
  cout<<"Min = ";
  cin>>fullMin;
  cout<<"Max = ";
  cin>>fullMax;
  c1->Close();
}

void gausFit::createGaus() {
  string full_func = "gaus(";
  for (int i=0; i<nGaussians; i++) {
    string name = "g";
    name += to_string(i);
    full_func += to_string(i*3);
    if(i<(nGaussians-1)) {
      full_func += ")+gaus(";
    } else {
      full_func += ")";
    }
    char name_i[name.length()+1];
    strcpy(name_i, name.c_str());
    TF1 *gaus_i = new TF1(name_i,"gaus", min[i], max[i]);
    gaussians.push_back(gaus_i);
  }
  char full_func_name[full_func.length()+1];
  strcpy(full_func_name, full_func.c_str());
  multigaus = new TF1("complete_fit",full_func_name,fullMin,fullMax);
}

void gausFit::fitIndividuals() {
  for (int i=0; i<nGaussians; i++) {
    histo->Fit(gaussians[i], "R0+");
    gaussians[i]->GetParameters(&params[i*3]);
  }
  multigaus->SetParameters(params);
  histo->Fit(multigaus, "R0+");
  chisq = multigaus->GetChisquare();
  ndf = multigaus->GetNDF();
  r_chisq = chisq/((Double_t)ndf);
  cout<<"Chi-Squared value for fit: "<<chisq<<endl;  
  cout<<"Degrees of freedom: "<<ndf<<endl;
  cout<<"Reduced Chi-square value: "<<r_chisq<<endl;
}

bool gausFit::drawFit() {
  TCanvas *c1 = new TCanvas();
  new_params = new Double_t[nGaussians*3];
  multigaus->GetParameters(&new_params[0]);
  histo->Draw();
  multigaus->SetLineColor(kBlue);
  multigaus->Draw("same");
  for(int i = 0; i<nGaussians; i++) {
    TF1 *gaus = gaussians[i];
    gaus->SetParameters(&new_params[i*3]);
    gaus->SetLineColor(kGreen);
    gaus->DrawF1(fullMax, fullMin, "same");
  }
  while(c1->WaitPrimitive()) {}
  string answer;
  cout<<"Do you want to try again?(y/n) ";
  cin>>answer;
  c1->Close();
  if(answer == "y") {
    return false;
  } else if (answer == "n") {
    return true;
  } else {
    cout<<"That wasn't y or n so guess you're going again"<<endl;
    return false;
  }
}

void gausFit::saveResults(string filename) {
  ofstream outfile(filename.c_str());
  if(outfile.is_open()) {
    outfile<<fixed<<showpoint<<setprecision(8);
    outfile<<"Chi-square = "<<chisq<<"\t"<<"DoF = "<<ndf<<"\t"<<"Reduced chi-square = "
           <<r_chisq<<endl;
    outfile<<endl;
    outfile<<setw(10)<<"Gaussian"<<"\t"<<"Amplitude"<<"\t"<<"Centroid"<<"\t"<<"Std. Dev."
           <<"\t"<<"Area"<<endl;
    for(int i=0; i<nGaussians; i++) {
      TF1 *gaus = gaussians[i];
      gaus->SetParameters(&new_params[i*3]);
      Double_t centroid = gaus->GetParameter(1);
      Double_t amplitude = gaus->GetParameter(0);
      Double_t sigma = gaus->GetParameter(2);
      Double_t area = amplitude*sigma*sqrt(2.0*TMath::Pi());
      outfile<<setw(10)<<gaus->GetName()<<"\t"<<amplitude<<"\t"<<centroid<<"\t"<<sigma<<"\t"
             <<area<<endl;
    }
  } else {
    cout<<"Error when writing fit results to file, could not open output file!"<<endl;
  }
}

