/* PeakFit()
 *
 * Class to make complex fits of peaks in spectra
 * Currently can provide a fit over a determined range with a preassigned number of peaks
 * in that range. Will return a reduced chi-square value as an initial test of goodness of fit
 * Current peak shapes allowed are gaussians and breit-wigner distributions
 *
 * Gordon M. -- July 2019
 */

#include "PeakFit.h"
#include "TApplication.h"

using namespace std;

PeakFit::PeakFit() {
  spec = new TSpectrum();
}

PeakFit::~PeakFit() {
  delete[] params;
  delete[] new_params;
}

/* Pulls histogram from a file. If the histogram doesnt exist returns false 
 * and program should get terminated in main. Also gets the bin width from the original histogram
 * to be used with integration results (bin width is assumed constant)
 */
bool PeakFit::getHisto(char* filename, char* histoname) {
  TFile *file = new TFile(filename, "READ");
  if(file->GetListOfKeys()->Contains(histoname)) {
    raw_histo = (TH1F*) file->Get(histoname);
    BIN_WIDTH = raw_histo->GetXaxis()->GetBinWidth(1);
    return true;
  } else {
    cout<<"Error in PeakFit::GetHisto!! Either the file isnt valid or the histogram name "
        <<"does not refer to a valid histogram."<<endl;
    params = new Double_t [1];
    new_params = new Double_t[1];
    return false;
  }
}

/* Uses TSpectrum to remove background from the histogram. More efficient and accurate
 * than trying to provide a polynomial fit of the background
 */
void PeakFit::bckgndRemoval() {
  TCanvas *c1 = new TCanvas();
  cout<<"Remvoing background..."<<endl;
  bool done = false;
  while(!done) {
    Int_t iters;
    //Number of iterations for the background esitmation. 20 is the usual program default.
    //Should note that too many iters ove a small range can wipe the entire spectrum
    cout<<"Enter number of iterations (more equals smoother and slower): ";
    cin>>iters;
    bckgnd = spec->Background(raw_histo, iters);
    histo = new TH1F("clean","clean", 1200, -300, 300);
    histo->GetXaxis()->SetRangeUser(fullMin, fullMax);
    histo->Add(raw_histo, bckgnd, 1, -1);
    histo->Draw();
    while(c1->WaitPrimitive()) {}
    string answer;
    cout<<"Is this satisfactory?(y/n)";
    cin>>answer;
    if(answer == "y") {done = true;}
    else {delete histo;}
  }
}

/* NOT CURRENTLY USED
 * Attempt to have TSpectrum automatically find the peaks in the histogram
 * Unfortunately, it doesnt seem to work well with peaks that are overlapping a lot
 * which kinda defeats the purpose of this whole thing
 */
void PeakFit::findPeaks() {
  Double_t sigma = 0.5;
  nGaussians = spec->Search(histo,sigma,"nobackgroundnoMarkov",0.05);
  Double_t *xpeaks = spec->GetPositionX();
  for (int i=0; i<nGaussians; i++) {
    Float_t min_i =  xpeaks[i]-2*sigma;
    Float_t max_i = xpeaks[i]+2*sigma;
    g_min.push_back(min_i);
    g_max.push_back(max_i);
  }
  params = new Double_t[nGaussians*3];
}

/* NOT CURRENTLY USED
 * If using FindPeaks, no need for user to enter in range of each individual
 */
void PeakFit::getRange() {
  TCanvas *c1 = new TCanvas();
  cout<<"Enter in the range of the histogram to fit over: "<<endl;
  raw_histo->Draw();
  while(c1->WaitPrimitive()) {}
  cout<<"Min = ";
  cin>>fullMin;
  cout<<"Max = ";
  cin>>fullMax;
  c1->Close();
  raw_histo->GetXaxis()->SetRangeUser(fullMin, fullMax);
}

/* Gets the type and range of each peak (and the entire fit) from the user
 * For more complicated functions (not gaus or polN) also need initial guess of params
 */
void PeakFit::getRanges() {
  TCanvas *c1 = new TCanvas();
  cout<<"Enter in the range of the histogram to fit over and the number of peaks: "<<endl;
  raw_histo->Draw();
  while(c1->WaitPrimitive()) {} //this means wait for canvas double click
  cout<<"Min = ";
  cin>>fullMin;
  cout<<"Max = ";
  cin>>fullMax;
  cout<<"Number of peaks: ";
  cin>>nPeaks;
  nGaussians=0; nBW=0;
  for(int i=0; i<nPeaks; i++) {
    Double_t min_i, max_i;
    string answer;
    cout<<"Is this a gaussian or a breit-wigner distribution?(g/b) ";
    cin>>answer;
    if(answer == "g") {
      cout<<"Enter in the range for gaus"<<to_string(nGaussians)<<": "<<endl;
      gaus_id.push_back(i);
      nGaussians++;
      raw_histo->Draw(); 
      while(c1->WaitPrimitive()) {}
      cout<<"Min = ";
      cin>>min_i;
      cout<<"Max = ";
      cin>>max_i;
      g_min.push_back(min_i);
      g_max.push_back(max_i);
    } else if(answer == "b") {
      Double_t mean, width;
      cout<<"Enter in range for breit-wigner"<<to_string(nBW)<<":"<<endl;
      bw_id.push_back(i);
      nBW++;
      raw_histo->Draw(); 
      while(c1->WaitPrimitive()) {}
      cout<<"Min = ";
      cin>>min_i;
      cout<<"Max = ";
      cin>>max_i;
      cout<<"BW requires initial parameter (mean,width) guess"<<endl;
      cout<<"Mean: ";
      cin>>mean;
      cout<<"FWHM: ";
      cin>>width;
      bw_params.push_back(mean);
      bw_params.push_back(width);
      bw_min.push_back(min_i);
      bw_max.push_back(max_i);
    }
  }
  c1->Close();
  raw_histo->GetXaxis()->SetRangeUser(fullMin, fullMax);//restrict range for fit
  totalParams = nGaussians*3+nBW*2; //total Parameters is (num of type) * (params of type)
  params = new Double_t[totalParams];
  //Set number of types for the global fit class
  func.nGaussians = nGaussians;
  func.nBW = nBW;
}

/* Creates all of the TF1's to be used later. Individuals are stored in vectors
 * To make a good estimate for the full function we need both individuals and the 
 * full function
 */
void PeakFit::createFunctions() {
  for (int i=0; i<nGaussians; i++) {
    string name = "g";
    name += to_string(i);
    char name_i[name.length()+1];
    strcpy(name_i, name.c_str());
    TF1 *gaus_i = new TF1(name_i,"gaus", g_min[i], g_max[i]);
    gaussians.push_back(gaus_i);
  }
  for(int i=0; i<nBW; i++) {
    string name = "bw"+to_string(i);
    char n_i[name.length()+1];
    strcpy(n_i, name.c_str());
    TF1 *bw_i = new TF1(n_i,"TMath::BreitWigner(x,[0],[1])",bw_min[i],bw_max[i]);
    bw_i->SetParameter(0, bw_params[i]);
    bw_i->SetParameter(1, bw_params[i+1]);
    breitwigners.push_back(bw_i);
  }
  multigaus = new TF1("complete_fit",func,fullMin,fullMax,totalParams);
}

/* Fits the individual functions and then takes the resulting parameters from the 
 * individual fits and gives them to the full fit as initial parameters. Results in a
 * much stronger fit than if initial guesses are used. (This method is particularly strong for
 * gaussians and other ROOT functions since there is no need for any user guessing)
 */
void PeakFit::fitHisto() {
  for (int i=0; i<nGaussians; i++) {
    histo->Fit(gaussians[i], "R0+");
    gaussians[i]->GetParameters(&params[i*3]);
    multigaus->SetParameter(i*3, params[i*3]);
    multigaus->SetParameter(i*3+1, params[i*3+1]);
    multigaus->SetParameter(i*3+2, params[i*3+2]);
  }
  for (int i=0; i<nBW; i++) {
    histo->Fit(breitwigners[i], "R0+");
    int bwi = nGaussians*3+i*2;
    breitwigners[i]->GetParameters(&params[bwi]);
    multigaus->SetParameter(bwi, params[bwi]);
    multigaus->SetParameter(bwi+1, params[bwi+1]);
  }
  histo->Fit(multigaus, "R0+");
  //Returns a reduced chi square value as an inital estimate of goodness of fit
  chisq = multigaus->GetChisquare();
  ndf = multigaus->GetNDF();
  r_chisq = chisq/((Double_t)ndf);
  cout<<"Chi-Squared value for fit: "<<chisq<<endl;  
  cout<<"Degrees of freedom: "<<ndf<<endl;
  cout<<"Reduced Chi-square value: "<<r_chisq<<endl;
}

/* Draws fit as both the single global function and the individuals with the parameters from
 * the global fit. User then has the option to either indicate desire to try the fit again or
 * to accept the current fit (to be handled by main)
 */
bool PeakFit::drawFit() {
  TCanvas *c1 = new TCanvas();
  new_params = new Double_t[totalParams];
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
  for(int i=0; i<nBW; i++) {
    TF1 *bw = breitwigners[i];
    int bwi = nGaussians*3+i*2;
    bw->SetParameters(&new_params[bwi]);
    bw->SetLineColor(kRed);
    bw->DrawF1(fullMax, fullMin, "same");
  }
  while(c1->WaitPrimitive()) {}
  string answer;
  cout<<"Do you want to try again?(y/n) ";
  cin>>answer;
  c1->Close();
  //indicate try again or not. the true or false value should then be used by main to 
  //either re-run or finish
  if(answer == "y") {
    return false;
  } else if (answer == "n") {
    return true;
  } else {
    cout<<"That wasn't y or n so guess you're going again"<<endl;
    return false;
  }
}

//Stores the results of the fit in a txt file (chi square, params, and integrated area)
void PeakFit::saveResults(char* filename) {
  ofstream outfile(filename);
  if(outfile.is_open()) {
    outfile<<fixed<<showpoint<<setprecision(7);
    outfile<<"Chi-square = "<<chisq<<"\t"<<"DoF = "<<ndf<<"\t"<<"Reduced chi-square = "
           <<r_chisq<<endl;
    outfile<<endl;
    outfile<<setw(10)<<"Gaussian"<<"\t"<<setw(10)<<"Amplitude"<<"\t"<<setw(10)<<"Centroid"
           <<"\t"<<setw(10)<<"Std. Dev."<<"\t"<<"Area"<<endl;
    //Note that integrals are divided by BIN_WIDTH... Numerical integration of the functions
    //is dependent on the bin width while a true integral wouldn't be impacted by the size
    //of dx since it would be assumed infinitesimal
    for(int i=0; i<nGaussians; i++) {
      TF1 *gaus = gaussians[i];
      gaus->SetParameters(&new_params[i*3]);
      Double_t centroid = gaus->GetParameter(1);
      Double_t amplitude = gaus->GetParameter(0);
      Double_t sigma = fabs(gaus->GetParameter(2));
      Double_t area = gaus->Integral(centroid-3*sigma, centroid+3*sigma)/BIN_WIDTH;
      outfile<<setw(10)<<gaus->GetName()<<"\t"<<amplitude<<"\t"<<centroid<<"\t"<<sigma<<"\t"
             <<area<<endl;
    }
    outfile<<endl;
    outfile<<setw(10)<<"BreitWigner"<<"\t"<<setw(10)<<"Mean"<<"\t"<<setw(10)<<"FWHM"<<"\t"
           <<setw(10)<<"Area"<<endl;
    for(int i=0; i<nBW; i++) {
      TF1 *bw = breitwigners[i];
      int bwi = nGaussians*3+i*2;
      bw->SetParameters(&new_params[bwi]);
      Double_t mean = bw->GetParameter(0);
      Double_t width = bw->GetParameter(1);
      Double_t area = bw->Integral(mean-width, mean+width)/BIN_WIDTH;
      outfile<<setw(10)<<bw->GetName()<<"\t"<<mean<<"\t"<<width<<"\t"<<area<<endl;
    }
  } else {
    cout<<"Error when writing fit results to file, could not open output file!"<<endl;
  }
  outfile.close();
}

/* Currently has own main which asks for 2 inputs (infile and outfile) to make it not 
 * restricted to one specific setup, but could be very easily folded into a larger analysis
 * process if desired
 */
int main(int argc, char **argv) {
  if(argc == 3) {
    TApplication *app = new TApplication("app", &argc, argv);
    argc = app->Argc();
    argv = app->Argv();
    cout<<"Input file: "<<argv[1]<<" Output file: "<<argv[2]<<endl;
    bool goodfit = false;
    while(!goodfit) {
      PeakFit pf;
      string answer;
      cout<<"Enter name of histogram to be fitted: ";
      cin>>answer;
      char histoname[answer.length()+1];
      strcpy(histoname, answer.c_str());
      if(pf.getHisto(argv[1], histoname)) {
        pf.getRanges();
        pf.bckgndRemoval();
        pf.createFunctions();
        pf.fitHisto();
        goodfit = pf.drawFit();
        if(goodfit) {
          cout<<"Writing fit results to "<<argv[2]<<endl;
          pf.saveResults(argv[2]);
        }
      } else {
        break;
      }
    }  
  } else {
    cout<<"Incorrect number of arguments! Name of input file, and of output file required!"<<endl;
    cout<<"Terminating abnormally"<<endl;
  }
  return 0;
}
