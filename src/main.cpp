/*main.cpp
 *Main function for sps analysis program
 *4 modes: -r run everything, -a only standard analysis, -f only aberration corrections, -b only background removal
 *Takes mode flag and then the data name (data file name w/o .root)
 *data name should be 20 characters or less
 *
 * Gordon M. Feb 2019
 * Updated by G.M. April 2019 for background removal
 */

#include "analysis.h"
#include "fit.h"
#include "background.h"
#include "TROOT.h"
#include "TApplication.h"
#include <iostream>
#include <unistd.h>
#include <string>

using namespace std;

//argument flags
struct options {
  int onlyFit; // -f
  int onlyAnalyze; // -a
  int runAll; // -r
  int cleanBackground; // -b
} options;

//flag string for getopt; if expecting value with flag use : after flag letter
static const char *optString = "farb";

int main(int argc, char* argv[]) {
  int opt = 0;
  options.onlyFit = 0;
  options.onlyAnalyze = 0;
  options.runAll = 0;
  
  char data[strlen(argv[2])+5]; //data name plus five for .root
  char histo[strlen(argv[2])+11]; //plus 11 for _histo.root
  char corr[strlen(argv[2])+10]; //plus 10 for _corr.root
  char clean[strlen(argv[2])+11]; //plus 11 for _clean.root

  strcpy(data, Form("%s.root", argv[2]));
  strcpy(histo, Form("%s_histo.root", argv[2]));
  strcpy(corr, Form("%s_corr.root", argv[2]));
  strcpy(clean, Form("%s_clean.root", argv[2]));

  char *pdata = data; char *phisto = histo; char *pcorr = corr; char *pclean = clean;
 
  opt =  getopt(argc, argv, optString); // 1 = found arg, -1 = no more valid args
  while( opt != -1) {
    switch(opt) {
      case 'f':
        options.onlyFit = 1;
        break;
      case 'a':
        options.onlyAnalyze = 1;
        break;
      case 'r':
        options.runAll = 1;
        break;
      case 'b':
        options.cleanBackground = 1;
        break;
    }
    opt =  getopt(argc, argv, optString); // iterate to next arg
  }
 
  TApplication app("app", &argc, argv);
  if ((options.runAll || options.onlyAnalyze)) {
    cout<<"Running SPS analysis..."<<endl;
    cout<<"Data: "<<pdata<<" Histograms: "<<phisto<<endl;
    cout<<"Sorting data..."<<endl;
    analysis a;
    a.run(pdata, phisto);
    cout<<"Sorting complete."<<endl;
  } if (options.runAll || options.onlyFit) {
    int nfuncs;
    cout<<"Running aberration corrections..."<<endl;
    cout<<"Data: "<<pdata<<" Histograms: "<<phisto<<" Corrections: "<<pcorr<<endl;
    cout<<"Enter number of polynomials to be fitted: ";
    cin>>nfuncs;
    cout<<"Performing x|theta corrections"<<endl;
    fit f(nfuncs);
    f.run(phisto, pcorr);
    cout<<"Corrections complete."<<endl;
  } if (options.runAll || options.cleanBackground) {
    cout<<"Cleaning up the corrected position histogram..."<<endl;
    cout<<"Corrections: "<<pcorr<<" Cleaned: "<<pclean<<endl;
    Backgnd destroy;
    destroy.run(pcorr, pclean);
    cout<<"Background annihiliated."<<endl;
  } 
  cout<<"SPS analysis complete."<<endl;
  return 0;
}
