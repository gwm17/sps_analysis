/*main.cpp
 *Main function for sps analysis program
 *3 modes: -r run everything, -a only standard analysis, -f only aberration corrections
 *-r & -f -> datafile, histogram file, fit file
 *-a -> datafile, histogram file
 *
 * Gordon M. Feb 2019
 */

#include "analysis.h"
#include "fit.h"
#include "TROOT.h"
#include "TApplication.h"
#include <iostream>
#include <unistd.h>

using namespace std;

//argument flags
struct options {
  int onlyFit; // -f
  int onlyAnalyze; // -a
  int runAll; // -r
} options;

//flag string for getopt; if expecting value with flag use : after flag letter
static const char *optString = "far";

int main(int argc, char* argv[]) {
  int opt = 0;
  options.onlyFit = 0;
  options.onlyAnalyze = 0;
  options.runAll = 0;
  
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
    }
    opt =  getopt(argc, argv, optString); // iterate to next arg
  }
 
  if ((options.runAll || options.onlyFit) && argc == 5) {
    int nfuncs;
    cout<<"Running SPS analysis with aberration correction..."<<endl;
    cout<<"Data: "<<argv[2]<<" Histograms: "<<argv[3]<<" Corrections: "<<argv[4]<<endl;
    TApplication app("app", &argc, argv);
    if(options.runAll) {
      cout<<"Sorting data..."<<endl;
      analysis a;
      a.run(app.Argv(2), app.Argv(3));
    }
    cout<<"Enter number of polynomials to be fitted: ";
    cin>>nfuncs;
    cout<<"Performing x|theta corrections"<<endl;
    fit f(nfuncs);
    f.run(app.Argv(2), app.Argv(4), app.Argv(3));
    cout<<"Finished"<<endl;
    return 0;
  } else if (options.onlyAnalyze && argc == 4) {
    cout<<"Running SPS analysis..."<<endl;
    cout<<"Data: "<<argv[2]<<" Histograms: "<<argv[3]<<endl;
    TApplication app("app", &argc, argv);
    analysis a;
    a.run(app.Argv(2), app.Argv(3));
    return 0;
  } else {
    cout << "Error: command line arguments not used correctly" << endl;
    cout << "See README.md for details"<< endl;
    return 0;
  }
}
