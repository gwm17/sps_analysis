#include "analysis.h"
#include "fit.h"
#include "TROOT.h"
#include "TApplication.h"
#include <iostream>

using namespace std;

int main(int argc, char* argv[]) {
  if (argc == 4) {
    int nfuncs;
    cout<<"Running SPS analysis with aberration correction..."<<endl;
    cout<<"Data: "<<argv[1]<<" Histograms: "<<argv[2]<<" Corrections: "<<argv[3]<<endl;
    TApplication app("app", &argc, argv);
    analysis a;
    a.run(app.Argv(1), app.Argv(2));
    cout<<"Enter number of polynomials to be fitted: ";
    cin>>nfuncs;
    cout<<"Performing x|theta corrections"<<endl;
    fit f(nfuncs);
    f.run(app.Argv(1), app.Argv(3), app.Argv(2));
    cout<<"Finished"<<endl;
    return 0;
  } else if (argc == 3) {
    cout<<"Running SPS analysis..."<<endl;
    cout<<"Data: "<<argv[1]<<" Histograms: "<<argv[2]<<endl;
    TApplication app("app", &argc, argv);
    analysis a;
    a.run(app.Argv(1), app.Argv(2));
    return 0;
  } else {
    cout << "Error: Application analysis requires a minimum 2 names, at most 3" << endl;
    cout << "1. Data name 2. Histogram name 3. Fit name"<< endl;
    return 0;
  }
}
