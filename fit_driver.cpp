#include <TApplication.h>
#include "fit.h"

using namespace std;

int main(int argc, char* argv[]) {

  TApplication myApp("myApp",&argc,argv);

  fit bob;
  bob.run("output_files/12C_dp_20deg_bb_centered.root",
	  "output_files/12C_dp_20deg_bb_centered_poly.root",
	  "output_files/12C_dp_20deg_bb_centered.root");

  return 0;

}
