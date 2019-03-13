#include <TApplication.h>
#include "fit.h"

using namespace std;

int main(int argc, char* argv[]) {

  TApplication myApp("myApp",&argc,argv);

  fit bob(1,4,3); //psvar order numpolys

  bob.run("output_files/12C_dp_20deg_th3ph2.root",
	  "corr_tree",
	  "output_files/12C_dp_20deg_th3ph2th4.root",
	  "corr_tree");

  return 0;

}
