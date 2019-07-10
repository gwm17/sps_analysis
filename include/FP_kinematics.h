/*

  Functions for the calculation of the kinematic shift of the FP
  for the SESPS @ FSU.


>>>  Get_Mass(Z, A) returns the mass (in amu) of an atom (Z, A).
     ***this requires a text file mass_info.txt w/ tabulated mass
     values.


>>>  Delta_Z(int...) returns the shift of the FP in the z-direction in
     cm. A negative (<0) delta-z is defined as a shift towards the
     magnet.
     ***NOTE: this uses Get_Mass(...), hence requires a tabulated
     mass file.

     Arguments: Delta_Z(int ZT, int AT, int ZP, int AP, int ZE, int AE,
	                double EP, double angle, double B),
        where Z,A are atomic number and mass number for each particle,
        EP is the KE of the projectile (i.e. beam energy in MeV),
        angle is the spectrograph angle (in degrees),
        B is the field in Gauss.


>>>  Delta_Z(double...) does the same, only takes in masses (in amu)
     explicitly, rather than (Z, A) for each particle (e.g. if you
     don't have a text file of tabulated masses).
     ***Make sure these are NUCLEAR masses, w/ electrons taken out.

     Arguments: Delta_Z(double MT, double MP, double ME, double MR,
   	                double EP, int ZE, double angle, double B),
        where NUCLEAR masses M are in amu, and the rest are the same
        as above. ***NOTE the required input of the ejectile Z for
	charge determination.


>>>  Wire_Dist() returns the distance (in cm) between the wires for
     the calculation of relative weights between FP1 and FP2.

  //format: T(P,E)R
  //   i.e., T = target,
  //         P = projectile,
  //         E = ejectile,
  //         R = residual;
  //expects angle in degs, B in G, masses and KE in MeV

  KGH -- Jul19

*/

#ifndef FP_KINEMATICS
#define FP_KINEMATICS

#include <iostream>
#include <fstream>
#include <cmath>

using namespace std;

double Get_Mass(int iZ, int iA) {

  //get masses
  ifstream infile;
  infile.open("mass_info.txt");

  if (!infile) {
    cerr << "***WARNING: cannot find file of masses (mass_info.txt); returning 0\n";
    return 0;
  }

  int Zcheck=0, Acheck=0;
  double mass=0;

  char junk[100];

  infile.getline(junk,100); //header

  while (Zcheck != -1) {

    infile >> Zcheck >> Acheck >> mass;

    if (Zcheck == iZ && Acheck == iA) {
      infile.close();
      return mass;
    }

  }

  //if we reach here, we did not find a mass for (iZ, iA)
  return 0;

}

//requires (Z,A) for T, P, and E, as well as energy of P,
// spectrograph angle of interest, and field value
double Delta_Z(int ZT, int AT, int ZP, int AP, int ZE, int AE,
	       double EP, double angle, double B) {

  /* CONSTANTS */
  const double UTOMEV = 931.4940954; //MeV per u;
  const double MEVTOJ = 1.60218E-13; //J per MeV
  const double RESTMASS_ELECTRON = 0.000548579909; //amu
  const double UNIT_CHARGE = 1.602E-19; //Coulombs
  const double C = 2.9979E8; //m/s

  /* SESPS-SPECIFIC */
  const double DISP = 1.96; //dispersion (x/rho)
  const double MAG = 0.39; //magnification in x
  const double DEGTORAD = M_PI/180.;

  int ZR = ZT + ZP - ZE, AR = AT + AP - AE;
  double EE=0; //ejectile energy

  double MT=0, MP=0, ME=0, MR=0; //masses (MeV)

  B /= 10000; //convert to tesla
  angle *= DEGTORAD;

  MT = (Get_Mass(ZT, AT) - ZT*RESTMASS_ELECTRON)*UTOMEV;
  MP = (Get_Mass(ZP, AP) - ZP*RESTMASS_ELECTRON)*UTOMEV;
  ME = (Get_Mass(ZE, AE) - ZE*RESTMASS_ELECTRON)*UTOMEV;
  MR = (Get_Mass(ZR, AR) - ZR*RESTMASS_ELECTRON)*UTOMEV;
  
  if (MT*MP*ME*MR == 0) {
    cerr << "***WARNING: error loading one or more masses; returning 0\n";
    return 0;
  }

  double Q = MT + MP - ME - MR; //Q-value
  
  //kinematics a la Iliadis p.590
  double term1 = sqrt(MP*ME*EP)/(ME + MR)*cos(angle);
  double term2 = (EP*(MR - MP) + MR*Q)/(ME + MR);

  EE = term1 + sqrt(term1*term1 + term2);
  EE *= EE;

  //momentum
  double PE = sqrt(EE*(EE+2*ME));

  //calculate rho from B a la B*rho = (proj. momentum)/(proj. charge)
  double rho = (PE*MEVTOJ)/(ZE*UNIT_CHARGE*C*B)*100; //in cm

  double K;

  K  = sqrt(MP*ME*EP/EE);
  K *= sin(angle);

  double denom = ME + MR - sqrt(MP*ME*EP/EE)*cos(angle);

  K /= denom;

  return -1*rho*DISP*MAG*K; //delta-Z in cm

}

//requires masses (in amu) for all particles, as well as energy of P,
// Z of E, spectrograph angle of interest, and rho of interest
//***make sure this is NUCLEAR mass, w/ electrons taken out
double Delta_Z(double MT, double MP, double ME, double MR,
	       double EP, int ZE, double angle, double B) {

  /* CONSTANTS */
  const double UTOMEV = 931.4940954; //MeV per u;
  const double MEVTOJ = 1.60218E-13; //J per MeV
  const double DEGTORAD = M_PI/180.;
  const double UNIT_CHARGE = 1.602E-19; //Coulombs
  const double C = 2.9979E8; //m/s

  /* SESPS-SPECIFIC */
  const double DISP = 1.96; //dispersion (x/rho)
  const double MAG = 0.39; //magnification in x

  double EE=0; //energy of ejectile

  //convert to units for calculation
  MT *= UTOMEV;
  MP *= UTOMEV;
  ME *= UTOMEV;
  MR *= UTOMEV;
  B /= 10000; //tesla
  angle *= DEGTORAD;

  double Q = MT + MP - ME - MR; //Q-value
  
  //kinematics for EE a la Iliadis p.590
  double term1 = sqrt(MP*ME*EP)/(ME + MR)*cos(angle);
  double term2 = (EP*(MR - MP) + MR*Q)/(ME + MR);

  EE = term1 + sqrt(term1*term1 + term2);
  EE *= EE;

  //momentum
  double PE = sqrt(EE*(EE+2*ME));

  //calculate rho from B a la B*rho = (proj. momentum)/(proj. charge)
  double rho = (PE*MEVTOJ)/(ZE*UNIT_CHARGE*C*B)*100; //in cm

  //now for actual delta-z calculation
  double K; //kinematic factor defined in Enge1979

  K  = sqrt(MP*ME*EP/EE);
  K *= sin(angle);

  double denom = ME + MR - sqrt(MP*ME*EP/EE)*cos(angle);

  K /= denom;

  return -1*rho*DISP*MAG*K; //delta-Z in cm

}

double Wire_Dist() {return 4.28625;} //cm

#endif
