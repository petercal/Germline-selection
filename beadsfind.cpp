#include "beads.h"

int main(int argc, char *argv[]) {
  int seed;
  seed = atoi(argv[1]); // random number seed
  srand48(seed);

  int targetage,nreps;
  double lambda,p,overdeath;
  targetage = atoi(argv[2]); // age of individual (in years)
  lambda = atof(argv[3]);    // mutation rate per cell division
  p = atof(argv[4]);         // selection parameter
  overdeath = atof(argv[5]); // delta parameter
  nreps = atoi(argv[6]);     // number simulations

  double symm;
  int i,j,model;
  
  model = atoi(argv[7]);
  // 3 symmetric selection model, 2 asymmetric selection model incorporating cell death,
  // 1 symmetric hotspot model, 0 asymmetric hotspot/selection model

  
  symm = atof(argv[8]);   // symmetric parameter q
  

  double targetlambda,targetmx;
  targetlambda = atof(argv[9]);
  targetmx = atof(argv[10]);

  double targetf50,targetp95;
  targetf50 = atof(argv[11]);
  targetp95 = atof(argv[12]);

  double fp;
  fp = atof(argv[13]);
  double nuid;
  nuid = atof(argv[14]);

  // VARY BELOW
  double veclambda[40], vecp[40], vecsymm[20];
  veclambda[0]   = 0.000000000054;
  veclambda[1]   = 0.000000021;
  veclambda[2]   = 0.000000022;
  veclambda[3]   = 0.000000023;
  veclambda[4]   = 0.000000024;
  veclambda[5]   = 0.000000025;
  veclambda[6]   = 0.000000026;
  veclambda[7]   = 0.000000027;
  veclambda[8]   = 0.000000028;
  veclambda[9]   = 0.00000000051;
  veclambda[10]  = 0.00000000052;
  veclambda[11]  = 0.00000000053;
  veclambda[12]  = 0.00000000054;
  veclambda[13]  = 0.00000000055;
  veclambda[14]  = 0.00000000056;
  veclambda[15]  = 0.00000000057;
  veclambda[16]  = 0.00000000058;

  veclambda[17]  = 0.0000000500;
  veclambda[19]  = 0.000000000070;
  veclambda[20]  = 0.000000000071;
  veclambda[21]  = 0.000000000072;
  veclambda[22]  = 0.000000000073;
  veclambda[23]  = 0.000000000074;
  veclambda[24]  = 0.000000000075;
  veclambda[25]  = 0.0000000000036;
  veclambda[26]  = 0.0000000000037;
  veclambda[27]  = 0.0000000000038;
  veclambda[28]  = 0.0000000000039;
  veclambda[29]  = 0.0000000000040;
  veclambda[30]  = 0.0000000000041;
  veclambda[31]  = 0.0000000000042;
  veclambda[32]  = 0.0000000000043;
  veclambda[33]  = 0.0000000000044;
  veclambda[34]  = 0.0000000000045;

  vecp[0] = 0.0021;
  vecp[1] = 0.0022;
  vecp[2] = 0.0023;
  vecp[3] = 0.0024;
  vecp[4] = 0.0025;
  vecp[5] = 0.0026;
  vecp[6] = 0.0027;
  vecp[7] = 0.0028;
  vecp[8] = 0.0029;
  vecp[9] = 0.0030;
  vecp[10] = 0.0031;
  vecp[11] = 0.0032;
  vecp[12] = 0.0033;
  vecp[13] = 0.0034;
  vecp[14] = 0.0035;
  vecp[15] = 0.0036;
  vecp[16] = 0.0037;
  vecp[17] = 0.0038;
  vecp[18] = 0.0039;

  vecsymm[0] = symm;

  for (i=0;i<19;i++) {
    for (j=0;j<1;j++) {
      serialgetparam(vecp[i],veclambda[j],targetage,model,vecsymm[0],overdeath,nreps,targetlambda,targetmx,targetf50,targetp95,fp,nuid);
    }
  }

  return 1;
}

