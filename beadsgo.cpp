#include "beads.h"

int main(int argc, char *argv[]) {

  int seed;
  seed = atoi(argv[1]); // random number seed
  srand48(seed);

  int yrs,nreps;
  double lambda,p,overdeath;
  yrs = atoi(argv[2]); // age of individual (in years)
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
  
  double fp;
  fp = atof(argv[9]);
  double nuid;
  nuid = atof(argv[10]);

  serialmanysimm(p,lambda,yrs,model,symm,overdeath,nreps,fp,nuid);

  return 1;
}

