
#include <iostream>
#include <fstream>
#include <iomanip>
#include <cmath>
#include <stdlib.h>
using namespace std;

#ifndef BEAD5_H
#define BEAD5_H

#define DECREASEAGE 0 // 1 incorporates cell death, 0 does not

double rlogunif(double lo,double hi) {
  double xx;
  xx = pow(10.0,log10(lo) + ((log10(hi)-log10(lo))*drand48()));
  return xx;
}

int rpoisson(double mean) {
  double temp;
  int n,flag;
  n = 0;
  temp = 0.0;
  flag = 1;
  while (flag) {
    n += 1;
    temp += -log(drand48());
    if (temp > mean) {
      n -= 1;
      flag = 0;
    }
  }
  return n;
}

double rnormal(double mean,double sd) {
  double y1,y2,x;
  int flag;
  flag = 1;
  while (flag) {
    y1 = -log(drand48());
    y2 = -log(drand48());
    if (y2 >= (0.5*(y1-1.0)*(y1-1.0))) {
      if (drand48() < 0.5) { x = y1; }
      else { x = -y1; }
      flag = 0;
    }
  }
  x = mean + (x*sd);
  return x;
}

double myrand(double n,double p) {
  double lambda,x;
  int i;

  lambda = n*p;
  if ((n>=100.0) && (lambda<=10.0)) { // POISSON APPROXIMATION
    x = ((double)(rpoisson(lambda)));
  }
  else {
    lambda = (double(n))*p*(1-p);
    if (lambda>=10.0) { // NORMAL APPROXIMATION
      x = ((double(n))*p) + (sqrt(lambda)*rnormal(0.0,1.0));
    }
    else { // NO APPROXIMATION
      x = 0.0;
      for (i=0;i<n;i++) {
	if (drand48() < p) { x += 1.0; }
      }
    }
  }

  if (x<0.0) { x = 0.0; }
  if (x>n) { x = n; }

  return x;
}

double grow(double sta,int stime,int etime,double p,double symm,double overdeath,int model) { 
// OUTPUT NUMBER OF MUTANTS IN THE CLUSTER
// INPUT INITIAL NUMBER MUTANTS (1 FOR ADULT PHASE, MORE FOR GROWTH PHASE)
// stime and etime in cell generations not years

  int i;
  double x,y,z,d,z2;
  x = sta;

  if (model == 0) { // ASYMMETRIC HOTSPOT/SELECTION MODELS
    if (p > 0.0) {
      for (i=stime;i<etime;i++) {
	x += myrand(x,p);
      }    
    }
  }
  if (model == 1) { // SYMMETRIC HOTSPOT MODEL
    if (symm > 0.0) {
      for (i=stime;i<etime;i++) {
	y = myrand(x,symm);
	z = myrand(x,symm);
	x = x + y - z;
	if ( x <= 0.0) {
	  x = 0.0;
	  break;
	}
      }
    }
  }
  if (model == 2) { // ASYMMETRIC SELECTION MODEL INCORPORATING CELL DEATH
    for (i=stime;i<etime;i++) {
      d = 0.0; // CAN SPEED UP BY ONLY CHECKING EVERY ~10 YEARS
      if (i >= ((35-13)*23)) { d = 0.00013; }
      if (i >= ((45-13)*23)) { d = 0.00086; }
      if (i >= ((55-13)*23)) { d = 0.00032; }
      if (i >= ((65-13)*23)) { d = 0.00170*overdeath; }
      if (i >= ((75-13)*23)) { d = 0.00150*overdeath; }

      y = myrand(x,p);
      z = myrand(x,d);
      x = x + y - z;
      if (x <= 0.0) {
	x = 0.0;
	break;
      }
    }
  }
  if (model == 3) { // SYMMETRIC SELECTION MODEL INCORPORATING CELL DEATH
    for (i=stime;i<etime;i++) {
      d = 0.0; // CAN SPEED UP BY ONLY CHECKING EVERY ~10 YEARS
      if (i >= ((35-13)*23)) { d = 0.00013; }
      if (i >= ((45-13)*23)) { d = 0.00086; }
      if (i >= ((55-13)*23)) { d = 0.00032; }
      if (i >= ((65-13)*23)) { d = 0.00170*overdeath; }
      if (i >= ((75-13)*23)) { d = 0.00150*overdeath; }

      y = myrand(x,(0.50+p));
      if (DECREASEAGE == 1) { z = myrand(x,d); }
      if (DECREASEAGE == 0) { z = 0.0; }
      z2 = myrand(x,(0.50-p));
      x = x + y - z - z2;
      if (x <= 0.0) {
	x = 0.0;
	break;
      }
    }
  }
  
  return x;
}

void space(double testis[6][8][4],double mclump,double mtot) {
  // OUTPUT DISTRIBUTES MUTATION CLUSTER INTO 3-DIMENSIONAL TESTIS ARRAY
  // INPUT 3-DIMENSIONAL TESTIS ARRAY,  
  // NUMBER OF MUTANTS IN CLUMP, TOTAL NUMBER OF CELLS IN TESTIS

  double x,y,z,x1,x2,y1,y2,z1,z2,leng,dtemp;
  int i,j,k;

  x = 6.*drand48(); 
  y = 8.*drand48(); 
  z = 4.*drand48();
  leng = pow( 192.*mclump/mtot, 0.333);

  if (leng > 4) {
    dtemp = round(mclump/192.);
    for (i=0;i<6;i++) {
      for (j=0;j<8;j++) {
	for (k=0;k<4;k++) {
	  testis[i][j][k] += dtemp;
	  if (testis[i][j][k] > (mtot/192.)) {
	    testis[0][0][0] = -9;
	  }
	}
      }
    }
  }
  else {
    if (x-(0.5*leng) < 0) { x += (0.5*leng)-x; }
    if (x+(0.5*leng) > 6) { x -= x+(0.5*leng)-6; }
    if (y-(0.5*leng) < 0) { y += (0.5*leng)-y; }
    if (y+(0.5*leng) > 8) { y -= y+(0.5*leng)-8; }
    if (z-(0.5*leng) < 0) { z += (0.5*leng)-z; }
    if (z+(0.5*leng) > 4) { z -= z+(0.5*leng)-4; }
    for (i=0;i<6;i++) {
      for (j=0;j<8;j++) {
	for (k=0;k<4;k++) {
	  x1 = max( ((double)(i)),  (x-(0.5*leng)) );
	  x2 = min( ((double)(i+1)),(x+(0.5*leng)) );
	  y1 = max( ((double)(j)),  (y-(0.5*leng)) );
	  y2 = min( ((double)(j+1)),(y+(0.5*leng)) );
	  z1 = max( ((double)(k)),  (z-(0.5*leng)) );
	  z2 = min( ((double)(k+1)),(z+(0.5*leng)) );
	  if ((x2 > x1) && (y2 > y1) && (z2 > z1)) {
	    dtemp = (x2-x1)*(y2-y1)*(z2-z1);
	    testis[i][j][k] += round( mclump*dtemp/pow(leng,3.) );
	    if (testis[i][j][k] > (mtot/192.) ) {
	      testis[0][0][0] = -9; 
	    }
	  }
	}
      }
    }
  }
}

void simm(double p,double lambda,int yrs,double testis[6][8][4],int model,double symm,double overdeath,double fp,double nuid) {
  int growgens = 30;
  
  int iii,jjj,cnt,divyr,adultgens;
  adultgens = (yrs-13)*23;

  int i,j,k,itime;
  double mtot,extra,dtemp,muts,x,totmut,d,fpm,rm;
  totmut = 0.0; 
  mtot = round(pow(2,((double)(growgens))));

  // FORMAT
  extra = 0;
  for (i=0;i<6;i++) {
    for (j=0;j<8;j++) {
      for (k=0;k<4;k++) {
	testis[i][j][k] = 0.0;
      }
    }
  }
      
  // GROWTH PHASE
  for (i=0;i<growgens;i++) {
    dtemp = pow(2,((double)(i+1)));
    muts = myrand(dtemp,lambda);
    if (muts > 0) {
      dtemp = pow(2,((double)(growgens-i-1)));
      for (j=0;j<muts;j++) {
	x = grow(dtemp,0,adultgens,p,symm,overdeath,model);
	totmut += x; 
	extra += (x-dtemp); 
	space(testis,x,mtot);
	if (testis[0][0][0] < 0.0) { break; }
      }
    }
  }
  
  // ADULT PHASE
  if (testis[0][0][0] >= 0.0) {
    for (i=0;i<adultgens;i++) {
      if (DECREASEAGE == 1) {
	d = 0.0; // CAN SPEED UP BY ONLY CHECKING EVERY ~10 YEARS
	if (i >= ((35-13)*23)) { d = 0.00013; }
	if (i >= ((45-13)*23)) { d = 0.00086; }
	if (i >= ((55-13)*23)) { d = 0.00032; }
	if (i >= ((65-13)*23)) { d = 0.00170; }
	if (i >= ((75-13)*23)) { d = 0.00150; }

	mtot -= myrand(mtot,d);
      }
      muts = myrand(mtot,lambda);
      if (muts > 0) {
	for (j=0;j<muts;j++) {
	  x = grow(1.0,i,adultgens,p,symm,overdeath,model);
	  totmut += x;
	  extra += x; 
	  space(testis,x,mtot);
	  if (testis[0][0][0] < 0.0) { break; }
	}
      }
    }
  }

  if (mtot > 0.0) {
    for (i=0;i<6;i++) {
      for (j=0;j<8;j++) {
	for (k=0;k<4;k++) {
	
	  testis[i][j][k] = testis[i][j][k]/(mtot/192.0);
	  fpm = myrand(nuid,fp);
	  rm = myrand(nuid,testis[i][j][k]);
	  testis[i][j][k] = (fpm+rm)/nuid;

	}
      }
    }
  }

}

void getx5(double testis[6][8][4],double x5[]) {
  // output: 1 mutants in top 5 pieces, 2 mutants in rest of testis
  int i,j,k,cnt,ibig;
  double big,lin[192];

  cnt = -1;
  for (i=0;i<6;i++) {
    for (j=0;j<8;j++) {
      for (k=0;k<4;k++) {
	cnt += 1;
	lin[cnt] = testis[i][j][k];
      }
    }
  }

  for (i=0;i<191;i++) {
    big = lin[i];
    for (j=(i+1);j<192;j++) {
      if (lin[j] > big) {
	big = lin[j];
	ibig = j;
      }
    }
    if (big > lin[i]) {
      lin[ibig] = lin[i];
      lin[i] = big;
    }
  }

  j = 5; // can change here
  x5[0] = 0.0;
  x5[1] = 0.0;
  for (i=0;i<j;i++) {
    if (lin[i] >= 0) {
      x5[0] += lin[i];
    }
  }
  for (i=j;i<192;i++) {
    if (lin[i] >= 0) {
      x5[1] += lin[i];
    }
  }
}

double getp95(double testis[6][8][4]) {
  // OUTPUT THE PERCENTAGE OF TESTIS PIECES REQUIRED TO COMPRISE 95% OF THE MUTANTS
  // INPUT 3-DIMENSIONAL TESTIS ARRAY

  int i,j,k,cnt,ibig;
  double big,tot,lin[192],good; good = 0.0;

  tot = 0.0;
  cnt = -1;
  for (i=0;i<6;i++) {
    for (j=0;j<8;j++) {
      for (k=0;k<4;k++) {
	cnt += 1;
	lin[cnt] = testis[i][j][k];
	if (testis[i][j][k] >= 0) {
	  tot += testis[i][j][k];
	  good += 1.0;
	}
      }
    }
  }
  for (i=0;i<191;i++) {
    big = lin[i];
    for (j=(i+1);j<192;j++) {
      if (lin[j] > big) {
	big = lin[j];
	ibig = j;
      }
    }
    if (big > lin[i]) {
      lin[ibig] = lin[i];
      lin[i] = big;
    }
  }

  cnt = -1;
  if (tot > 0) {
    big = 0.;
    while ((big/tot) < 0.95) {
      cnt += 1; if (cnt > 191) { cnt = 191; break; }
      big += lin[cnt];
    }
  }

  return ((double)(cnt+1))/good;
}
     
void stat4(double testis[6][8][4],double xstat[]) {
  int i,j,k;
  double good;
  double myx5[2];
  good = 0.0;

  if (testis[0][0][0] < 0.0) {
    xstat[0] = -9;
    xstat[1] = -9;
    xstat[2] = -9;
    xstat[3] = -9;
    xstat[4] = -9;
    xstat[5] = -9;
    xstat[6] = -9;
  }
  else {
    xstat[0] = 0.0;
    xstat[2] = testis[0][0][0];
    xstat[3] = testis[0][0][0];
    xstat[4] = 0.0;
    for (i=0;i<6;i++) {
      for (j=0;j<8;j++) {
	for (k=0;k<4;k++) {
	  if (testis[i][j][k] >= 0) {
	    good += 1.0;
	    xstat[0] += testis[i][j][k];
	    if (testis[i][j][k] > xstat[2]) {
	      xstat[2] = testis[i][j][k];
	    }
	    if (testis[i][j][k] < xstat[3]) {
	      xstat[3] = testis[i][j][k];
	    }
	    if ((testis[i][j][k]*1000000.) < 50.) {
	      xstat[4] += 1.0;
	    }
	  }
	}
      }
    }
    xstat[0] = (xstat[0]/good)*1000000.; // FREQ
    xstat[2] = xstat[2]*1000000.; // MAX
    xstat[3] = xstat[3]*1000000.; // MIN
    xstat[1] = getp95(testis);
    xstat[4] = xstat[4]/good; // frac pieces <= 50 mutants per million molecules
    getx5(testis,myx5);
    xstat[5] = myx5[0];
    xstat[6] = myx5[1];
  }
}

void onesimm(double p,double lambda,int yrs,int model,double symm,double overdeath,double fp,double nuid) {
  double testis[6][8][4];
  int i,j,k;

  simm(p,lambda,yrs,testis,model,symm,overdeath,fp,nuid);
  for (i=0;i<6;i++) {
    for (j=0;j<8;j++) {
      for (k=0;k<4;k++) {
	cout << testis[i][j][k] << endl;
      }
    }
  }
}

void serialmanysimm(double p,double lambda,int yrs,int model,double symm,double overdeath,int nreps,double fp,double nuid) {

  double testis[6][8][4],xstat[7];
  int flag;

  for (flag=0;flag<nreps;flag++) {
      simm(p,lambda,yrs,testis,model,symm,overdeath,fp,nuid);
      stat4(testis,xstat);
      // OUTPUT
      cout << yrs << " " << model << " " << p << " " << lambda << " " << xstat[0] << " " << xstat[2] << " " << fp << " " << nuid << endl;
  }
}

void allpieces(double testis[6][8][4],double allp[]) {
  int i,j,k,ii;
  ii = 0;
  for (i=0;i<6;i++) {
    for (j=0;j<8;j++) {
      for (k=0;k<4;k++) {
	allp[ii] = testis[i][j][k];
	ii = ii + 1;
      }
    }
  }
}

void serialgetparam(double p,double lambda,int yrs,int model,double symm,double overdeath,int nreps,double targetav,double targetmx,double targetf50,double targetp95,double fp,double nuid) {

  // output fraction of simulations within 50% and 5% of Av target value
  
  double testis[6][8][4],xstat[7];
  int i,deli,flag,cont,half,top;
  double m50,m05;
  m50 = 0.0;
  m05 = 0.0;
  double xhh,x55,xh5,x5h;
  xhh = 0.0;
  x55 = 0.0;
  xh5 = 0.0;
  x5h = 0.0;
  double xuu,xud,xdu,xdd,bad,ratdat,ratsim;
  bad = 0.0;
  xuu = 0.0;
  xud = 0.0;
  xdu = 0.0;
  xdd = 0.0;
  double x5hhz,x5hzh,x5hhh;
  x5hhz = 0.0;
  x5hzh = 0.0;
  x5hhh = 0.0;
  double x55h,x5h5,x555;
  x55h = 0.0;
  x5h5 = 0.0;
  x555 = 0.0;
  double mydist,myavedist;
  myavedist = 0.0;

  double mxavedist,mx50,mx05;
  mxavedist = 0.0;
  mx50 = 0.0;
  mx05 = 0.0;

  for (flag=0;flag<nreps;flag++) {
      simm(p,lambda,yrs,testis,model,symm,overdeath,fp,nuid);
      stat4(testis,xstat);
      
      if (xstat[0] >= 0.0) {
	if ( ((targetav*(1-0.5)) <= xstat[0]) && (xstat[0] <= (targetav*(1+0.5))) ) { 
	  m50 += 1.0; 
	}
	if ( ((targetav*(1-0.05)) <= xstat[0]) && (xstat[0] <= (targetav*(1+0.05))) ) { 
	  m05 += 1.0; 
	}
	mydist = targetav - xstat[0];
	if (mydist < 0.0) { mydist = (-1) * mydist; } 
	myavedist = myavedist + mydist;


	if ( ((targetmx*(1-0.5)) <= xstat[2]) && (xstat[2] <= (targetmx*(1+0.5))) ) { 
	  mx50 += 1.0; 
	}
	if ( ((targetmx*(1-0.05)) <= xstat[2]) && (xstat[2] <= (targetmx*(1+0.05))) ) { 
	  mx05 += 1.0; 
	}
	mydist = targetmx - xstat[2];
	if (mydist < 0.0) { mydist = (-1) * mydist; } 
	mxavedist = mxavedist + mydist;
 
      }
      else { 
	bad += 1.0;
      } 
  }

  m50 = m50 / (double(nreps));
  m05 = m05 / (double(nreps));
  myavedist = myavedist / (double(nreps));

  mx50 = mx50 / (double(nreps));
  mx05 = mx05 / (double(nreps));
  mxavedist = mxavedist / (double(nreps));

  // OUTPUT
  cout << yrs << " " << model << " " << p << " " << lambda << " " << symm << " " << overdeath << " HERE " << m50 << " " << m05 << " " << myavedist << " " << fp << " " << nuid << " " << mx50 << " " << mx05 << " " << mxavedist << " " << endl;
}

#endif
