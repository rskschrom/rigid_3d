#include <iostream>
#include <fstream>
#include <string>
#include <math.h>
#include <vector>
#include "params.h"
#include "state.h"

// initialize body mass and particles
void init_body(Params p, State &s, int nx, int ny){

  int i, j;
  int npar = nx*ny;
  std::vector<float> x(npar);
  std::vector<float> y(npar);
    
  // set state variables
  s.mass = pow(p.ds, 2)*s.rho*npar;
    
  // arrange body particles on x-y grid
  for (i = 0; i < nx; i++){
    for (j = 0; j < ny; j++){
        x[i*ny+j] = (i-(nx-1)/2.)*p.ds;
        y[i*ny+j] = (j-(ny-1)/2.)*p.ds;
    }
  }
    
  s.x = x;
  s.y = y;
}

// initialize body mass and particles from text file
void init_body_file(Params p, State &s, std::string fname){

  int i, j;
  int npar;
  std::vector<float> x;
  std::vector<float> y;
  int xv, yv;
  float mnx, mny;
    
  // read position data from a file (indices; no header)
  std::ifstream ifile;
  ifile.open(fname);
  
  mnx = 0.;
  mny = 0.;
  while (!ifile.eof()){
    ifile >> xv;
    ifile >> yv;
    x.push_back(xv*p.ds);
    y.push_back(yv*p.ds);
    mnx += xv*p.ds;
    mny += yv*p.ds;
  }
  
  ifile.close();
  npar = x.size();
  mnx = mnx/npar;
  mny = mny/npar;
  
  // center particle at origin
  for (i = 0; i < npar; i++){
    x[i] = x[i]-mnx;
    y[i] = y[i]-mny;
  }
    
  // set state variables
  s.mass = pow(p.ds, 2)*s.rho*npar;
  s.x = x;
  s.y = y;
  s.xt.resize(npar);
  s.yt.resize(npar);
}
