#include <iostream>
#include <fstream>
#include <string>
#include <math.h>
#include <vector>
#include "params.h"
#include "state.h"
#include <Eigen/Dense>

// initialize body mass and particles from text file
void init_body_file(Params p, State &s, std::string fname){
  int npar;
  std::vector<float> x, y, z;
  int xv, yv, zv;
  float mnx, mny, mnz;
  Eigen::Matrix3f inerm, iinerm;
    
  // read position data from a file (indices; no header)
  std::ifstream ifile;
  ifile.open(fname);
  
  // skip 3 comment lines at beginning of file
  constexpr auto max_size = std::numeric_limits<std::streamsize>::max();
  ifile.ignore(max_size, '\n');
  ifile.ignore(max_size, '\n');
  ifile.ignore(max_size, '\n');
  
  mnx = 0.;
  mny = 0.;
  mnz = 0.;
  
  while (!ifile.eof()){
    ifile >> xv;
    ifile >> yv;
    ifile >> zv;
    x.push_back(xv*p.ds);
    y.push_back(yv*p.ds);
    z.push_back(zv*p.ds);
    
    mnx += xv*p.ds;
    mny += yv*p.ds;
    mnz += zv*p.ds;
  }
  
  ifile.close();
  npar = x.size();
  mnx = mnx/npar;
  mny = mny/npar;
  mnz = mnz/npar;
  s.r.resize(3*npar);
  
  // center particle com at origin
  for (int i = 0; i < npar; i++){
    x[i] = x[i]-mnx;
    y[i] = y[i]-mny;
    z[i] = z[i]-mnz;
    
    s.r[3*i] = x[i];
    s.r[3*i+1] = y[i];
    s.r[3*i+2] = z[i];
  }
  
  // set state variables and allocate vectors
  s.dm = pow(p.ds, 2)*s.rho;
  s.mass = s.dm*npar;
  s.rt.resize(3*npar);
  s.com.resize(3);
  s.vel.resize(3);
  s.inerm.resize(9);
  s.iinerm.resize(9);
  s.orient.resize(4);
  s.omega.resize(4);
  
  s.orient = {0.,0.,0.,1.};
  
  // calculate inertia moment tensor
  for (int i = 0; i < 3; i++){
    for (int j = 0; j < 3; j++){
      
      // tensor off diagonal
      if (i!=j){
        for (int k = 0; k < npar; k++){
          s.inerm[3*i+j] += -s.r[3*k+i]*s.r[3*k+j];
        }
      }
      
      // tensor diagonal
      else{
        for (int k = 0; k < npar; k++){
          s.inerm[3*i+j] += s.r[3*k]*s.r[3*k]+
                            s.r[3*k+1]*s.r[3*k+1]+
                            s.r[3*k+2]*s.r[3*k+2]-
                            s.r[3*k+i]*s.r[3*k+i];
        }
      }
      
      s.inerm[3*i+j] = s.dm*s.inerm[3*i+j];
      inerm(i,j) = s.inerm[3*i+j];
    }
  }
  
  // inverse of inertia tensor
  Eigen::EigenSolver<Eigen::Matrix3f> es(inerm);  
  Eigen::Matrix3f iD;
  Eigen::Matrix3f eig_vecs = es.eigenvectors().real();
  Eigen::Array<float, 1, 3> ieigv = es.eigenvalues().real().array().inverse();
  
  iD = ieigv.matrix().asDiagonal();
  iinerm = eig_vecs*iD*eig_vecs.transpose();
  
  for (int i = 0; i < 3; i++){
    for (int j = 0; j < 3; j++){
      s.iinerm[3*i+j] = iinerm(i,j);
    }
  }

}
