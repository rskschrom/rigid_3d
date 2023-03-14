#include <iostream>
#include <iomanip>
#include <fstream>
#include <iterator>
#include <math.h>
#include <vector>
#include "params.h"
#include "state.h"
#include "init.h"
#include "physics.h"
#include "geometry.h"

int main(void){

  int i, j, k, npar1, npar2, cpi1, cpi2;
  float pair_dist, mass1, mass2, fx1_sum, fy1_sum, fx2_sum, fy2_sum;
  float t1_sum, t2_sum, ax1, ay1, ax2, ay2, alp1, alp2;
  float inerm1, inerm2, ke1, ke2, rke1, rke2;

  // initial body properties
  struct Params p;
  struct State s1, s2;
  std::vector<float> fx1, fy1, fx2, fy2;
  std::vector<float> xtr1_out, ytr1_out, xtr2_out, ytr2_out;
  std::vector<float> fx1_out, fy1_out, fx2_out, fy2_out, dis_vec;
  std::vector<int> ip1, ip2, hi1, hi2;

  s1.theta = 0.;
  s1.omega = 0.2;
  s1.vx = -0.3;
  s1.vy = 0.;
  s1.xcom = 0.75;
  s1.ycom = 3.;
  s1.rho = 920.;
  s1.dm = pow(p.ds, 2)*s1.rho;
  s1.nhx = 9;
  s1.nhy = 9;
    
  s2.theta = p.pi/6.;
  s2.omega = 0.;
  s2.vx = 0.;
  s2.vy = 0.;
  s2.xcom = 0.;
  s2.ycom = -2.;
  s2.rho = 920000.;
  s2.dm = pow(p.ds, 2)*s2.rho;
  s2.nhx = 9;
  s2.nhy = 9;
  
  //init_body(p, s1, 30, 10);
  //init_body(p, s2, 10, 10);
  
  // define particle bodies and transform them
  init_body_file(p, s1, "par_data2.txt");
  init_body_file(p, s2, "par_data2.txt");
  transform_particles(p, s1);
  transform_particles(p, s2);

  npar1 = s1.x.size();
  npar2 = s2.x.size();
    
  spatial_hash(p, s1);
  spatial_hash(p, s2);
  
  printf("%d %d\n", npar1, npar2);
  
  // write out initial particle base positions
  std::ofstream file_x1("x1.txt");
  std::ostream_iterator<float> fiterx1(file_x1, "\n");
  std::copy(s1.x.begin(), s1.x.end(), fiterx1);
 
  std::ofstream file_y1("y1.txt");
  std::ostream_iterator<float> fitery1(file_y1, "\n");
  std::copy(s1.y.begin(), s1.y.end(), fitery1);

  std::ofstream file_x2("x2.txt");
  std::ostream_iterator<float> fiterx2(file_x2, "\n");
  std::copy(s2.x.begin(), s2.x.end(), fiterx2);

  std::ofstream file_y2("y2.txt");
  std::ostream_iterator<float> fitery2(file_y2, "\n");
  std::copy(s2.y.begin(), s2.y.end(), fitery2);
  
  std::ofstream fileh1("hash1.txt");
  std::ostream_iterator<int> fiterh1(fileh1, "\n");
  std::copy(s1.hi.begin(), s1.hi.end(), fiterh1);
  
  std::ofstream fileh2("hash2.txt");
  std::ostream_iterator<int> fiterh2(fileh2, "\n");
  std::copy(s2.hi.begin(), s2.hi.end(), fiterh2);
  
  std::ofstream filehx1("hx1.txt");
  std::ostream_iterator<float> fiterhx1(filehx1, "\n");
  std::copy(s1.xh.begin(), s1.xh.end(), fiterhx1);
  
  std::ofstream filehy1("hy1.txt");
  std::ostream_iterator<float> fiterhy1(filehy1, "\n");
  std::copy(s1.yh.begin(), s1.yh.end(), fiterhy1);
  
  std::ofstream filehx2("hx2.txt");
  std::ostream_iterator<float> fiterhx2(filehx2, "\n");
  std::copy(s2.xh.begin(), s2.xh.end(), fiterhx2);
  
  std::ofstream filehy2("hy2.txt");
  std::ostream_iterator<float> fiterhy2(filehy2, "\n");
  std::copy(s2.yh.begin(), s2.yh.end(), fiterhy2);
  
  // output files for rigid body position and orientation
  std::ofstream file_r1("rigid1.txt");
  std::ofstream file_r2("rigid2.txt");
  
  file_r1 << std::setw(9) << std::setprecision(4) << std::fixed << s1.xcom
            << std::setw(9) << std::setprecision(4) << std::fixed << s1.ycom 
            << std::setw(9) << std::setprecision(4) << std::fixed << s1.theta << std::endl;
            
  file_r2 << std::setw(9) << std::setprecision(4) << std::fixed << s2.xcom
            << std::setw(9) << std::setprecision(4) << std::fixed << s2.ycom 
            << std::setw(9) << std::setprecision(4) << std::fixed << s2.theta << std::endl;
  
  // calculate inertia moments and mass
  inerm1 = 0.;
  inerm2 = 0.;
  mass1 = 0.;
  mass2 = 0.;
  
  for (i = 0; i < npar1; i++){
    mass1 += s1.dm;
    inerm1 += s1.dm*(pow(s1.x[i], 2)+pow(s1.y[i], 2));
  }
  
  for (i = 0; i < npar2; i++){
    mass2 += s2.dm;
    inerm2 += s2.dm*(pow(s2.x[i], 2)+pow(s2.y[i], 2));
  }
  
  s1.mass = mass1;
  s2.mass = mass2;
  s1.inerm = inerm1;
  s2.inerm = inerm2;
  
  // time loop  
  for (k = 0; k < p.nt; k++){
    printf("%d\n", k);    
    // detect collision
    lj_collision(p, s1, s2, 0., 0., 0.1*p.g, 0., 0., 0.);
    
    evolve_motion(p, s1, 0., 0., 0.1*p.g, p.dt),
    evolve_motion(p, s2, 0., 0., 0., p.dt);
    
    // file output
    file_r1 << std::setw(9) << std::setprecision(4) << std::fixed << s1.xcom
            << std::setw(9) << std::setprecision(4) << std::fixed << s1.ycom 
            << std::setw(9) << std::setprecision(4) << std::fixed << s1.theta << std::endl;
    
    file_r2 << std::setw(9) << std::setprecision(4) << std::fixed << s2.xcom
            << std::setw(9) << std::setprecision(4) << std::fixed << s2.ycom 
            << std::setw(9) << std::setprecision(4) << std::fixed << s2.theta << std::endl;
    
    // assess energy conservation
    ke1 = 0.5*mass1*(s1.vx*s1.vx+s1.vy*s1.vy);
    ke2 = 0.5*mass2*(s2.vx*s2.vx+s2.vy*s2.vy);
    rke1 = 0.5*inerm1*s1.omega*s1.omega;
    rke2 = 0.5*inerm2*s2.omega*s2.omega;
    //printf("%8.6f %8.6f %8.5f\n", ke1, ke2, ke1+ke2);
    //printf("%8.6f %8.6f %8.5f\n", rke1, rke2, rke1+rke2);
    //printf("%8.6f\n", ke1+ke2+rke1+rke2);

  }
  file_r1.close();
  file_r2.close();

}
