#include <iostream>
#include <fmt/core.h>
#include <iomanip>
#include <fstream>
#include <iterator>
#include <math.h>
#include <vector>
#include "params.h"
#include "state.h"
#include "init.h"
#include "quat.h"
#include "physics.h"
#include "geometry.h"

int main(void){

  int npar, ncol, rcount;
  std::vector<float> alpha, acc, dis_vec, f, bbox(6), con(3);
  std::vector<float> dv(6), iitens_rot(9), traj;
  std::vector<int> ip;
  
  bool attached = false;

  // initial body properties
  struct Params p;
  struct State s;

  s.rho = 920.;
  
  // define particle bodies and transform them
  init_body_file(p, s, "crystal_0000_1.5.txt");

  npar = s.r.size()/3;

  s.omega = {0.,0.,0.,-0.5};
  s.com = {0.,0.05,0.4};
  s.vel = {0.,0.,-0.1};
  
  // write out initial particle base positions
  std::ofstream file_r("r.txt");
  std::ostream_iterator<float> fiter_r(file_r, "\n");
  std::copy(s.r.begin(), s.r.end(), fiter_r);
  
  
  // output files for rigid body position and orientation
  std::ofstream file_o("orient.txt");
  std::ofstream file_p("pos.txt");
  
  file_o << fmt::format("{:10.6f}{:10.6f}{:10.6f}{:10.6f}\n",
                         s.orient[0], s.orient[1], s.orient[2], s.orient[3]);          
  file_p << fmt::format("{:9.4f}{:9.4f}{:9.4f}\n",
                         s.com[0], s.com[1], s.com[2]);
  
  // initial impulse
  alpha = {0.,6.,-4,-1.};
  acc = {0.,0.,0.};
  evolve_motion(p, s, alpha, acc, p.dt);
  ncol = 0;
  
  
  // time loop  
  for (int k = 0; k < p.nt; k++){
    printf("%d\n", k);

    //evolve_motion_intertial(p, s1, p.dt);
    evolve_motion(p, s, {0.,0.,0.,0.}, {0.,0.,0.}, p.dt);
    
    transform_particles(p, s);
    
    rcount = 0;
    printf("%8.6f %8.6f %8.6f %8.6f %8.6f %8.6f\n", s.vel[0], s.vel[1], s.vel[2]);
    
    // file output
    file_o << fmt::format("{:10.6f}{:10.6f}{:10.6f}{:10.6f}\n",
                         s.orient[0], s.orient[1], s.orient[2], s.orient[3]);          
    file_p << fmt::format("{:12.4f}{:12.4f}{:12.4f}\n",
                         s.com[0], s.com[1], s.com[2]);
  }
  
  file_o.close();
  file_p.close();
    
  
  std::ofstream file_rt("rt.txt");
  std::ostream_iterator<float> fiter_rt(file_rt, "\n");
  std::copy(s.rt.begin(), s.rt.end(), fiter_rt);
    
  printf("%d\n", int(dis_vec.size()));
  
}
