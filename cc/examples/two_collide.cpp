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

  int npar1, npar2, ncol, rcount;
  std::vector<float> alpha, acc, dis_vec, f1, f2, bbox1(6), bbox2(6), con(3);
  std::vector<float> dv1(6), dv2(6), iitens_rot(9), traj1, traj2;
  std::vector<int> ip1, ip2;
  
  bool attached = false;

  // initial body properties
  struct Params p;
  struct State s1, s2;

  s1.rho = 920.;
  s2.rho = 920.;
  
  // define particle bodies and transform them
  init_body_file(p, s1, "crystal_0000_1.5.txt");
  init_body_file(p, s2, "crystal_0000_1.5.txt");

  npar1 = s1.r.size()/3;
  npar2 = s2.r.size()/3;

  s1.omega = {0.,0.,0.,-0.5};
  s2.omega = {0.,0.,0.,0.1};
  s1.com = {0.,0.05,0.4};
  s2.com = {0.,0.,0.-0.2};
  //s1.orient = {sqrt(2.)/2.,0.,sqrt(2.)/2.,0.};
  s1.vel = {0.,0.,-0.1};
  
  // write out initial particle base positions
  std::ofstream file_r1("r1.txt");
  std::ostream_iterator<float> fiter_r1(file_r1, "\n");
  std::copy(s1.r.begin(), s1.r.end(), fiter_r1);
  
  std::ofstream file_r2("r2.txt");
  std::ostream_iterator<float> fiter_r2(file_r2, "\n");
  std::copy(s2.r.begin(), s2.r.end(), fiter_r2);
  
  // output files for rigid body position and orientation
  std::ofstream file_o1("orient1.txt");
  std::ofstream file_p1("pos1.txt");
  std::ofstream file_o2("orient2.txt");
  std::ofstream file_p2("pos2.txt");
  
  file_o1 << fmt::format("{:10.6f}{:10.6f}{:10.6f}{:10.6f}\n",
                         s1.orient[0], s1.orient[1], s1.orient[2], s1.orient[3]);          
  file_p1 << fmt::format("{:9.4f}{:9.4f}{:9.4f}\n",
                         s1.com[0], s1.com[1], s1.com[2]);
                         
  file_o2 << fmt::format("{:10.6f}{:10.6f}{:10.6f}{:10.6f}\n",
                         s2.orient[0], s2.orient[1], s2.orient[2], s2.orient[3]);          
  file_p2 << fmt::format("{:9.4f}{:9.4f}{:9.4f}\n",
                         s2.com[0], s2.com[1], s2.com[2]);
  

  // initial impulse
  alpha = {0.,6.,-4,-1.};
  //alpha = {0.,0.,0.,0.};
  acc = {0.,0.,0.};
  evolve_motion(p, s1, alpha, acc, p.dt);
  ncol = 0;
  
  /*
  // test particle trajectory
  transform_particles(p, s1);
  transform_particles(p, s2);
    
  traj1.resize(3*2000);
  traj2.resize(3*2000);
  particle_trajectory(p, s1, 0, 0.01, traj1);
  particle_trajectory(p, s2, 0, 0.01, traj2);
  
  std::ofstream file_tr1("traj1.txt");
  std::ostream_iterator<float> fiter_tr1(file_tr1, "\n");
  std::copy(traj1.begin(), traj1.end(), fiter_tr1);
  
  std::ofstream file_tr2("traj2.txt");
  std::ostream_iterator<float> fiter_tr2(file_tr2, "\n");
  std::copy(traj2.begin(), traj2.end(), fiter_tr2);
  */
  
  // time loop  
  for (int k = 0; k < p.nt; k++){
    printf("%d\n", k);

    //evolve_motion_intertial(p, s1, p.dt);
    evolve_motion(p, s1, {0.,0.,0.,0.}, {0.,0.,0.}, p.dt);
    evolve_motion(p, s2, {0.,0.,0.,0.}, {0.,0.,0.}, p.dt);
    
    transform_particles(p, s1);
    transform_particles(p, s2);
    
    //body_bbox(p, s1, bbox1);
    //body_bbox(p, s2, bbox2);
    rcount = 0;
    //if (!attached){
    test_collision(p, s1, s2, {0.,0.,0.,0.}, {0.,0.,0.}, {0.,0.,0.,0.}, {0.,0.,0.},
                   attached, rcount);
    //}
                   
    printf("%8.6f %8.6f %8.6f %8.6f %8.6f %8.6f\n", s1.vel[0], s1.vel[1], s1.vel[2], s2.vel[0], s2.vel[1], s2.vel[2]);
    
    // file output
    file_o1 << fmt::format("{:10.6f}{:10.6f}{:10.6f}{:10.6f}\n",
                         s1.orient[0], s1.orient[1], s1.orient[2], s1.orient[3]);          
    file_p1 << fmt::format("{:12.4f}{:12.4f}{:12.4f}\n",
                         s1.com[0], s1.com[1], s1.com[2]);
                         
    file_o2 << fmt::format("{:10.6f}{:10.6f}{:10.6f}{:10.6f}\n",
                         s2.orient[0], s2.orient[1], s2.orient[2], s2.orient[3]);          
    file_p2 << fmt::format("{:12.4f}{:12.4f}{:12.4f}\n",
                         s2.com[0], s2.com[1], s2.com[2]);

  }
  
  file_o1.close();
  file_p1.close();
  file_o2.close();
  file_p2.close();
    
  
  std::ofstream file_rt1("rt1.txt");
  std::ostream_iterator<float> fiter_rt1(file_rt1, "\n");
  std::copy(s1.rt.begin(), s1.rt.end(), fiter_rt1);
    
  std::ofstream file_rt2("rt2.txt");
  std::ostream_iterator<float> fiter_rt2(file_rt2, "\n");
  std::copy(s2.rt.begin(), s2.rt.end(), fiter_rt2);
    
  //std::ofstream file_f1("f1.txt");
  //std::ostream_iterator<float> fiter_f1(file_f1, "\n");
  //std::copy(f1.begin(), f1.end(), fiter_f1);
    
  //std::ofstream file_f2("f2.txt");
  //std::ostream_iterator<float> fiter_f2(file_f2, "\n");
  //std::copy(f2.begin(), f2.end(), fiter_f2);
    
  printf("%d\n", int(dis_vec.size()));
  
}
