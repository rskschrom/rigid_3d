#include <iostream>
#include <iomanip>
#include <fstream>
#include <math.h>
#include <vector>
#include "params.h"
#include "state.h"
#include "quat.h"
#include "util.h"
#include "geometry.h"

// evolve motion of body without accelerations
void evolve_motion_intertial(Params p, State &s, float dt){
  std::vector<float> dq(3);
  
  multiply(s.omega, s.orient, dq);

  // update body position and velocity
  for (int i = 0; i < 4; i++){
    s.orient[i] += dt*0.5*dq[i];
  }
  
  for (int i = 0; i < 3; i++){
    s.com[i] += dt*s.vel[i];
  }
  
  // update particle positions and velocities
  //transform_particles(p, s);
}

// evolve motion with RK4
void evolve_motion(Params p, State &s, std::vector<float> alpha, std::vector<float> acc, float dt){
  std::vector<float> k1(4), k2(4), k3(4), k4(4);
  std::vector<float> o2n(4), o3n(4), o4n(4);
  std::vector<float> q2n(4), q3n(4), q4n(4);
  
  // step 1
  multiply(s.omega, s.orient, k1);
  
  // step 2
  for (int i = 0; i < 4; i++){
    o2n[i] = s.omega[i]+alpha[i]*dt/2.;
    q2n[i] = s.orient[i]+k1[i]*dt/2.;
  }
  multiply(o2n, q2n, k2);
  
  // step 3
  for (int i = 0; i < 4; i++){
    o3n[i] = s.omega[i]+alpha[i]*dt/2.;
    q3n[i] = s.orient[i]+k2[i]*dt/2.;
  }
  multiply(o3n, q3n, k3);
  
  // step 4
  for (int i = 0; i < 4; i++){
    o4n[i] = s.omega[i]+alpha[i]*dt;
    q4n[i] = s.orient[i]+k3[i]*dt;
  }
  multiply(o4n, q4n, k4);

  // update orientation
  for (int i = 0; i < 4; i++){
    s.orient[i] += 1./6.*(k1[i]+2.*k2[i]+2.*k3[i]+k4[i])*dt;
    s.omega[i] += alpha[i]*dt;
  }
  
  // update velocity and center or mass
  for (int i = 0; i < 3; i++){
    s.com[i] += dt*s.vel[i];
    s.vel[i] += dt*acc[i];
  }
  
  // update particle positions and velocities
  //transform_particles(p, s);
}

// apply leonard jones potential to particles of each body
void lj_potential(Params p, State s1, State s2, std::vector<float> con,
                  std::vector<float> &f1, std::vector<float> &f2){
                  
  float rij_mag, rcon, damp;
  int npar1 = s1.r.size()/3;
  int npar2 = s2.r.size()/3;
    
  f1.resize(3*npar1);
  f2.resize(3*npar2);
    
  damp = 0.1;
    
  // force of body 1 particles on body 2 particles
  for (int i = 0; i < npar1; i++){
    f1[3*i] = 0.;
    f1[3*i+1] = 0.;
    f1[3*i+2] = 0.;
    
    // only calculate forces if near contact point
    rcon = pow(s1.rt[3*i]-con[0], 2)+
           pow(s1.rt[3*i+1]-con[1], 2)+
           pow(s1.rt[3*i+2]-con[2], 2);
           
    //if (rcon<9.*p.ds*p.ds){        
    for (int j = 0; j < npar2; j++){
      rij_mag = pow(s1.rt[3*i]-s2.rt[3*j], 2)+
                pow(s1.rt[3*i+1]-s2.rt[3*j+1], 2)+
                pow(s1.rt[3*i+2]-s2.rt[3*j+2], 2);

      if ((rij_mag>0.)&(rij_mag<4.*p.ds*p.ds)){
        rij_mag = std::min(rij_mag, p.ds*p.ds);
        f1[3*i] += -damp*p.g*(pow(p.ds, 4)/(rij_mag*rij_mag)-
                     p.ds*p.ds/rij_mag)*(s1.rt[3*i]-s2.rt[3*j])/rij_mag;
        f1[3*i+1] += -damp*p.g*(pow(p.ds, 4)/(rij_mag*rij_mag)-
                       p.ds*p.ds/rij_mag)*(s1.rt[3*i+1]-s2.rt[3*j+1])/rij_mag;
        f1[3*i+2] += -damp*p.g*(pow(p.ds, 4)/(rij_mag*rij_mag)-
                     p.ds*p.ds/rij_mag)*(s1.rt[3*i+2]-s2.rt[3*j+2])/rij_mag;
      }
    }
    //}
  }
    
  // force of body 2 particles on body 1 particles
  for (int i = 0; i < npar2; i++){
    f2[3*i] = 0.;
    f2[3*i+1] = 0.;
    f2[3*i+2] = 0.;
    
    // only calculate forces if near contact point
    rcon = pow(s2.rt[3*i]-con[0], 2)+
           pow(s2.rt[3*i+1]-con[1], 2)+
           pow(s2.rt[3*i+2]-con[2], 2);
           
    //if (rcon<9.*p.ds*p.ds){        
    for (int j = 0; j < npar1; j++){
      rij_mag = pow(s2.rt[3*i]-s1.rt[3*j], 2)+
                pow(s2.rt[3*i+1]-s1.rt[3*j+1], 2)+
                pow(s2.rt[3*i+2]-s1.rt[3*j+2], 2);

      if ((rij_mag>0.)&(rij_mag<4.*p.ds*p.ds)){
        rij_mag = std::min(rij_mag, p.ds*p.ds);
        f2[3*i] += -damp*p.g*(pow(p.ds, 4)/(rij_mag*rij_mag)-
                   p.ds*p.ds/rij_mag)*(s2.rt[3*i]-s1.rt[3*j])/rij_mag;
        f2[3*i+1] += -damp*p.g*(pow(p.ds, 4)/(rij_mag*rij_mag)-
                     p.ds*p.ds/rij_mag)*(s2.rt[3*i+1]-s1.rt[3*j+1])/rij_mag;
        f2[3*i+2] += -damp*p.g*(pow(p.ds, 4)/(rij_mag*rij_mag)-
                     p.ds*p.ds/rij_mag)*(s2.rt[3*i+2]-s1.rt[3*j+2])/rij_mag;
      }
    }
    //}
  }
}


// use leonard jones potential to simulate collision impulse
void lj_collision(Params p, State &s1, State &s2, std::vector<float> con, std::vector<float> &dvel1, std::vector<float> &dvel2){
  int i;
  float rij_mag, f1_mag, f2_mag, imp, vrel, pair_dist;
  float f1_mag_sum, f2_mag_sum, tcol, denom1, denom2;
  float xcc1, ycc1, zcc1, xcc2, ycc2, zcc2;
  
  int npar1 = s1.r.size()/3;
  int npar2 = s2.r.size()/3;
  std::vector<float> f1, f2;
  std::vector<float> nvec(3), f1_sum(3), f2_sum(3), v1(3), v2(3), b1(3), b2(3);
  std::vector<float> iitens_rot1(9), iitens_rot2(9);
  
  // transform inertia moment tensors
  rotate_itens(p, s1, iitens_rot1);
  rotate_itens(p, s2, iitens_rot2);

  // get potential forces to use for surface normal vectors
  lj_potential(p, s1, s2, con, f1, f2);
    
  for (i = 0; i < npar1; i++){
    f1_sum[0] += f1[3*i];
    f1_sum[1] += f1[3*i+1];
    f1_sum[2] += f1[3*i+2];
  }
  
  f1_mag = sqrt(f1_sum[0]*f1_sum[0]+f1_sum[1]*f1_sum[1]+f1_sum[2]*f1_sum[2]);
  //printf("f1mag %8.6f\n", f1_mag);
  
  // avoid divide by zero
  if (f1_mag>0.){
    nvec[0] = f1_sum[0]/f1_mag;
    nvec[1] = f1_sum[1]/f1_mag;
    nvec[2] = f1_sum[2]/f1_mag;
  }
  else{
    nvec[0] = 0.;
    nvec[1] = 0.;
    nvec[2] = 0.;
  }
  
  //printf("%8.6f %8.6f %8.6f\n", nvec[0], nvec[1], nvec[2]);

  // calculate impulse
  xcc1 = con[0]-s1.com[0];
  ycc1 = con[1]-s1.com[1];
  zcc1 = con[2]-s1.com[2];
  
  xcc2 = con[0]-s2.com[0];
  ycc2 = con[1]-s2.com[1];
  zcc2 = con[2]-s2.com[2];
  
  v1[0] = s1.vel[0]+s1.omega[2]*zcc1-s1.omega[3]*ycc1;
  v1[1] = s1.vel[1]+s1.omega[3]*xcc1-s1.omega[1]*zcc1;
  v1[2] = s1.vel[2]+s1.omega[1]*ycc1-s1.omega[2]*xcc1;
  
  v2[0] = s2.vel[0]+s2.omega[2]*zcc2-s2.omega[3]*ycc2;
  v2[1] = s2.vel[1]+s2.omega[3]*xcc2-s2.omega[1]*zcc2;
  v2[2] = s2.vel[2]+s2.omega[1]*ycc2-s2.omega[2]*xcc2;

  vrel = nvec[0]*(v1[0]-v2[0])+nvec[1]*(v1[1]-v2[1])+nvec[2]*(v1[2]-v2[2]);
  
  //printf("%8.6f %8.6f %8.6f %8.6f\n", vrel, s1.vel[0], s1.vel[1], s1.vel[2]);
  b1[0] = iitens_rot1[0]*(ycc1*nvec[2]-zcc1*nvec[1])+
          iitens_rot1[1]*(zcc1*nvec[0]-xcc1*nvec[2])+
          iitens_rot1[2]*(xcc1*nvec[1]-ycc1*nvec[0]);
  b1[1] = iitens_rot1[3]*(ycc1*nvec[2]-zcc1*nvec[1])+
          iitens_rot1[4]*(zcc1*nvec[0]-xcc1*nvec[2])+
          iitens_rot1[5]*(xcc1*nvec[1]-ycc1*nvec[0]);
  b1[2] = iitens_rot1[6]*(ycc1*nvec[2]-zcc1*nvec[1])+
          iitens_rot1[7]*(zcc1*nvec[0]-xcc1*nvec[2])+
          iitens_rot1[8]*(xcc1*nvec[1]-ycc1*nvec[0]);
          
  b2[0] = iitens_rot2[0]*(ycc2*nvec[2]-zcc2*nvec[1])+
          iitens_rot2[1]*(zcc2*nvec[0]-xcc2*nvec[2])+
          iitens_rot2[2]*(xcc2*nvec[1]-ycc2*nvec[0]);
  b2[1] = iitens_rot2[3]*(ycc2*nvec[2]-zcc2*nvec[1])+
          iitens_rot2[4]*(zcc2*nvec[0]-xcc2*nvec[2])+
          iitens_rot2[5]*(xcc2*nvec[1]-ycc2*nvec[0]);
  b2[2] = iitens_rot2[6]*(ycc2*nvec[2]-zcc2*nvec[1])+
          iitens_rot2[7]*(zcc2*nvec[0]-xcc2*nvec[2])+
          iitens_rot2[8]*(xcc2*nvec[1]-ycc2*nvec[0]);
  
  denom1 = nvec[0]*(b1[1]*zcc1-b1[2]*ycc1)+
           nvec[1]*(b1[2]*xcc1-b1[0]*zcc1)+
           nvec[2]*(b1[0]*ycc1-b1[1]*xcc1);
           
  denom2 = nvec[0]*(b2[1]*zcc2-b2[2]*ycc2)+
           nvec[1]*(b2[2]*xcc2-b2[0]*zcc2)+
           nvec[2]*(b2[0]*ycc2-b2[1]*xcc2);
  imp = -(1.+p.eps)*vrel/(1./s1.mass+1./s2.mass+denom1+denom2);
  //imp = -(1.+p.eps)*vrel/(1./s1.mass+1./s2.mass);
    
  //printf("collision impulse %8.6f, vrel %8.6f, nvec %8.6f %8.6f %8.6f\n", imp, vrel, nvec[0], nvec[1], nvec[2]);
  
  // change in linear velocities from collision
  dvel1[0] = imp*nvec[0]/s1.mass;
  dvel1[1] = imp*nvec[1]/s1.mass;
  dvel1[2] = imp*nvec[2]/s1.mass;
  
  dvel2[0] = -imp*nvec[0]/s2.mass;
  dvel2[1] = -imp*nvec[1]/s2.mass;
  dvel2[2] = -imp*nvec[2]/s2.mass;
  
  // change in angular velocities from collision
  dvel1[3] = imp*b1[0];
  dvel1[4] = imp*b1[1];
  dvel1[5] = imp*b1[2];
  
  dvel2[3] = -imp*b2[0];
  dvel2[4] = -imp*b2[1];
  dvel2[5] = -imp*b2[2];
  //printf("%8.6f %8.6f %8.6f\n", iitens_rot2[0], ycc2, nvec[2]);

}

// summed result of multiple collisions
void multi_collision(Params p, State s1, State s2, std::vector<float> alpha1, std::vector<float> acc1,
                     std::vector<float> alpha2, std::vector<float> acc2,
                     std::vector<float> &dv1_tot, std::vector<float> &dv2_tot){
  int ncol, ncol_act;  
  int npar1 = s1.r.size()/3;
  int npar2 = s2.r.size()/3;
  std::vector<float> dis_vec;
  std::vector<float> dv1(6), dv2(6), con(3);
  std::vector<int> ip1, ip2;
  
  // clear prior dv
  dv1_tot.clear();
  dv2_tot.clear();
  
  // determine if collision has occured
  inter_particle_pair(p, s1, s2, ip1, ip2, dis_vec);
  ncol = dis_vec.size();
  ncol_act = 0;
  
  for (int k = 0; k < ncol; k++){       
    con[0] = 0.5*(s1.rt[3*ip1[k]]+s2.rt[3*ip2[k]]);
    con[1] = 0.5*(s1.rt[3*ip1[k]+1]+s2.rt[3*ip2[k]+1]);
    con[2] = 0.5*(s1.rt[3*ip1[k]+2]+s2.rt[3*ip2[k]+2]);
    
    // get velocity adjustments to first collision
    lj_collision(p, s1, s2, con, dv1, dv2);
    
    if ((dv1[0]*dv1[0]+dv1[1]*dv1[1]+dv1[2]*dv1[2])>0.){
      ncol_act++;
    }
    
    for (int i = 0; i < 6; i++){
      dv1_tot[i] += dv1[i];
      dv2_tot[i] += dv2[i];
    }
  }
  
  // scale dv by number of collisions
  for (int i = 0; i < 6; i++){
    dv1_tot[i] = dv1_tot[i]/ncol_act;
    dv2_tot[i] = dv2_tot[i]/ncol_act;
  }
}

// get first collision (or mean of multiple collisions on plane)
void first_collision(Params p, State s1, State s2, std::vector<float> alpha1, std::vector<float> acc1,
                     std::vector<float> alpha2, std::vector<float> acc2,
                     std::vector<float> &dv1_tot, std::vector<float> &dv2_tot){
  int ncol, ncol_old, niter, niter_max;
  int npar1 = s1.r.size()/3;
  int npar2 = s2.r.size()/3;
  std::vector<float> dis_vec;
  std::vector<float> dv1(6), dv2(6), con(3);
  std::vector<int> ip1, ip2;
  float wgt, wgt_sum, dt_step, dt_old;
  
  // clear prior dv
  dv1_tot.clear();
  dv2_tot.clear();
  
  // determine if collision has occured
  inter_particle_pair(p, s1, s2, ip1, ip2, dis_vec);
  ncol = dis_vec.size();
  
  // get more accurate order of collisions with microstepping
  niter = 0;
  niter_max = 10;
  dt_step = p.dt/(niter_max-1);
  ncol_old = ncol;
  
  if (ncol > 0){
  
    while ((niter < niter_max)&&(ncol > 0)){
    
      evolve_motion(p, s1, alpha1, acc1, -dt_step);
      evolve_motion(p, s2, alpha2, acc2, -dt_step);
      
      transform_particles(p, s1);
      transform_particles(p, s2);
      
      inter_particle_pair(p, s1, s2, ip1, ip2, dis_vec);
      ncol = dis_vec.size();
      //dt_old = dt_step;
      
      /*           
      if ((ncol == 0)&&(ncol != ncol_old)){
        dt_step = -dt_step/2.;
        niter++;
      }
      
      if (ncol > ncol_old){
        dt_step = -dt_step/2.;
        niter++;
      }
      */
      niter++;
      
      //printf("ncol test %d %d\n", niter, ncol);
      
      ncol_old = ncol;
      
    }
    
    if (ncol==0){
      evolve_motion(p, s1, alpha1, acc1, dt_step);
      evolve_motion(p, s2, alpha2, acc2, dt_step);
      
      transform_particles(p, s1);
      transform_particles(p, s2);
      
      inter_particle_pair(p, s1, s2, ip1, ip2, dis_vec);
      ncol = dis_vec.size();
    }
    
    // get average impulse over all collusions
    wgt_sum = 0.;
    printf("ncol %d\n", ncol);
  
    for (int i = 0; i < ncol; i++){
      con[0] = 0.5*(s1.rt[3*ip1[i]]+s2.rt[3*ip2[i]]);
      con[1] = 0.5*(s1.rt[3*ip1[i]+1]+s2.rt[3*ip2[i]+1]);
      con[2] = 0.5*(s1.rt[3*ip1[i]+2]+s2.rt[3*ip2[i]+2]);
      //printf("%8.6f %8.6f %8.6f\n", con[0], con[1], con[2]);
    
      // get velocity adjustments to first collision
      lj_collision(p, s1, s2, con, dv1, dv2);
    
      for (int i = 0; i < 6; i++){
        dv1_tot[i] += dv1[i]/ncol;
        dv2_tot[i] += dv2[i]/ncol;
      }
    }
  }
}

// single collisions
void single_collision(Params p, State s1, State s2, std::vector<float> alpha1, std::vector<float> acc1,
                     std::vector<float> alpha2, std::vector<float> acc2,
                     std::vector<float> &dv1_tot, std::vector<float> &dv2_tot){
  int ncol;
  int npar1 = s1.r.size()/3;
  int npar2 = s2.r.size()/3;
  std::vector<float> dis_vec;
  std::vector<float> dv1(6), dv2(6), con(3);
  std::vector<int> ip1, ip2;
  float wgt, wgt_sum;
  
  // clear prior dv
  dv1_tot.clear();
  dv2_tot.clear();
  
  // determine if collision has occured
  inter_particle_pair(p, s1, s2, ip1, ip2, dis_vec);
  ncol = dis_vec.size();
  
  // weighted average of contact points
  wgt_sum = 0.;
  
  if (ncol > 0){       
    for (int i = 0; i < ncol; i++){
      wgt = 1./pow(dis_vec[i], 2);
      //printf("%8.6f ", wgt);
      wgt_sum += wgt;
      con[0] += 0.5*wgt*(s1.rt[3*ip1[i]]+s2.rt[3*ip2[i]]);
      con[1] += 0.5*wgt*(s1.rt[3*ip1[i]+1]+s2.rt[3*ip2[i]+1]);
      con[2] += 0.5*wgt*(s1.rt[3*ip1[i]+2]+s2.rt[3*ip2[i]+2]);
    }
    //printf("\n");
    
    con[0] = con[0]/wgt_sum;
    con[1] = con[1]/wgt_sum;
    con[2] = con[2]/wgt_sum;
    
    printf("%8.6f %8.6f %8.6f\n", con[0], con[1], con[2]);
    
    // get velocity adjustments to first collision
    lj_collision(p, s1, s2, con, dv1, dv2);
    
    for (int i = 0; i < 6; i++){
      dv1_tot[i] += dv1[i];
      dv2_tot[i] += dv2[i];
    }
  }
}

// test collision method
void test_collision(Params p, State &s1, State &s2, std::vector<float> alpha1, std::vector<float> acc1,
                    std::vector<float> alpha2, std::vector<float> acc2, bool &attached, int &rcount){
  int ncol, step, nmulti;
  float dt_fine;
  
  int npar1 = s1.r.size()/3;
  int npar2 = s2.r.size()/3;
  std::vector<float> dis_vec;
  std::vector<float> dv1(6), dv2(6), dv1_tot(6), dv2_tot(6), con(3);
  std::vector<int> ip1, ip2;
  bool intersect = true;
  
  rcount++;
  
  // determine if collision has occured    
  transform_particles(p, s1);
  transform_particles(p, s2);
  
  inter_particle_pair(p, s1, s2, ip1, ip2, dis_vec);
  ncol = dis_vec.size();
  printf("collision points %d\n", ncol);
  
  //multi_collision(p, s1, s2, alpha1, acc1, alpha2, acc2, dv1_tot, dv2_tot);
  first_collision(p, s1, s2, alpha1, acc1, alpha2, acc2, dv1_tot, dv2_tot);
  //single_collision(p, s1, s2, alpha1, acc1, alpha2, acc2, dv1_tot, dv2_tot);
  
  // carry out collision via microsteps
  if (ncol>0){
    dt_fine = p.dt/10.;
    step = 0;
    nmulti = 0;
    //attached = true;
    /*
    // back integration to collision time
    while ((intersect)&&(step < 20)){
      step++;
      
            
    }
    */
    // motion adjustment post collision
    for (int i = 0; i < 3; i++){
      s1.vel[i] += dv1_tot[i];
      s2.vel[i] += dv2_tot[i];
      s1.omega[i+1] += dv1_tot[i+3];
      s2.omega[i+1] += dv2_tot[i+3];
    }
    
    printf("%8.6f %8.6f\n", s1.vel[2], s2.vel[2]);
    /*
    // test if another collision will occur after adjustment
    evolve_motion(p, s1, alpha1, acc1, p.dt);
    evolve_motion(p, s2, alpha2, acc2, p.dt);
      
    transform_particles(p, s1);
    transform_particles(p, s2);
  
    inter_particle_pair(p, s1, s2, ip1, ip2, dis_vec);
    ncol = dis_vec.size();
    printf("post collision %d %d\n", ncol, rcount);
    
    if ((ncol > 0)&&(rcount < 20)){
      evolve_motion(p, s1, alpha1, acc1, -p.dt);
      evolve_motion(p, s2, alpha2, acc2, -p.dt);
      
      transform_particles(p, s1);
      transform_particles(p, s2);
      test_collision(p, s1, s2, alpha1, acc1, alpha2, acc2, attached, rcount);
    }
    */
    
    //evolve_motion(p, s1, alpha1, acc1, -p.dt);
    //evolve_motion(p, s2, alpha2, acc2, -p.dt);
 
  }
  
  
}
