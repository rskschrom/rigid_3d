#include <iostream>
#include <iomanip>
#include <fstream>
#include <math.h>
#include <vector>
#include "params.h"
#include "state.h"
#include "geometry.h"

// rotate and translate body particles
void transform_particles(Params p, State &s){
  int i;
  int npar = s.x.size();
  
  for (i = 0; i < npar; i++){
      s.xt[i] = cos(s.theta)*s.x[i]-sin(s.theta)*s.y[i]+s.xcom;
      s.yt[i] = sin(s.theta)*s.x[i]+cos(s.theta)*s.y[i]+s.ycom;
  }
}

// evolve motion of body
void evolve_motion(Params p, State &s, float alpha, float ax, float ay, float dt){
  // update body position and velocity
  s.theta = s.theta+dt*s.omega;
  s.omega = s.omega+dt*alpha;
  s.xcom = s.xcom+dt*s.vx;
  s.ycom = s.ycom+dt*s.vy;
  s.vx = s.vx+dt*ax;
  s.vy = s.vy+dt*ay;
  
  // update particle positions and velocities
  transform_particles(p, s);
}

// apply leonard jones potential to particles of each body
void lj_potential(Params p, State s1, State s2, float xcon, float ycon,
                  std::vector<float> &fx1, std::vector<float> &fy1,
                  std::vector<float> &fx2, std::vector<float> &fy2){
  int i, j;
  float rij_mag, rcon, damp;
  int npar1 = s1.x.size();
  int npar2 = s2.x.size();
    
  fx1.resize(npar1);
  fy1.resize(npar1);
  fx2.resize(npar2);
  fy2.resize(npar2);
    
  damp = 0.1;
    
  // force of body 1 particles on body 2 particles
  for (i = 0; i < npar1; i++){
    fx1[i] = 0.;
    fy1[i] = 0.;
    
    // only calculate forces if near contact point
    rcon = sqrt(pow(s1.xt[i]-xcon, 2)+pow(s1.yt[i]-ycon, 2));
    if (rcon<3.*p.ds){        
      for (j = 0; j < npar2; j++){
        rij_mag = sqrt(pow(s1.xt[i]-s2.xt[j], 2)+pow(s1.yt[i]-s2.yt[j], 2));

        if ((rij_mag>0.)&(rij_mag<p.ds)){
                
          fx1[i] += -damp*p.g*(pow(p.ds/rij_mag, 4)-pow(p.ds/rij_mag, 2))*(s1.xt[i]-s2.xt[j])/pow(rij_mag, 2);
          fy1[i] += -damp*p.g*(pow(p.ds/rij_mag, 4)-pow(p.ds/rij_mag, 2))*(s1.yt[i]-s2.yt[j])/pow(rij_mag, 2);
        }
      }
    }
  }
    
  // force of body 2 particles on body 1 particles
  for (i = 0; i < npar2; i++){
    fx2[i] = 0.;
    fy2[i] = 0.;
    
    // only calculate forces if near contact point
    rcon = sqrt(pow(s2.xt[i]-xcon, 2)+pow(s2.yt[i]-ycon, 2));
    if (rcon<3.*p.ds){     
      for (j = 0; j < npar1; j++){
        rij_mag = sqrt(pow(s2.xt[i]-s1.xt[j], 2)+pow(s2.yt[i]-s1.yt[j], 2));

        if ((rij_mag>0.)&(rij_mag<p.ds)){
                
          fx2[i] += -damp*p.g*(pow(p.ds/rij_mag, 4)-pow(p.ds/rij_mag, 2))*(s2.xt[i]-s1.xt[j])/pow(rij_mag, 2);
          fy2[i] += -damp*p.g*(pow(p.ds/rij_mag, 4)-pow(p.ds/rij_mag, 2))*(s2.yt[i]-s1.yt[j])/pow(rij_mag, 2);
        }
      }
    }
  }
}

// calculate first collision
void first_collision(Params p, State &s1, State &s2, std::vector<int> &ip1, std::vector<int> &ip2, float tcol){
  int i, j, ip1_first, ip2_first;
  int ncol = ip1.size();
  std::vector<float> tcols(ncol);
  float x1, y1, x2, y2, vrx, vry, r12;
  float vx1, vy1, vx2, vy2, e12x, e12y, vr12;
  
  // get collision times
  //printf("%d\n", ncol);
  tcol = 0.1;
  ip1_first = ip1[0];
  ip2_first = ip2[0];
  
  for (i = 0; i < ncol; i++){
    x1 = s1.xt[ip1[i]];
    x2 = s2.xt[ip2[i]];
    y1 = s1.yt[ip1[i]];
    y2 = s2.yt[ip2[i]];
    
    // relative velocities
    vx1 = s1.vx-s1.omega*(y1-s1.ycom);
    vx2 = s2.vx-s2.omega*(y2-s2.ycom);
    vy1 = s1.vy+s1.omega*(x1-s1.xcom);
    vy2 = s2.vy+s2.omega*(x2-s2.xcom);
    vrx = vx1-vx2;
    vry = vy1-vy2;
    
    
    // velocity along vector direction from p1-p2
    r12 = sqrt(pow(x2-x1, 2)+pow(y2-y1, 2));
    e12x = (x1-x2)/r12;
    e12y = (y1-y2)/r12;
    vr12 = vrx*e12x+vry*e12y;
    tcols[i] = (p.ds-r12)/vr12;
    
    //printf("%d %8.6f %8.6f\n", i, r12, vry);
    //printf("%8.6f\n", s1.vy-s2.vy);
    //printf("tcol %8.6f\n", tcols[i]);
    
    // compare collision time to earliest collision time
    if (tcols[i]<tcol){
      tcol = tcols[i];
      ip1_first = ip1[i];
      ip2_first = ip2[i];
    } 
    
  }
  
  // return 1-length vectors of particle indices for first collision
  ip1.resize(1);
  ip2.resize(1);
  ip1[0] = ip1_first;
  ip2[0] = ip2_first;
}

// use leonard jones potential to simulate collision impulse
void lj_collision(Params p, State &s1, State &s2, float alpha1, float alpha2, float ax1, float ax2, float ay1, float ay2){
  int i, j, cpi1, cpi2;
  float rij_mag, fx1_sum, fy1_sum, fx2_sum, fy2_sum, f1_mag, f2_mag;
  float nx, ny, imp, vrel, pair_dist;
  float vx1, vy1, vx2, vy2, eps, xcon1, ycon1, xcon2, ycon2, xcon, ycon;
  float f1_mag_sum, f2_mag_sum, tcol;
  
  int npar1 = s1.x.size();
  int npar2 = s2.x.size();
  std::vector<float> fx1, fy1, fx2, fy2, dis_vec;
  std::vector<int> ip1, ip2;
  
  // determine if collision has occured
  inter_particle_pair(s1.xt, s1.yt, s2.xt, s2.yt, p.ds, ip1, ip2, dis_vec);
  //printf("%d\n", int(dis_vec.size()));  
  if (dis_vec.size()>0){
    //printf("%d contact(s)!, %8.6f\n", int(dis_vec.size()), dis_vec[0]);
    
    // find first collision time and location in timestep
    first_collision(p, s1, s2, ip1, ip2, tcol);
    xcon1 = s1.xt[ip1[0]];
    ycon1 = s1.yt[ip1[0]];
    xcon2 = s2.xt[ip2[0]];
    ycon2 = s2.yt[ip2[0]];
    xcon = 0.5*(xcon1+xcon2);
    ycon = 0.5*(ycon1+ycon2);
    
    lj_potential(p, s1, s2, xcon, ycon, fx1, fy1, fx2, fy2);
    
    // backward integration to initial collision time
    evolve_motion(p, s1, alpha1, ax1, ay1, -p.dt);
    evolve_motion(p, s2, alpha2, ax2, ay2, -p.dt);

     
    // get potential forces to use for surface normal vectors
    fx1_sum = 0.;
    fy1_sum = 0.;
    f1_mag_sum = 0.;
    
    for (i = 0; i < npar1; i++){
      fx1_sum += fx1[i];
      fy1_sum += fy1[i];
      f1_mag_sum += sqrt(pow(fx1[i], 2)+pow(fy1[i], 2));
    }
  
    f1_mag = sqrt(fx1_sum*fx1_sum+fy1_sum*fy1_sum);
    nx = fx1_sum/f1_mag;
    ny = fy1_sum/f1_mag;

    // calculate impulse
    eps = 0.1;
    vx1 = s1.vx-s1.omega*(ycon-s1.ycom);
    vx2 = s2.vx-s2.omega*(ycon-s2.ycom);
    vy1 = s1.vy+s1.omega*(xcon-s1.xcom);
    vy2 = s2.vy+s2.omega*(xcon-s2.xcom);
    vrel = nx*(vx1-vx2)+ny*(vy1-vy2);
    
    imp = -(1.+eps)*vrel/(1./s1.mass+1./s2.mass+pow(nx*(ycon-s1.ycom)-ny*(xcon-s1.xcom), 2)/s1.inerm+
                                                pow(nx*(ycon-s2.ycom)-ny*(xcon-s2.xcom), 2)/s2.inerm);
    

    // adjust angular and linear velocities to account for collision
    s1.vx += imp*nx/s1.mass;
    s1.vy += imp*ny/s1.mass;
    s2.vx += -imp*nx/s2.mass;
    s2.vy += -imp*ny/s2.mass;
    
    s1.omega += -imp*(nx*(ycon-s1.ycom)-ny*(xcon-s1.xcom))/s1.inerm;
    s2.omega += imp*(nx*(ycon-s2.ycom)-ny*(xcon-s2.xcom))/s2.inerm;

  }
}
