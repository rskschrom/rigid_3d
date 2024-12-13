#include <iostream>
#include <cmath>
#include <vector>
#include "params.h"
#include "state.h"
//#include "util.h"
#include "quat.h"
#include <Eigen/Dense>

// calculate bounding box of body
void body_bbox(Params p, State s, std::vector<float> &bbox){
  float xmin, xmax, ymin, ymax, zmin, zmax;
  int npar = s.rt.size()/3;
  
  xmin = s.rt[0];
  xmax = s.rt[0];
  ymin = s.rt[1];
  ymax = s.rt[1];
  zmin = s.rt[2];
  zmax = s.rt[2];
  
  for (int i = 0; i < npar; i++){
    if (s.rt[3*i]<xmin){
      xmin = s.rt[3*i];
    }
    if (s.rt[3*i]>xmax){
      xmax = s.rt[3*i];
    }
    if (s.rt[3*i+1]<ymin){
      ymin = s.rt[3*i+1];
    }
    if (s.rt[3*i+1]>ymax){
      ymax = s.rt[3*i+1];
    }
    if (s.rt[3*i+2]<zmin){
      zmin = s.rt[3*i+2];
    }
    if (s.rt[3*i+2]>zmax){
      zmax = s.rt[3*i+2];
    }
  }
  
  bbox[0] = xmin;
  bbox[1] = xmax;
  bbox[2] = ymin;
  bbox[3] = ymax;
  bbox[4] = zmin;
  bbox[5] = zmax;
  
}

// rotate and translate body particles
void transform_particles(Params p, State &s){
  int npar = s.r.size()/3;
  float a, b, c;
  std::vector<float> prot(3);

  a = s.orient[0]*s.orient[0]-
      s.orient[1]*s.orient[1]-
      s.orient[2]*s.orient[2]-
      s.orient[3]*s.orient[3];
      
  c = -2.*s.orient[0];
  
  for (int i = 0; i < npar; i++){
    b = 2.*(s.r[3*i]*s.orient[1]+
            s.r[3*i+1]*s.orient[2]+
            s.r[3*i+2]*s.orient[3]);
            
    // prot = a*p+b*v+c*pxv
    prot[0] = a*s.r[3*i]+b*s.orient[1]+
              c*(s.r[3*i+1]*s.orient[3]-s.r[3*i+2]*s.orient[2]);
              
    prot[1] = a*s.r[3*i+1]+b*s.orient[2]+
              c*(s.r[3*i+2]*s.orient[1]-s.r[3*i]*s.orient[3]);
              
    prot[2] = a*s.r[3*i+2]+b*s.orient[3]+
              c*(s.r[3*i]*s.orient[2]-s.r[3*i+1]*s.orient[1]);
    
    s.rt[3*i] = prot[0]+s.com[0];
    s.rt[3*i+1] = prot[1]+s.com[1];
    s.rt[3*i+2] = prot[2]+s.com[2];
  }
}

// get particle trajectory
void particle_trajectory(Params p, State s, int ip, float dt, std::vector<float> &traj){
  int ntime = traj.size()/3;
  float a, b, c, omega_mag, xrel, yrel, zrel, ang;
  std::vector<float> prot(3), urot(3), qrot(4);
  
  // rotation axis
  omega_mag = sqrt(s.omega[1]*s.omega[1]+
                   s.omega[2]*s.omega[2]+
                   s.omega[3]*s.omega[3]);
  urot = {s.omega[1]/omega_mag,s.omega[2]/omega_mag,s.omega[3]/omega_mag};
  
  xrel = s.rt[3*ip]-s.com[0];
  yrel = s.rt[3*ip+1]-s.com[1];
  zrel = s.rt[3*ip+2]-s.com[2];
  
  // time loop
  for (int i = 0; i < ntime; i++){
    // rotational component
    ang = i*dt*omega_mag;
    qrot = {cos(float(ang/2.)),
            sin(float(ang/2.))*urot[0],
            sin(float(ang/2.))*urot[1],
            sin(float(ang/2.))*urot[2]};
  
    a = qrot[0]*qrot[0]-
        qrot[1]*qrot[1]-
        qrot[2]*qrot[2]-
        qrot[3]*qrot[3];
      
    c = -2.*qrot[0];
  
  
    b = 2.*(xrel*qrot[1]+
            yrel*qrot[2]+
            zrel*qrot[3]);
            
    // prot = a*p+b*v+c*pxv
    prot[0] = a*xrel+b*qrot[1]+c*(yrel*qrot[3]-zrel*qrot[2]);
              
    prot[1] = a*yrel+b*qrot[2]+c*(zrel*qrot[1]-xrel*qrot[3]);
              
    prot[2] = a*zrel+b*qrot[3]+c*(xrel*qrot[2]-yrel*qrot[1]);
  
    traj[3*i] = prot[0]+s.vel[0]*i*dt+s.com[0];
    traj[3*i+1] = prot[1]+s.vel[1]*i*dt+s.com[1];
    traj[3*i+2] = prot[2]+s.vel[2]*i*dt+s.com[2];
    
    //printf("%8.6f %8.6f\n", xrel, yrel);
  }
}

// rotate inverse inertia tensor to rotated reference frame
void rotate_itens(Params p, State s, std::vector<float> &iitens_rot){
  Eigen::Matrix3f arot, iit, iit_rot;
  float qr, qi, qj, qk;
  
  // create rotation matrix from quaternion
  qr = s.orient[0];
  qi = s.orient[1];
  qj = s.orient[2];
  qk = s.orient[3];
  
  arot(0,0) = 1.-2.*(qj*qj+qk*qk);
  arot(1,1) = 1.-2.*(qi*qi+qk*qk);
  arot(2,2) = 1.-2.*(qi*qi+qj*qj);
  
  arot(0,1) = 2*(qi*qj-qk*qr);
  arot(0,2) = 2*(qi*qk+qj*qr);
  arot(1,2) = 2*(qj*qk-qi*qr);
  
  arot(1,0) = 2*(qi*qj+qk*qr);
  arot(2,0) = 2*(qi*qk-qj*qr);
  arot(2,1) = 2*(qj*qk+qi*qr);
  
  // copy tensor vector to matrix
  for (int i = 0; i < 3; i++){
    for (int j = 0; j < 3; j++){
      iit(i,j) = s.iinerm[3*i+j];
    }
    //printf("%d %8.6f %8.6f %8.6f\n", i, s.iinerm[3*i], s.iinerm[3*i+1], s.iinerm[3*i+2]);
  }
  
  // perform transformation
  iit_rot = arot*(iit*arot.transpose());
  
  // copy matrix values to rotated tensor vector
  for (int i = 0; i < 3; i++){
    for (int j = 0; j < 3; j++){
      iitens_rot[3*i+j] = iit_rot(i,j);
    }
    //printf("%d %8.6f %8.6f %8.6f\n", i, iitens_rot[3*i], iitens_rot[3*i+1], iitens_rot[3*i+2]);
  }
}

// get all intersection pairs where distance < ds
void inter_particle_pair(Params p, State s1, State s2, std::vector<int> &ip1, std::vector<int> &ip2, std::vector<float> &dis_vec){
  int k, i, j, npair_bin, sti1, sti2;
  int npar1 = s1.r.size()/3;
  int npar2 = s2.r.size()/3;
  float ds = p.ds;
  float dist, dotp;
  std::vector<float> vp1(3), vp2(3), r12(3);
  
  // clear prior vector data
  ip1.clear();
  ip2.clear();
  dis_vec.clear();
    
  // naive all-pair algorithm
  for (int i = 0; i < npar1; i++){
    for (int j = 0; j < npar2; j++){
    
      dist = pow(s1.rt[3*i]-s2.rt[3*j], 2)+
             pow(s1.rt[3*i+1]-s2.rt[3*j+1], 2)+
             pow(s1.rt[3*i+2]-s2.rt[3*j+2], 2);
      
      if (dist<=4.*ds*ds){
      
        /*
        // calculate particle velocities
        vp1[0] = s1.vel[0]+s1.omega[2]*(s1.rt[3*i+2]-s1.com[2])
                          -s1.omega[3]*(s1.rt[3*i+1]-s1.com[1]);
        vp2[0] = s2.vel[0]+s2.omega[2]*(s2.rt[3*j+2]-s2.com[2])
                          -s2.omega[3]*(s2.rt[3*j+1]-s2.com[1]);
                          
        vp1[1] = s1.vel[1]+s1.omega[3]*(s1.rt[3*i]-s1.com[0])
                          -s1.omega[1]*(s1.rt[3*i+2]-s1.com[2]);
        vp2[1] = s2.vel[1]+s2.omega[3]*(s2.rt[3*j]-s2.com[0])
                          -s2.omega[1]*(s2.rt[3*j+2]-s2.com[2]);
        
        vp1[2] = s1.vel[2]+s1.omega[1]*(s1.rt[3*i+1]-s1.com[1])
                          -s1.omega[2]*(s1.rt[3*i]-s1.com[0]);
        vp2[2] = s2.vel[2]+s2.omega[1]*(s2.rt[3*j+1]-s2.com[1])
                          -s2.omega[2]*(s2.rt[3*j]-s2.com[0]);
                          
        // only consider pairs that are moving towards eachother
        dotp = (vp2[0]-vp1[0])*(s2.rt[3*j]-s1.rt[3*i])
              +(vp2[1]-vp1[1])*(s2.rt[3*j+1]-s1.rt[3*i+1])
              +(vp2[2]-vp1[2])*(s2.rt[3*j+2]-s1.rt[3*i+2]);
        
        if (dotp < 0.){
          ip1.push_back(i);
          ip2.push_back(j);
          dis_vec.push_back(sqrt(dist));
        }
        */
        ip1.push_back(i);
        ip2.push_back(j);
        dis_vec.push_back(sqrt(dist));
      }
    }
  }
}

// test intersection pairs for distance < ds
void test_particle_pair(Params p, State s1, State s2, std::vector<int> ip1, std::vector<int> ip2, bool &intersect){
  int k, i, j, npair_bin, sti1, sti2;
  int npair = ip1.size();
  float ds = p.ds;
  float dist;
  
  // loop over intersection pairs
  intersect = false;
  i = 0;
  
  while ((i < npair)&&(!intersect)){
    dist = pow(s1.rt[3*ip1[i]]-s2.rt[3*ip2[i]], 2)+
           pow(s1.rt[3*ip1[i]+1]-s2.rt[3*ip2[i]+1], 2)+
           pow(s1.rt[3*ip1[i]+2]-s2.rt[3*ip2[i]+2], 2);
           
    if (dist<=ds*ds){
      intersect = true;
    }
    i++;
  }
}
