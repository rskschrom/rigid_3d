#include <iostream>
#include <math.h>
#include <vector>
#include "params.h"
#include "state.h"
#include "util.h"

// rotate and translate body particles
void transform_particles(Params p, State &s){
  int i;
  int npar = s.x.size();
  
  for (i = 0; i < npar; i++){
    s.xt[i] = cos(s.theta)*s.x[i]-sin(s.theta)*s.y[i]+s.xcom;
    s.yt[i] = sin(s.theta)*s.x[i]+cos(s.theta)*s.y[i]+s.ycom;
  }
}

// rotate and translate hash center points; output 2d grids
void transform_hash(Params p, State &s, std::vector<float> &xht2d, std::vector<float> &yht2d){
  int i, j;
  
  for (i = 0; i < s.nhx; i++){
    for (j = 0; j < s.nhy; j++){
      xht2d[i*s.nhy+j] = cos(s.theta)*s.xh[i]-sin(s.theta)*s.yh[j]+s.xcom;
      yht2d[i*s.nhy+j] = sin(s.theta)*s.xh[i]+cos(s.theta)*s.yh[j]+s.ycom;
    }
  }
}

// spatial hashing of points
void spatial_hash(Params p, State &s){
  int i, j, hi, hj, svi;
  int npar = s.xt.size();
  float dhx, dhy;
  float xmin = s.x[0];
  float xmax = s.x[0];
  float ymin = s.y[0];
  float ymax = s.y[0];
  std::vector<int> sort_vec(npar);
  std::vector<int> chunk_ind(s.nhx*s.nhy); // starting index of each chunk
  //std::vector<int> chunk_offs(s.nhx*s.nhy);
  
  // resize arrays
  s.hi.resize(npar);
  s.xh.resize(s.nhx);
  s.yh.resize(s.nhy);
  s.hcount.resize(s.nhx*s.nhy);
  
  // calculate bin widths
  for (i = 0; i < npar; i++){
    if (s.x[i]<xmin){
      xmin = s.x[i];
    }
    if (s.y[i]<ymin){
      ymin = s.y[i];
    }
    if (s.x[i]>xmax){
      xmax = s.x[i];
    }
    if (s.y[i]>ymax){
      ymax = s.y[i];
    }
  }
  
  // set hash indices for each particle
  dhx = (xmax-xmin)/(s.nhx-1);
  dhy = (ymax-ymin)/(s.nhy-1);
  
  for (i = 0; i < npar; i++){
    hi = (s.x[i]-xmin)/dhx;
    hj = (s.y[i]-ymin)/dhy;    
    s.hi[i] = hj+hi*s.nhy;
    s.hcount[hj+hi*s.nhy]++;
    //printf("%d %d\n", hj+hi*s.nhy, int(s.hcount.size()));
  }
  
  // sort particles by hash index  
  for (i = 0; i < s.nhx*s.nhy-1; i++){
    chunk_ind[i+1] = chunk_ind[i]+s.hcount[i];
  }
  
  for (i = 0; i < npar; i++){
    sort_vec[i] = chunk_ind[s.hi[i]];
    chunk_ind[s.hi[i]]++;
    printf("%d %d\n", i, s.hi[i]);
  }
  
  sortWithVector<float>(s.x, sort_vec);
  sortWithVector<float>(s.y, sort_vec);
  sortWithVector<int>(s.hi, sort_vec);
  
  // test sort
  for (i = 0; i < npar; i++){
    printf("%d %d\n", i, s.hi[i]);
  }
  
  // calculate bin center positions
  for (i = 0; i < s.nhx; i++){
    s.xh[i] = xmin+(float(i)+0.5)*dhx;
  }
  for (i = 0; i < s.nhy; i++){
    s.yh[i] = ymin+(float(i)+0.5)*dhy;
  }
}

// get the index of the closest particle to a point
int close_particle_index(std::vector<float> x, std::vector<float> y, float xp, float yp){
  int i;
  int npar = x.size();
  int cpi = 0;
  float dist;
  float min_dist = sqrt(pow(x[0]-xp, 2)+pow(y[0]-yp, 2));
  
  for (i = 0; i < npar; i++){
    dist = sqrt(pow(x[i]-xp, 2)+pow(y[i]-yp, 2));
    if (dist<min_dist){
      min_dist = dist;
      cpi = i;
    }
  }
  
  return cpi;
}

// get closest particle pair
void close_particle_pair(std::vector<float> x1, std::vector<float> y1, std::vector<float> x2, std::vector<float> y2, int &cpi1, int &cpi2){
  int i, j;
  int npar1 = x1.size();
  int npar2 = x2.size();
  
  float dist;
  float min_dist = sqrt(pow(x1[0]-x2[0], 2)+pow(y1[0]-y2[0], 2));
  
  // naive all-pairs implementation
  cpi1 = 0;
  cpi2 = 0;
   
  for (i = 0; i < npar1; i++){
    for (j = 0; j < npar2; j++){
      dist = sqrt(pow(x1[i]-x2[j], 2)+pow(y1[i]-y2[j], 2));
      if (dist<min_dist){
        min_dist = dist;
        cpi1 = i;
        cpi2 = j;
      }
    }
  }
}

// get all intersection pairs where distance < ds
void inter_particle_pair(Params p, State &s1, State &s2, std::vector<int> &ip1, std::vector<int> &ip2, std::vector<float> &dis_vec){
  int i, j;
  int npar1 = s1.x.size();
  int npar2 = s2.x.size();
  int nh1 = s1.nhx*s1.nhy;
  int nh2 = s2.nhx*s2.nhy;
  float ds = p.ds;
  float dist;
  float dsh = std::max(abs(s1.xh[1]-s1.xh[0]), abs(s1.yh[1]-s1.yh[0]));
  float disth; 
  std::vector<float> xht1(nh1), yht1(nh1);
  std::vector<float> xht2(nh2), yht2(nh2);
  std::vector<int> iph1, iph2;
  
  // clear prior vector data
  ip1.clear();
  ip2.clear();
  dis_vec.clear();
  
  // get binned chunks that are close together
  transform_hash(p, s1, xht1, yht1);
  transform_hash(p, s2, xht2, yht2);

  printf("%8.6f\n", dsh);
  
  // get bins that are close to each other  
  for (i = 0; i < nh1; i++){
    for (j = 0; j < nh2; j++){
      disth = pow(xht1[i]-xht2[j], 2)+pow(yht1[i]-yht2[j], 2);
      if (disth<=dsh*dsh){
        iph1.push_back(i);
        iph2.push_back(j);        
      }
    }
  }
  
  printf("%d\n", int(iph1.size()));
  
  // naive all-pairs implementation  
  for (i = 0; i < npar1; i++){
    for (j = 0; j < npar2; j++){
      dist = pow(s1.xt[i]-s2.xt[j], 2)+pow(s1.yt[i]-s2.yt[j], 2);
      if (dist<=ds*ds){
        ip1.push_back(i);
        ip2.push_back(j);
        dis_vec.push_back(sqrt(dist));
      }
    }
  }
}

// get minimum distance between particles on each body
float min_pair_dist(State s1, State s2){
  int cpi1, cpi2;
  float pair_dist;
  
  close_particle_pair(s1.xt, s1.yt, s2.xt, s2.yt, cpi1, cpi2);
  pair_dist = sqrt(pow(s1.xt[cpi1]-s2.xt[cpi2], 2)+pow(s1.yt[cpi1]-s2.yt[cpi2], 2));
  
  return pair_dist;
}
