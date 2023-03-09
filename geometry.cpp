#include <iostream>
#include <math.h>
#include <vector>
#include "state.h"

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
void inter_particle_pair(std::vector<float> x1, std::vector<float> y1, std::vector<float> x2, std::vector<float> y2,
                         float ds, std::vector<int> &ip1, std::vector<int> &ip2, std::vector<float> &dis_vec){
  int i, j;
  int npar1 = x1.size();
  int npar2 = x2.size();
  float dist;
  
  // clear prior vector data
  ip1.clear();
  ip2.clear();
  dis_vec.clear();
  
  // naive all-pairs implementation  
  for (i = 0; i < npar1; i++){
    for (j = 0; j < npar2; j++){
      dist = sqrt(pow(x1[i]-x2[j], 2)+pow(y1[i]-y2[j], 2));
      if (dist<=ds){
        ip1.push_back(i);
        ip2.push_back(j);
        dis_vec.push_back(dist);
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
