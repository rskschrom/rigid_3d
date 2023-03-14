#include <vector>

#ifndef State_H
#define State_H

// state variables for simulation
struct State{
  // per particle vectors
  std::vector<float> x; // base reference frame particle positions
  std::vector<float> y;
  std::vector<float> xt; // transformed reference frame particle positions
  std::vector<float> yt;
  std::vector<int> hi; // hash index
  
  // per hash vectors
  std::vector<float> xh; // bin positions
  std::vector<float> yh;  
  std::vector<int> hcount; // count of particles in bin
  
  int nhx; // number of bins for hashing
  int nhy; 
  
  // body properties
  float mass;
  float inerm;
  float xcom;
  float ycom;
  float vx;
  float vy;
  float theta;
  float omega;
  float rho;
  float dm;
};
#endif
