#include <vector>

#ifndef State_H
#define State_H

// state variables for simulation
struct State{
  // per particle vectors
  std::vector<float> r; // base reference frame particle positions
  std::vector<float> rt; // transformed reference frame particle positions

  // body properties
  std::vector<float> com;
  std::vector<float> vel;
  std::vector<float> inerm;
  std::vector<float> iinerm;
  std::vector<float> orient;
  std::vector<float> omega;
  float mass;
  float rho;
  float dm;
};
#endif
