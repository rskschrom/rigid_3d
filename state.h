#include <vector>

#ifndef State_H
#define State_H

// state variables for simulation
struct State{
    std::vector<float> x; // base reference frame particle positions
    std::vector<float> y;
    std::vector<float> xt; // transformed reference frame particle positions
    std::vector<float> yt;
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
