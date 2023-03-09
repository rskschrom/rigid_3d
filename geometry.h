#include <iostream>
#include <math.h>
#include <vector>
#include "state.h"

int close_particle_index(std::vector<float> x, std::vector<float> y, float xp, float yp);
void close_particle_pair(std::vector<float> x1, std::vector<float> y1, std::vector<float> x2, std::vector<float> y2, int &cpi1, int &cpi2);
void inter_particle_pair(std::vector<float> x1, std::vector<float> y1, std::vector<float> x2, std::vector<float> y2,
                         float ds, std::vector<int> &ip1, std::vector<int> &ip2, std::vector<float> &dis_vec);
float min_pair_dist(State s1, State s2);
