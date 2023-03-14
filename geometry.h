#include <iostream>
#include <math.h>
#include <vector>
#include "params.h"
#include "state.h"

void transform_particles(Params p, State &s);
void transform_hash(Params p, State &s, std::vector<float> xh2d, std::vector<float> yh2d);
void spatial_hash(Params p, State &s);
int close_particle_index(std::vector<float> x, std::vector<float> y, float xp, float yp);
void close_particle_pair(std::vector<float> x1, std::vector<float> y1, std::vector<float> x2, std::vector<float> y2, int &cpi1, int &cpi2);
void inter_particle_pair(Params p, State &s1, State &s2, std::vector<int> &ip1, std::vector<int> &ip2, std::vector<float> &dis_vec);
float min_pair_dist(State s1, State s2);
