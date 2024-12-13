#include <iostream>
#include <math.h>
#include <vector>
#include "params.h"
#include "state.h"

void body_bbox(Params p, State s, std::vector<float> &bbox);
void transform_particles(Params p, State &s);
void particle_trajectory(Params p, State s, int ip, float dt, std::vector<float> &traj);
void rotate_itens(Params p, State s, std::vector<float> &itens_rot);
void inter_particle_pair(Params p, State s1, State s2, std::vector<int> &ip1, std::vector<int> &ip2, std::vector<float> &dis_vec);
void test_particle_pair(Params p, State s1, State s2, std::vector<int> ip1, std::vector<int> ip2, bool &intersect);
