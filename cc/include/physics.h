#include <vector>
#include "params.h"
#include "state.h"

void evolve_motion_intertial(Params p, State &s, float dt);
void evolve_motion(Params p, State &s, std::vector<float> alpha, std::vector<float> acc, float dt);
void lj_potential(Params p, State s1, State s2, std::vector<float> con,
                  std::vector<float> &f1, std::vector<float> &f2);
void lj_collision(Params p, State &s1, State &s2, std::vector<float> con, std::vector<float> &dvel1, std::vector<float> &dvel2);
void multi_collision(Params p, State s1, State s2, std::vector<float> alpha1, std::vector<float> acc1,
                     std::vector<float> alpha2, std::vector<float> acc2,
                     std::vector<float> &dv1_tot, std::vector<float> &dv2_tot);
void first_collision(Params p, State s1, State s2, std::vector<float> alpha1, std::vector<float> acc1,
                     std::vector<float> alpha2, std::vector<float> acc2,
                     std::vector<float> &dv1_tot, std::vector<float> &dv2_tot);
void single_collision(Params p, State s1, State s2, std::vector<float> alpha1, std::vector<float> acc1,
                     std::vector<float> alpha2, std::vector<float> acc2,
                     std::vector<float> &dv1_tot, std::vector<float> &dv2_tot);
void test_collision(Params p, State &s1, State &s2, std::vector<float> alpha1, std::vector<float> acc1,
                    std::vector<float> alpha, std::vector<float> acc, bool &attached, int &rcount);
