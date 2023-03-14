#include <vector>
#include "params.h"
#include "state.h"

void evolve_motion(Params p, State &s, float alpha, float ax, float ay, float dt);
void lj_potential(Params p, State s1, State s2, float xcon, float ycon,
                  std::vector<float> &fx1, std::vector<float> &fy1,
                  std::vector<float> &fx2, std::vector<float> &fy2);
void first_collision(Params p, State &s1, State &s2, std::vector<int> &ip1, std::vector<int> &ip2, float tcol);
void lj_collision(Params p, State &s1, State &s2, float alpha1, float alpha2, float ax1, float ax2, float ay1, float ay2);
