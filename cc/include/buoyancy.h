#include <vector>
#include "particle.h"
#include <Eigen/Dense>

#ifndef BUOYANCY_H
#define BUOYANCY_H

std::vector<float> calcBuoyancyTorque(Eigen::Matrix3f matInerm, std::vector<float> orient, float g);
Eigen::Vector3f calcBuoyancyTorque(Eigen::Matrix3f matInerm, Eigen::Vector4f orientV, float g);

#endif
