#include <iostream>
#include <fstream>
#include <math.h>
#include <vector>
#include <Eigen/Dense>

#ifndef SOLVER_H
#define SOLVER_H

Eigen::Vector4f orientF(Eigen::Vector3f omegaV, Eigen::Vector4f orientV);
Eigen::Vector4f omegaF(Eigen::Vector3f omegaV, Eigen::Vector3f torqueV,
                       Eigen::Matrix3f matInerm, Eigen::Matrix3f matIInerm);
std::vector<float> rigidMotionRK4(Eigen::Vector3f omegaV, Eigen::Vector4f orientV,
                                  Eigen::Vector3f torqueV,
                                  Eigen::Matrix3f matInerm, Eigen::Matrix3f matIInerm, float dt);
#endif
