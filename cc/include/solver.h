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
Eigen::Vector4f omegaBuoyF(Eigen::Vector3f omegaV, Eigen::Vector4f orientV,
                           Eigen::Matrix3f matInerm, Eigen::Matrix3f matIInerm,
                           float g, float rhofGrad, float rhob);
std::vector<float> rigidMotionRK4(Eigen::Vector3f omegaV, Eigen::Vector4f orientV,
                                  Eigen::Vector3f torqueV,
                                  Eigen::Matrix3f matInerm, Eigen::Matrix3f matIInerm,
                                  float dt, float g, float rhofGrad, float rhob);
std::vector<float> rigidMotionEB(Eigen::Vector3f omegaV, Eigen::Vector4f orientV,
                                 Eigen::Vector3f torqueV,
                                 Eigen::Matrix3f matInerm, Eigen::Matrix3f matIInerm,
                                 float dt, float g, float rhofGrad, float rhob);
#endif
