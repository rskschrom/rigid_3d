#include <iostream>
#include <math.h>
#include <vector>
#include <Eigen/Dense>

#ifndef QUAT_H
#define QUAT_H

std::vector<float> multiply(std::vector<float> q1, std::vector<float> q2);
std::vector<float> multiplyVecQuat(std::vector<float> v, std::vector<float> q);
Eigen::Vector4f multiplyVecQuat(Eigen::Vector3f v, Eigen::Vector4f q);
std::vector<float> conj(std::vector<float> q);
Eigen::Vector4f conj(Eigen::Vector4f q);
std::vector<float> vecRotate(std::vector<float> v, std::vector<float> q);
Eigen::Vector3f vecRotate(Eigen::Vector3f v, Eigen::Vector4f q);

#endif
