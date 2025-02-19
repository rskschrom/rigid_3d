#include <Eigen/Dense>

#ifndef MATRIX_H
#define MATRIX_H

Eigen::Matrix3f matrixInv(Eigen::Matrix3f mat);
Eigen::Matrix3f fvecMat(std::vector<float> fvec);
Eigen::Matrix3f quatToMatrix(std::vector<float> q);
Eigen::Matrix3f quatToMatrix(Eigen::Vector4f q);
//std::vector<float> transformPoints(std::vector<float> relPoints,
//                                   std::vector<float> orient, std::vector<float> comPos);
Eigen::Matrix3f sortedEigVecs(Eigen::Matrix3f tensor);
std::vector<float> centerPoints(std::vector<float> points);
std::vector<float> pointsToBodyFrame(std::vector<float> points, Eigen::Matrix3f inermPA);
#endif
