#include <Eigen/Dense>

#ifndef MATRIX_H
#define MATRIX_H

Eigen::Matrix3f matrixInv(Eigen::Matrix3f mat);
Eigen::Matrix3f fvecMat(std::vector<float> fvec);
Eigen::Matrix3f quatToMatrix(std::vector<float> q);

#endif
