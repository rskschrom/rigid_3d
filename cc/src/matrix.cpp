#include <iostream>
#include <fstream>
#include <string>
#include <math.h>
#include <vector>
#include "matrix.h"

// inverse of 3x3 matrix
Eigen::Matrix3f matrixInv(Eigen::Matrix3f mat)
{
    Eigen::Matrix3f imat;

    // inverse of inertia tensor
    Eigen::EigenSolver<Eigen::Matrix3f> es(mat);  
    Eigen::Matrix3f iD;
    Eigen::Matrix3f eig_vecs = es.eigenvectors().real();
    Eigen::Array<float, 1, 3> ieigv = es.eigenvalues().real().array().inverse();
  
    iD = ieigv.matrix().asDiagonal();
    imat = eig_vecs*iD*eig_vecs.transpose();
  
    return imat;
}

// convert flattened std::vector(9) to 3x3 matrix
Eigen::Matrix3f fvecMat(std::vector<float> fvec)
{
    Eigen::Matrix3f mat;

    for(int i = 0; i < 3; i++){
        for(int j = 0; j < 3; j++){
            mat(i,j) = fvec[3*i+j];
        }
    }
    return mat;
}
