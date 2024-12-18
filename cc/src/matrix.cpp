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

// convert quaternion to rotation matrix
Eigen::Matrix3f quatToMatrix(std::vector<float> q)
{
    Eigen::Matrix3f rmat;
    float w, x, y, z;
    
    w = q[0];
    x = q[1];
    y = q[2];
    z = q[3];
    
    // set matrix values
    rmat(0,0) = 1.-2.*(y*y+z*z);
    rmat(1,1) = 1.-2.*(x*x+z*z);
    rmat(2,2) = 1.-2.*(x*x+y*y);
    
    rmat(0,1) = 2.*x*y-2.*w*z;
    rmat(1,0) = 2.*x*y+2.*w*z;
    rmat(0,2) = 2.*w*y+2.*x*z;
    rmat(2,0) = 2.*x*z-2.*w*y;
    rmat(1,2) = 2.*y*z-2.*w*x;
    rmat(2,1) = 2.*w*x+2.*y*z;
    
    return rmat;
}
