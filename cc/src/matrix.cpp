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

// convert quaternion to rotation matrix
Eigen::Matrix3f quatToMatrix(Eigen::Vector4f q)
{
    Eigen::Matrix3f rmat;
    float w, x, y, z;
    
    w = q(0);
    x = q(1);
    y = q(2);
    z = q(3);
    
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

/*
// transform relative points to absolute points given orientation and com
std::vector<float> transformPoints(std::vector<float> relPoints,
                                   std::vector<float> orient, std::vector<float> comPos)
{
    int npar = relPoints.size()/3;
    float a, b, c;
    std::vector<float> prot(3);
    std::vector<float> absPoints(npar*3);

    // get rotation coefficients
    a = orient[0]*orient[0]-
        orient[1]*orient[1]-
        orient[2]*orient[2]-
        orient[3]*orient[3];
      
    c = -2.*orient[0];
  
    // loop over points
    for (int i = 0; i < npar; i++){
      b = 2.*(relPoints[3*i]*orient[1]+
              relPoints[3*i+1]*orient[2]+
              relPoints[3*i+2]*orient[3]);
            
      // prot = a*p+b*v+c*pxv
      prot[0] = a*relPoints[3*i]+b*orient[1]+
                c*(relPoints[3*i+1]*orient[3]-relPoints[3*i+2]*orient[2]);
              
      prot[1] = a*relPoints[3*i+1]+b*orient[2]+
                c*(relPoints[3*i+2]*orient[1]-relPoints[3*i]*orient[3]);
              
      prot[2] = a*relPoints[3*i+2]+b*orient[3]+
                c*(relPoints[3*i]*orient[2]-relPoints[3*i+1]*orient[1]);
    
      absPoints[3*i] = prot[0]+comPos[0];
      absPoints[3*i+1] = prot[1]+comPos[1];
      absPoints[3*i+2] = prot[2]+comPos[2];
    }
  
    return absPoints;
  
}
*/
// get the sorted eigenvectors of a tensor
Eigen::Matrix3f sortedEigVecs(Eigen::Matrix3f tensor)
{
    
    // get principal axes of inertia momentum tensor
    Eigen::EigenSolver<Eigen::Matrix3f> es(tensor);
    Eigen::Matrix3f eigVecs = es.eigenvectors().real();
    Eigen::Matrix3f permute = Eigen::Matrix3f::Zero();
    
    // sort eigenvectors by decreasing eigenvalue
    float a = es.eigenvalues().real()[0];
    float b = es.eigenvalues().real()[1];
    float c = es.eigenvalues().real()[2];
    float tmpVal;
    int tmpInd;
    std::vector<int> indSort = {0, 1, 2};
    std::vector<float> eigVals = {a, b, c};
        
    if (eigVals[0]>eigVals[1]){
        tmpVal = eigVals[0];
        tmpInd = indSort[0];
        eigVals[0] = eigVals[1];
        eigVals[1] = tmpVal;
        indSort[0] = indSort[1];
        indSort[1] = tmpInd;
    }
    if (eigVals[0]>eigVals[2]){
        tmpVal = eigVals[0];
        tmpInd = indSort[0];
        eigVals[0] = eigVals[2];
        eigVals[2] = tmpVal;
        indSort[0] = indSort[2];
        indSort[2] = tmpInd;
    }
    if (eigVals[1]>eigVals[2]){
        tmpVal = eigVals[1];
        tmpInd = indSort[1];
        eigVals[1] = eigVals[2];
        eigVals[2] = tmpVal;
        indSort[1] = indSort[2];
        indSort[2] = tmpInd;
    }
        
    // permute matrix
    for (int i = 0; i<3; i++){
        permute(i,indSort[i]) = 1.;
    }
    
    eigVecs = eigVecs * permute;
    
    return eigVecs;
}

// subtract center of mass from points
std::vector<float> centerPoints(std::vector<float> points)
{
    int npar = points.size()/3;
    std::vector<float> relPoints(npar*3);
    float comX, comY, comZ;
    
    // calculate center of mass
    comX = 0.;
    comY = 0.;
    comZ = 0.;
    
    for (int i = 0; i < npar; i++){
        comX = comX+points[3*i];
        comY = comY+points[3*i+1];
        comZ = comZ+points[3*i+2];
    }
    
    comX = comX/float(npar);
    comY = comY/float(npar);
    comZ = comZ/float(npar);
  
    // rotate relative points
    for (int i = 0; i < npar; i++){
        relPoints[3*i] = points[3*i]-comX;
        relPoints[3*i+1] = points[3*i+1]-comY;
        relPoints[3*i+2] = points[3*i+2]-comZ;
    }
  
    return relPoints;
  
}

// rotate points centered at origin to the body-relative reference frame along the inertia moment tensor principal axes
std::vector<float> pointsToBodyFrame(std::vector<float> relPoints, Eigen::Matrix3f inermPA)
{
    int npar = relPoints.size()/3;
    std::vector<float> rotPoints(npar*3);
  
    // rotate relative points
    for (int i = 0; i < npar; i++){
        rotPoints[3*i] = inermPA(0,0)*relPoints[3*i]+
                         inermPA(0,1)*relPoints[3*i+1]+
                         inermPA(0,2)*relPoints[3*i+2];
        rotPoints[3*i+1] = inermPA(1,0)*relPoints[3*i]+
                           inermPA(1,1)*relPoints[3*i+1]+
                           inermPA(1,2)*relPoints[3*i+2];
        rotPoints[3*i+2] = inermPA(2,0)*relPoints[3*i]+
                           inermPA(2,1)*relPoints[3*i+1]+
                           inermPA(2,2)*relPoints[3*i+2];
    }
  
    return rotPoints;
  
}