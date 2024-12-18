#include <iostream>
#include <fstream>
#include <math.h>
#include <vector>
#include "quat.h"
#include "matrix.h"

/*
// AB2 solver for omega and orient
std::vector<float> rigidMotionAB2(Eigen::Vector3f omegaV, Eigen::Vector4f orientV,
                                  Eigen::Vector3f torqueV,
                                  Eigen::Matrix3f matInerm, Eigen::Matrix3f matIInerm, float dt)
{
    std::vector<float> solVector(7);
    Eigen::Vector4f dq, orientV1, orientV2;

    // step 1
    dq = multiplyVecQuatEigen(omegaV, orientV);
    orientV1 = orientV+dt*0.5*dq;
    omegaV1 = omegaV+dt*matIInerm*(torqueV-omegaV.cross(matInerm*omegaV));
    
    // step 2
    dq = multiplyVecQuatEigen(omegaV1, orientV1);
    orientV2 = orientV1+dt*0.5*dq;
    omegaV2 = omegaV1+dt*matIInerm*(torqueV-omegaV1.cross(matInerm*omegaV1));
    
    // combine results
    omegaV = 
}
*/

// operator for orientation
Eigen::Vector4f orientF(Eigen::Vector3f omegaV, Eigen::Vector4f orientV)
{
    Eigen::Vector4f dOrient;
    
    dOrient = 0.5*multiplyVecQuatEigen(omegaV, orientV);
    
    return dOrient;
}

// operator for angular velocity
Eigen::Vector3f omegaF(Eigen::Vector3f omegaV, Eigen::Vector3f torqueV,
                       Eigen::Matrix3f matInerm, Eigen::Matrix3f matIInerm)
{
    Eigen::Vector3f dOmega;
    
    dOmega = matIInerm*(torqueV-omegaV.cross(matInerm*omegaV));
    
    return dOmega;
}

// RK4 solver for omega and orient
std::vector<float> rigidMotionRK4(Eigen::Vector3f omegaV, Eigen::Vector4f orientV,
                                  Eigen::Vector3f torqueV,
                                  Eigen::Matrix3f matInerm, Eigen::Matrix3f matIInerm, float dt)
{
    std::vector<float> solVector(7);
    Eigen::Vector4f orientV1;
    Eigen::Vector3f omegaV1;
    Eigen::Vector3f ko1, ko2, ko3, ko4;
    Eigen::Vector4f kq1, kq2, kq3, kq4;
    Eigen::Vector3f o2n, o3n, o4n;
    Eigen::Vector4f q2n, q3n, q4n;

    // step 1
    kq1 = orientF(omegaV, orientV);
    ko1 = omegaF(omegaV, torqueV, matInerm,  matIInerm);
    
    // step 2
    q2n = orientV+dt/2.*kq1;
    o2n = omegaV+dt/2.*ko1;
    kq2 = orientF(o2n, q2n);
    ko2 = omegaF(o2n, torqueV, matInerm,  matIInerm);
    
    // step 3
    q3n = orientV+dt/2.*kq2;
    o3n = omegaV+dt/2.*ko2;
    kq3 = orientF(o3n, q3n);
    ko3 = omegaF(o3n, torqueV, matInerm,  matIInerm);
    
    // step 4
    q4n = orientV+dt*kq3;
    o4n = omegaV+dt*ko3;
    kq4 = orientF(o4n, q4n);
    ko4 = omegaF(o4n, torqueV, matInerm,  matIInerm);
    
    // combine results
    omegaV1 = omegaV+dt/6.*(ko1+2*ko2+2*ko3+ko4);
    orientV1 = orientV+dt/6.*(kq1+2*kq2+2*kq3+kq4);
    
    for (int i = 0; i < 3; i++){
        solVector[i] = omegaV1(i);
        solVector[3+i] = orientV1(i);
    }
    solVector[6] = orientV1(3);
    
    return solVector;
}
