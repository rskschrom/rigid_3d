#include <iostream>
#include <fstream>
#include <math.h>
#include <vector>
#include "quat.h"
#include "matrix.h"
#include "buoyancy.h"

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
    
    dOrient = 0.5*multiplyVecQuat(omegaV, orientV);
    
    return dOrient;
}

// operator for angular velocity
Eigen::Vector3f omegaF(Eigen::Vector3f omegaV, Eigen::Vector3f torqueV,
                       Eigen::Matrix3f matInerm, Eigen::Matrix3f matIInerm)
{
    Eigen::Vector3f dOmega;
    
    //dOmega = matIInerm * (torqueV-omegaV.cross(matInerm * omegaV));
    dOmega = matIInerm * torqueV;
    
    return dOmega;
}

// operator for angular velocity with dynamic buoyancy
Eigen::Vector3f omegaBuoyF(Eigen::Vector3f omegaV, Eigen::Vector4f orientV,
                           Eigen::Matrix3f matInerm, Eigen::Matrix3f matIInerm, float g)
{
    Eigen::Vector3f dOmega, torqueWV, torqueV;
    
    // get buoyancy torque and transform to body reference frame
    torqueWV = calcBuoyancyTorque(matInerm, orientV, g);
    torqueV = vecRotate(torqueWV, conj(orientV));
    
    //dOmega = matIInerm * (torqueV-omegaV.cross(matInerm * omegaV));
    dOmega = matIInerm * torqueV;
    
    return dOmega;
}

// RK4 solver for omega and orient
std::vector<float> rigidMotionRK4(Eigen::Vector3f omegaV, Eigen::Vector4f orientV,
                                  Eigen::Vector3f torqueV,
                                  Eigen::Matrix3f matInerm, Eigen::Matrix3f matIInerm,
                                  float dt, float g)
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
    //ko1 = omegaF(omegaV, torqueV, matInerm,  matIInerm);
    ko1 = omegaBuoyF(omegaV, orientV, matInerm, matIInerm, g);
    
    // step 2
    q2n = orientV+dt/2.*kq1;
    o2n = omegaV+dt/2.*ko1;
    kq2 = orientF(o2n, q2n);
    //ko2 = omegaF(o2n, torqueV, matInerm,  matIInerm);
    ko2 = omegaBuoyF(o2n, q2n, matInerm, matIInerm, g);
    
    // step 3
    q3n = orientV+dt/2.*kq2;
    o3n = omegaV+dt/2.*ko2;
    kq3 = orientF(o3n, q3n);
    //ko3 = omegaF(o3n, torqueV, matInerm,  matIInerm);
    ko3 = omegaBuoyF(o3n, q3n, matInerm, matIInerm, g);
    
    // step 4
    q4n = orientV+dt*kq3;
    o4n = omegaV+dt*ko3;
    kq4 = orientF(o4n, q4n);
    //ko4 = omegaF(o4n, torqueV, matInerm,  matIInerm);
    ko4 = omegaBuoyF(o4n, q4n, matInerm, matIInerm, g);
    
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

// Euler Backward solver for omega and orient
std::vector<float> rigidMotionEB(Eigen::Vector3f omegaV, Eigen::Vector4f orientV,
                                 Eigen::Vector3f torqueV,
                                 Eigen::Matrix3f matInerm, Eigen::Matrix3f matIInerm,
                                 float dt, float g)
{
    std::vector<float> solVector(7);
    Eigen::Vector4f orientV1;
    Eigen::Vector3f omegaV1;
    float qNorm;
    
    // initial step
    orientV1 = orientV+dt*orientF(omegaV, orientV);
    omegaV1 = omegaV+dt*omegaBuoyF(omegaV, orientV, matInerm, matIInerm, g);
    //std::cout << omegaV1 << std::endl;
    
    // iterative solution
    for (int k = 0; k < 1; k++){
        orientV1 = orientV+dt*orientF(omegaV1, orientV1);
        omegaV1 = omegaV+dt*omegaBuoyF(omegaV1, orientV1, matInerm, matIInerm, g);
        //std::cout << omegaV1 << std::endl;
    }
    
    // renorm orientation if necessary
    qNorm = orientV1(0)*orientV1(0)+orientV1(1)*orientV1(1)+
            orientV1(2)*orientV1(2)+orientV1(3)*orientV1(3);
            
    if (abs(qNorm-1.)>1.e-4){
        orientV1 = orientV1/sqrt(qNorm);
    }
    
    for (int i = 0; i < 3; i++){
        solVector[i] = omegaV1(i);
        solVector[3+i] = orientV1(i);
    }
    solVector[6] = orientV1(3);
    
    return solVector;
}
