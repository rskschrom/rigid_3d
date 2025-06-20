#include <iostream>
#include <math.h>
#include <vector>
#include "quat.h"
#include "matrix.h"
#include "particle.h"

// calculate torque based on buoyancy on inertia momentum tensor
std::vector<float> calcBuoyancyTorque(Eigen::Matrix3f matInerm, std::vector<float> orient, float g,
                                      float rhofGrad, float rhob)
{
    std::vector<float> torque(3);
    Eigen::Matrix3f matInermWorld, matRot;
    
    // get inertia mometum tensor in world frame
    matRot = quatToMatrix(orient);
    matInermWorld = matRot * (matInerm * matRot.transpose());
    
    //std::cout << matRot << "\nrotation" << std::endl;
    //std::cout << matInerm << "\nbody" << std::endl;
    //std::cout << matInermWorld << "\nworld" << std::endl;
    
    torque[0] = -g*rhofGrad/rhob*matInermWorld(1,2);
    torque[1] = g*rhofGrad/rhob*matInermWorld(0,2);
    torque[2] = 0.;
    
    //std::cout << torque[1] << "\ttorque y" << par.omega[1] << "\tomega y" << std::endl;
    
    return torque;
}

// calculate torque based on buoyancy on inertia momentum tensor
Eigen::Vector3f calcBuoyancyTorque(Eigen::Matrix3f matInerm, Eigen::Vector4f orientV, float g,
                                   float rhofGrad, float rhob)
{
    Eigen::Vector3f torqueV;
    Eigen::Matrix3f matInermWorld, matRot;
    
    // get inertia mometum tensor in world frame
    matRot = quatToMatrix(orientV);
    matInermWorld = matRot * (matInerm * matRot.transpose());
    
    torqueV(0) = -g*rhofGrad/rhob*matInermWorld(1,2);
    torqueV(1) = g*rhofGrad/rhob*matInermWorld(0,2);
    torqueV(2) = 0.;
    
    return torqueV;
}
