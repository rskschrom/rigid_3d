#include <iostream>
#include <math.h>
#include <vector>
#include "quat.h"
#include "matrix.h"
#include "particle.h"

// calculate torque based on buoyancy on inertia momentum tensor
std::vector<float> calcBuoyancyTorque(Particle par, float g)
{
    std::vector<float> torque(3);
    Eigen::Matrix3f matInerm = par.getMatInerm();
    Eigen::Matrix3f matInermWorld, matRot;
    
    // get inertia mometum tensor in world frame
    matRot = quatToMatrix(par.orient);
    matInermWorld = matRot * (matInerm * matRot.transpose());
    
    //std::cout << matRot << "\nrotation" << std::endl;
    //std::cout << matInerm << "\nbody" << std::endl;
    //std::cout << matInermWorld << "\nworld" << std::endl;
    
    torque[0] = -g*matInermWorld(1,2);
    torque[1] = g*matInermWorld(0,2);
    torque[2] = 0.;
    
    //std::cout << torque[1] << "\ttorque y" << par.omega[1] << "\tomega y" << std::endl;
    
    return torque;
}
