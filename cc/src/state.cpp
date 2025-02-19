#include "state.h"
#include "matrix.h"

float State::rotationalPotentialEnergy(float g, float rhofGrad, float rhob){
    Eigen::Matrix3f matInerm = getMatInerm();
    Eigen::Matrix3f matInermWorld, matRot;
    float ixx, iyy, izz;
    float rotPE;
    
    // get inertia mometum tensor in world frame
    matRot = quatToMatrix(orient);
    matInermWorld = matRot * (matInerm * matRot.transpose());
    
    ixx = matInermWorld(0,0);
    iyy = matInermWorld(1,1);
    izz = matInermWorld(2,2);
    
    rotPE = -g*rhofGrad/(4.*rhob)*(ixx+iyy-izz);
    
    return rotPE;
}

float State::rotationalKineticEnergy(){
    Eigen::Matrix3f matInerm = getMatInerm();
    float rotKE;
        
    rotKE = 0.5*(matInerm(0,0)*omega[0]*omega[0]+
                 matInerm(1,1)*omega[1]*omega[1]+
                 matInerm(2,2)*omega[2]*omega[2]);
        
    return rotKE;
}

void State::initialize()
{
    // get inertia tensor and inverse
    Eigen::Matrix3f matInerm = getMatInerm();
    Eigen::Matrix3f matIInerm = matrixInv(matInerm);
    
    setMatIInerm(matIInerm);
    
}

void State::write()
{
    writeVector(comPos, "position.txt");
    writeVector(orient, "orientation.txt");
}
