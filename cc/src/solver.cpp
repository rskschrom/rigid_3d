#include <iostream>
#include <fstream>
#include <math.h>
#include <vector>
#include "quat.h"
#include "matrix.h"
#include "buoyancy.h"

// operator for orientation
Eigen::Vector4f orientF(Eigen::Vector3f omegaV, Eigen::Vector4f orientV)
{
    Eigen::Vector4f dOrient;
    Eigen::Vector3f omegaWV;
    
    // get omega in world reference frame
    omegaWV = vecRotate(omegaV, orientV);
    
    dOrient = 0.5*multiplyVecQuat(omegaWV, orientV);
    
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
                           Eigen::Matrix3f matInerm, Eigen::Matrix3f matIInerm, float g,
                           float rhofGrad, float rhob)
{
    Eigen::Vector3f dOmega, torqueWV, torqueV;
    
    // get buoyancy torque and transform to body reference frame
    torqueWV = calcBuoyancyTorque(matInerm, orientV, g, rhofGrad, rhob);
    torqueV = vecRotate(torqueWV, conj(orientV));
    
    dOmega = matIInerm * (torqueV-omegaV.cross(matInerm * omegaV));
    //dOmega = matIInerm * torqueV;
    
    return dOmega;
}

// RK4 solver for omega and orient
std::vector<float> rigidMotionRK4(Eigen::Vector3f omegaV, Eigen::Vector4f orientV,
                                  Eigen::Vector3f torqueV,
                                  Eigen::Matrix3f matInerm, Eigen::Matrix3f matIInerm,
                                  float dt, float g, float rhofGrad, float rhob)
{
    std::vector<float> solVector(7);
    Eigen::Vector4f orientV1;
    Eigen::Vector3f omegaV1;
    Eigen::Vector3f ko1, ko2, ko3, ko4;
    Eigen::Vector4f kq1, kq2, kq3, kq4;
    Eigen::Vector3f o2n, o3n, o4n;
    Eigen::Vector4f q2n, q3n, q4n;
    float qNorm;

    // step 1
    kq1 = orientF(omegaV, orientV);
    //ko1 = omegaF(omegaV, torqueV, matInerm,  matIInerm);
    ko1 = omegaBuoyF(omegaV, orientV, matInerm, matIInerm, g, rhofGrad, rhob);
    
    // step 2
    q2n = orientV+dt/2.*kq1;
    o2n = omegaV+dt/2.*ko1;
    kq2 = orientF(o2n, q2n);
    //ko2 = omegaF(o2n, torqueV, matInerm,  matIInerm);
    ko2 = omegaBuoyF(o2n, q2n, matInerm, matIInerm, g, rhofGrad, rhob);
    
    // step 3
    q3n = orientV+dt/2.*kq2;
    o3n = omegaV+dt/2.*ko2;
    kq3 = orientF(o3n, q3n);
    //ko3 = omegaF(o3n, torqueV, matInerm,  matIInerm);
    ko3 = omegaBuoyF(o3n, q3n, matInerm, matIInerm, g, rhofGrad, rhob);
    
    // step 4
    q4n = orientV+dt*kq3;
    o4n = omegaV+dt*ko3;
    kq4 = orientF(o4n, q4n);
    //ko4 = omegaF(o4n, torqueV, matInerm,  matIInerm);
    ko4 = omegaBuoyF(o4n, q4n, matInerm, matIInerm, g, rhofGrad, rhob);
    
    // combine results
    omegaV1 = omegaV+dt/6.*(ko1+2*ko2+2*ko3+ko4);
    orientV1 = orientV+dt/6.*(kq1+2*kq2+2*kq3+kq4);
    
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

// Euler Backward solver for omega and orient
std::vector<float> rigidMotionEB(Eigen::Vector3f omegaV, Eigen::Vector4f orientV,
                                 Eigen::Vector3f torqueV,
                                 Eigen::Matrix3f matInerm, Eigen::Matrix3f matIInerm,
                                 float dt, float g, float rhofGrad, float rhob)
{
    std::vector<float> solVector(7);
    Eigen::Vector4f orientV1;
    Eigen::Vector3f omegaV1;
    float qNorm;
    
    // initial step
    orientV1 = orientV+dt*orientF(omegaV, orientV);
    omegaV1 = omegaV+dt*omegaBuoyF(omegaV, orientV, matInerm, matIInerm, g, rhofGrad, rhob);
    //std::cout << omegaV1 << std::endl;
    
    // iterative solution
    for (int k = 0; k < 1; k++){
        orientV1 = orientV+dt*orientF(omegaV1, orientV1);
        omegaV1 = omegaV+dt*omegaBuoyF(omegaV1, orientV1, matInerm, matIInerm, g, rhofGrad, rhob);
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

// Predictor corrector (Seelen et al. 2016) solver for omega and orient
// note that this solver assumes orientV and omegaV are defined on half steps
//
std::vector<float> rigidMotionPCDM(Eigen::Vector3f omegaBVn12, Eigen::Vector4f orientVn12,
                                   Eigen::Vector3f dOmegaBVn,
                                   Eigen::Matrix3f matInerm, Eigen::Matrix3f matIInerm,
                                   float dt, float g, float rhofGrad, float rhob)
{
    std::vector<float> solVector(10);
    Eigen::Vector4f orientVn1, orientVn32, dOrient;
    Eigen::Vector3f omegaBVn34, omegaWVn34, omegaBVn1, omegaWVn1,
                    torqueWVn1, torqueBVn1, dOmegaBVn1, omegaBVn32;
    float qNorm, omegaMag;
    
    // predict omega at timestep n+3/4
    omegaBVn34 = omegaBVn12+0.25*dt*dOmegaBVn;
    omegaWVn34 = vecRotate(omegaBVn34, orientVn12);
    
    // predict orientation quaternion at timestep n+1
    omegaMag = omegaWVn34(0)*omegaWVn34(0)+
               omegaWVn34(1)*omegaWVn34(1)+
               omegaWVn34(2)*omegaWVn34(2);
    if (omegaMag>1.e-16){
        dOrient(0) = cos(omegaMag*0.25*dt);
        dOrient(1) = sin(omegaMag*0.25*dt)*omegaWVn34(0)/omegaMag;
        dOrient(2) = sin(omegaMag*0.25*dt)*omegaWVn34(1)/omegaMag;
        dOrient(3) = sin(omegaMag*0.25*dt)*omegaWVn34(2)/omegaMag;
    }
    else{
        dOrient(0) = 1.;
        dOrient(1) = 0.;
        dOrient(2) = 0.;
        dOrient(3) = 0.;
    }
    
    orientVn1 = multiply(dOrient, orientVn12);
    //orientVn1(0) = cos(omegaMag*0.25*dt)*orientVn12(0);
    //orientVn1(1) = sin(omegaMag*0.25*dt)*omegaWVn34(0)/omegaMag*orientVn12(1);
    //orientVn1(2) = sin(omegaMag*0.25*dt)*omegaWVn34(1)/omegaMag*orientVn12(2);
    //orientVn1(3) = sin(omegaMag*0.25*dt)*omegaWVn34(2)/omegaMag*orientVn12(3);
    
    // predict angular velocity at timestep n+1
    omegaBVn1 = omegaBVn12+0.5*dt*dOmegaBVn;
    omegaWVn1 = vecRotate(omegaBVn1, orientVn1);
    
    // calculate torque and angular acceleration
    torqueWVn1 = calcBuoyancyTorque(matInerm, orientVn1, g, rhofGrad, rhob);
    torqueBVn1 = vecRotate(torqueWVn1, conj(orientVn1));
    dOmegaBVn1 = matIInerm * (torqueBVn1-omegaBVn1.cross(matInerm * omegaBVn1));
        
    // calculate angular velocity and orientation at timestep n+3/2
    omegaBVn32 = omegaBVn12+dt*dOmegaBVn1;
    
    omegaMag = omegaWVn1(0)*omegaWVn1(0)+
               omegaWVn1(1)*omegaWVn1(1)+
               omegaWVn1(2)*omegaWVn1(2);
    
    if (omegaMag>1.e-16){
        dOrient(0) = cos(omegaMag*0.5*dt);
        dOrient(1) = sin(omegaMag*0.5*dt)*omegaWVn1(0)/omegaMag;
        dOrient(2) = sin(omegaMag*0.5*dt)*omegaWVn1(1)/omegaMag;
        dOrient(3) = sin(omegaMag*0.5*dt)*omegaWVn1(2)/omegaMag;
    }
    else{
        dOrient(0) = 1.;
        dOrient(1) = 0.;
        dOrient(2) = 0.;
        dOrient(3) = 0.;
    }
    
    orientVn32 = multiply(dOrient, orientVn12);
    
    //orientVn32(0) = cos(omegaMag*0.5*dt)*orientVn12(0);
    //orientVn32(1) = sin(omegaMag*0.5*dt)*omegaWVn1(0)/omegaMag*orientVn12(1);
    //orientVn32(2) = sin(omegaMag*0.5*dt)*omegaWVn1(1)/omegaMag*orientVn12(2);
    //orientVn32(3) = sin(omegaMag*0.5*dt)*omegaWVn1(2)/omegaMag*orientVn12(3);
    
    // renorm orientation if necessary
    qNorm = orientVn32(0)*orientVn32(0)+orientVn32(1)*orientVn32(1)+
            orientVn32(2)*orientVn32(2)+orientVn32(3)*orientVn32(3);
    orientVn32 = orientVn32/sqrt(qNorm);       
    //if (abs(qNorm-1.)>1.e-5){
    //    orientVn32 = orientVn32/sqrt(qNorm);
    //}
    
    // set solution vector
    for (int i = 0; i < 3; i++){
        solVector[i] = omegaBVn32(i);
        solVector[3+i] = orientVn32(i);
        solVector[7+i] = dOmegaBVn1(i);
    }
    solVector[6] = orientVn32(3);
    
    return solVector;
}