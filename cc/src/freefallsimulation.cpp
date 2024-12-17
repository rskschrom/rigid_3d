#include <iostream>
#include <iomanip>
#include <fstream>
#include <math.h>
#include <vector>
#include "freefallsimulation.h"
#include "quat.h"
#include "matrix.h"

void FreeFallSimulation::evolveMotionInertial(int nstepUser)
{
    std::vector<float> dq;
    int nstep;
    
    // use simulation nsteps if no value is provided
    if (nstepUser==0){
        nstep = nt;
    }
    else{
        nstep = nstepUser;
    }
  
    for (int i = 0; i < nstep; i++){
        dq = multiplyVecQuat(par.omega, par.orient);

        // update body position and velocity
        for (int i = 0; i < 4; i++){
            par.orient[i] += dt*0.5*dq[i];
        }
  
        for (int i = 0; i < 3; i++){
            par.comPos[i] += dt*par.comVel[i];
        }
        
        // save motion history
        posHistory.insert(posHistory.end(), std::begin(par.comPos), std::end(par.comPos));
        orientHistory.insert(orientHistory.end(), std::begin(par.orient), std::end(par.orient));
        incrSimStep();
    }
}

void FreeFallSimulation::evolveMotion(std::vector<float> torque, int nstepUser)
{
    std::vector<float> dq;
    std::vector<float> inerm;
    Eigen::Matrix3f matInerm, matIInerm;
    int nstep;
    
    // get inertia moment tensor and inverse
    matInerm = par.getMatInerm();
    matIInerm = par.getMatIInerm();

    // rotate torque and angular to body reference frame
    std::vector<float> omegaBody, torqueBody;
    omegaBody = vecRotate(par.omega, par.orient);
    torqueBody = vecRotate(torque, par.orient);
    Eigen::Vector3f omegaBV(omegaBody.data()), torqueBV(torqueBody.data());
    Eigen::Vector3f dOmegaBV;
        
    // use simulation nsteps if no value is provided
    if (nstepUser==0){
        nstep = nt;
    }
    else{
        nstep = nstepUser;
    }
  
    for (int i = 0; i < nstep; i++){
        std::cout << i << "time step." << std::endl;
        std::cout << omegaBody[0] << " " << omegaBody[1] << " " << omegaBody[2] << " " << std::endl;
        dq = multiplyVecQuat(par.omega, par.orient);

        // update body position and velocity
        for (int i = 0; i < 4; i++){
            par.orient[i] += dt*0.5*dq[i];
        }
  
        for (int i = 0; i < 3; i++){
            par.comPos[i] += dt*par.comVel[i];
        }
        
        // update angular velocity
        dOmegaBV = matIInerm*(torqueBV-omegaBV.cross(matInerm*omegaBV));
        std::cout << dOmegaBV << std::endl;
        
        for (int i = 0; i < 3; i++){
            omegaBody[i] += dt*dOmegaBV[i];
        }

        // rotate new angular velocity to world reference frame
        par.setOmega(vecRotate(omegaBody, conj(par.orient)));

        //std::cout << omegaBody[0] << " " << omegaBody[1] << " " << omegaBody[2] << " " << std::endl;
        std::cout << torqueBV << std::endl;
        
        // save motion history
        posHistory.insert(posHistory.end(), std::begin(par.comPos), std::end(par.comPos));
        orientHistory.insert(orientHistory.end(), std::begin(par.orient), std::end(par.orient));
        incrSimStep();
    }
}

void FreeFallSimulation::writeHistories()
{
    writeVector(posHistory, "position.txt");
    writeVector(orientHistory, "orientation.txt");
}
