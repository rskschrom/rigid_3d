#include <iostream>
#include <iomanip>
#include <fstream>
#include <math.h>
#include <vector>
#include "freefallsimulation.h"
#include "quat.h"

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
        dq = multiply(par.omega, par.orient);

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
    int nstep;
    
    // use simulation nsteps if no value is provided
    if (nstepUser==0){
        nstep = nt;
    }
    else{
        nstep = nstepUser;
    }
  
    for (int i = 0; i < nstep; i++){
        dq = multiply(par.omega, par.orient);

        // update body position and velocity
        for (int i = 0; i < 4; i++){
            par.orient[i] += dt*0.5*dq[i];
        }
  
        for (int i = 0; i < 3; i++){
            par.comPos[i] += dt*par.comVel[i];
        }
        
        // update angular velocity with torque
        for (int i = 0; i < 4; i++){
            par.omega[i] += dt*torque[i];
        }
        
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
