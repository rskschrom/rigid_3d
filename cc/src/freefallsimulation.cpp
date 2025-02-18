#include <iostream>
#include <iomanip>
#include <fstream>
#include <math.h>
#include <vector>
#include "freefallsimulation.h"
#include "quat.h"
#include "matrix.h"
#include "buoyancy.h"
#include "solver.h"

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
        dq = multiplyVecQuat(st.omega, st.orient);

        // update body position and velocity
        for (int i = 0; i < 4; i++){
            st.orient[i] += dt*0.5*dq[i];
        }
  
        for (int i = 0; i < 3; i++){
            st.comPos[i] += dt*st.comVel[i];
        }
        
        // save motion history
        posHistory.insert(posHistory.end(), std::begin(st.comPos), std::end(st.comPos));
        orientHistory.insert(orientHistory.end(), std::begin(st.orient), std::end(st.orient));
        incrSimStep();
    }
}

void FreeFallSimulation::evolveMotion(std::vector<float> torque, int nstepUser)
{
    std::vector<float> dq;
    std::vector<float> inerm;
    Eigen::Matrix3f matInerm, matIInerm;
    Eigen::Vector3f dOmegaBV = Eigen::Vector3f::Zero();
    int nstep;
    float qNorm;
    
    // get inertia moment tensor and inverse
    matInerm = st.getMatInerm();
    matIInerm = st.getMatIInerm();

    // rotate torque and angular to body reference frame
    std::vector<float> omegaBody, torqueBody;
    std::vector<float> orientWorld(4);
    std::vector<float> solVector(7);
    
    omegaBody = vecRotate(st.omega, conj(st.orient));
    torqueBody = vecRotate(torque, conj(st.orient));
    Eigen::Vector3f omegaBV(omegaBody.data()), torqueBV(torqueBody.data());
    Eigen::Vector4f orientV(st.orient.data());
        
    // use simulation nsteps if no value is provided
    if (nstepUser==0){
        nstep = nt;
    }
    else{
        nstep = nstepUser;
    }
  
    for (int i = 0; i < nstep; i++){

        solVector = rigidMotionPCDM(omegaBV, orientV, dOmegaBV,
                                    matInerm, matIInerm, dt, g, rhofGrad, rhob);
        //solVector = rigidMotionRK4(omegaBV, orientV, torqueBV,
        //                           matInerm, matIInerm, dt, g, rhofGrad, rhob);
        //solVector = rigidMotionEB(omegaBV, orientV, torqueBV,
        //                          matInerm, matIInerm, dt, g, rhofGrad, rhob);
                                   
        for (int i = 0; i < 3; i++){
            // apply damping to angular velocity
            omegaBody[i] = (1.-dampFrac)*solVector[i];
            omegaBV(i) = omegaBody[i];
            dOmegaBV(i) = solVector[7+i];
        }
        for (int i = 0; i < 4; i++){
            orientWorld[i] = solVector[3+i];
            orientV(i) = solVector[3+i];
        }        
        
        st.setOrient(orientWorld);
        
        // test quaternion norm
        qNorm = orientWorld[0]*orientWorld[0]+orientWorld[1]*orientWorld[1]+
                orientWorld[2]*orientWorld[2]+orientWorld[3]*orientWorld[3];
        //std::cout << qNorm << "\tqNorm" << std::endl;

        // rotate new angular velocity to world reference frame
        st.setOmega(vecRotate(omegaBody, st.orient));
        //par.setOmega(omegaBody);
        
        // save motion history
        posHistory.insert(posHistory.end(), std::begin(st.comPos), std::end(st.comPos));
        orientHistory.insert(orientHistory.end(), std::begin(st.orient), std::end(st.orient));
        incrSimStep();
    }
}

void FreeFallSimulation::evolveMotionBuoyancy(int nstepUser)
{
    int nstep;

    // rotate torque and angular to body reference frame
    std::vector<float> torqueBuoyancy;
        
    // use simulation nsteps if no value is provided
    if (nstepUser==0){
        nstep = nt;
    }
    else{
        nstep = nstepUser;
    }
  
    for (int i = 0; i < nstep; i++){
    
        // apply buoyancy torque
        torqueBuoyancy = calcBuoyancyTorque(st.getMatInerm(), st.orient, g, rhofGrad, rhob);
        evolveMotion(torqueBuoyancy, 1);

    }
}

void FreeFallSimulation::writeHistories()
{
    writeVector(posHistory, "position.txt");
    writeVector(orientHistory, "orientation.txt");
}
