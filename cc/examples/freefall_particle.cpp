#include <iostream>
#include <iomanip>
#include <fstream>
#include <iterator>
#include <math.h>
#include <vector>
#include "freefallsimulation.h"
#include "buoyancy.h"
#include "quat.h"

int main(void){

    std::vector<float> com;
    std::vector<float> vel;
    std::vector<float> orient;
    std::vector<float> omega, omegaB;
    std::vector<float> torque;
    
    float g = -9.81;
    float rhob = 920.;
    float rhofGrad = 0.1;
    float pmass = rhob*pow(0.1*1.e-3, 3.);
    int nstep = 2500;
    float theta = 10.;
    float pi = 3.14159265;
    
    // perform freefall simulation with initial velocity and angular velocity
    com = {0.,0.05,0.4};
    vel = {0.,0.,0.};
    orient = {cosf(theta*pi/180./2.),0.,sinf(theta*pi/180./2.),0.};
    //omega = {0.,0.,10.};
    omega = vecRotate({0.,0.,0.2}, orient);
    omegaB = vecRotate(omega, conj(orient));

    std::cout << omega[0] << "\t" << omega[1] << "\t" << omega[2] << "\t"
              << "\tomega world" << std::endl;
              
    std::cout << omegaB[0] << "\t" << omegaB[1] << "\t" << omegaB[2] << "\t"
              << "\tomega body" << std::endl;

    Particle par = Particle("crystal_points.txt", pmass);
    par.setComPos(com);
    par.setComVel(vel);
    par.setOrient(orient);
    par.setOmega(omega);
    
    writeVector(par.relPoints, "r.txt");
    
    FreeFallSimulation ffsim = FreeFallSimulation(par, nstep, 1.e-3, g, rhofGrad, rhob);

    //torque = {0.,0.,-0.000004};
    //ffsim.evolveMotion(torque, nstep=1);
    ffsim.evolveMotionBuoyancy(nstep=20000);
    //ffsim.evolveMotionInertial(nstep=999);
    //ffsim.evolveMotion(torque, nstep=1);
    //ffsim.evolveMotionInertial(nstep=20000);
    ffsim.writeHistories();
    
    std::cout << ffsim.getSimStep() << "\n" << std::endl;
}
