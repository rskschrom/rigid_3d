#include <iostream>
#include <iomanip>
#include <fstream>
#include <iterator>
#include <math.h>
#include <vector>
#include "freefallsimulation.h"

int main(void){

    std::vector<float> com;
    std::vector<float> vel;
    std::vector<float> orient;
    std::vector<float> omega;
    std::vector<float> torque;
    
    float pmass = 920.*pow(0.1*1.e-3, 3.);
    int nstep = 2500;
    float theta = 0.;
    float pi = 3.14159265;
    
    // perform freefall simulation with initial velocity and angular velocity
    com = {0.,0.05,0.4};
    vel = {0.,0.,-0.1};
    orient = {cosf(theta*pi/180./2.),0.,sinf(theta*pi/180./2.),0.};
    //omega = {0.,0.,-0.5};
    omega = {0.,0.,0.};

    Particle par = Particle("crystal_points.txt", pmass);
    par.setComPos(com);
    par.setComVel(vel);
    par.setOrient(orient);
    par.setOmega(omega);
    
    torque = {0.,0.,-0.004};
    
    FreeFallSimulation ffsim = FreeFallSimulation(par, nstep, 1.e-3, -9.81);
    
    ffsim.evolveMotion(torque, nstep=1);
    ffsim.evolveMotionInertial(nstep=2499);
    ffsim.writeHistories();
    
    std::cout << ffsim.getSimStep() << "\n" << std::endl;
}
