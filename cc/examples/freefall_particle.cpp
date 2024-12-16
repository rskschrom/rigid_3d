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
    
    float pmass = 920.;
    int nstep = 2500;
    
    // perform freefall simulation with initial velocity and angular velocity
    com = {0.,0.05,0.4};
    vel = {0.,0.,-0.1};
    orient = {0.,0.,0.,1.};
    omega = {0.,0.,0.,-0.5};
    Particle par = Particle("crystal_0000_1.5.txt", com, vel, orient, omega, pmass);
    
    FreeFallSimulation ffsim = FreeFallSimulation(par, nstep, 1.e-3, -9.81);
    
    ffsim.evolveMotionInertial();
    std::cout << ffsim.getSimStep() << "\n";
    ffsim.writeHistories();
}
