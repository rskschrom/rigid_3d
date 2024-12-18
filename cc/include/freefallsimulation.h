#include <iostream>
#include <iomanip>
#include <fstream>
#include <math.h>
#include <vector>
#include "simulation.h"
#include "particle.h"

#ifndef FREEFALLSIMULATION_H
#define FREEFALLSIMULATION_H

/// Class for free-fall simulations
class FreeFallSimulation: public Simulation
{   
    public:
    
        /// particle object
        Particle par;
        
        /// vectors for the history of positions and orientations
        std::vector<float> posHistory;
        std::vector<float> orientHistory;
                 
        /*!
         * Set values for constructor.
         *
         * \param par the particle object.
         * \param nt the number of timesteps.
         * \param dt the timestep (seconds).
         * \param g the gravitational acceleration
         */
        FreeFallSimulation(Particle par, int nt, float dt, float g):
            Simulation(nt, dt, g),
            par(par) {}
            
        /*!
         * Set values for constructor with defaults for simulation
         *
         * \param par the particle object.
         */
        FreeFallSimulation(Particle par):
            Simulation(),
            par(par) {}
            
        /*!
         * Integrate simulation forward with applied torque.
         *
         * \param torque the applied torque
         * \param nstep the number of time steps to integrate
         */    
        void evolveMotion(std::vector<float> torque, int nstep = 0);
        
        /*!
         * Integrate simulation forward with inertial motions.
         *
         * \param nstep the number of time steps to integrate
         */    
        void evolveMotionInertial(int nstep = 0);
        
        /*!
         * Integrate simulation forward with buoyancy torque.
         *
         * \param nstep the number of time steps to integrate
         */    
        void evolveMotionBuoyancy(int nstep = 0);
        
        /*!
         * Write history arrays.
         *
         */
        void writeHistories();
};
#endif
