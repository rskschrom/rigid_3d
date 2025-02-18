#include <iostream>
#include <iomanip>
#include <fstream>
#include <math.h>
#include <vector>
#include "simulation.h"
#include "state.h"

#ifndef FREEFALLSIMULATION_H
#define FREEFALLSIMULATION_H

/// Class for free-fall simulations
class FreeFallSimulation: public Simulation
{   
    public:
    
        /// object parameters
        float rhofGrad;
        float rhob;
        float dampFrac;
        State st;
        
        /// vectors for the history of positions and orientations
        std::vector<float> posHistory;
        std::vector<float> orientHistory;
                 
        /*!
         * Set values for constructor.
         *
         * \param st the State object.
         * \param nt the number of timesteps.
         * \param dt the timestep (seconds).
         * \param g the gravitational acceleration
         * \param rhofGrad the linear fluid density gradient
         * \param rhob the body density
         * \param dampFrac the damping fraction on angular velocity
         */
        FreeFallSimulation(State st, int nt, float dt, float g, float rhofGrad, float rhob, float dampFrac):
            Simulation(nt, dt, g),
            rhofGrad(rhofGrad),
            rhob(rhob),
            dampFrac(dampFrac),
            st(st) {}
            
        /*!
         * Set values for constructor with defaults for simulation
         *
         * \param st the State object.
         */
        FreeFallSimulation(State st, float rhofGrad, float rhob, float dampFrac):
            Simulation(),
            rhofGrad(rhofGrad),
            rhob(rhob),
            dampFrac(dampFrac),
            st(st) {}
            
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
