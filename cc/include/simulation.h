#include <iostream>
#include <iomanip>
#include <fstream>
#include <math.h>
#include <vector>
#include "io.h"

/// Class for generic simulations
class Simulation
{   
    private:
    
        int simStep = 0;
        
    public:
    
        /// parameters
        int nt;
        float dt;
        float g;
    
        /*!
         * Default constructor.
         *
         */    
        Simulation():
            nt(500),
            dt(1.e-2),
            g(-9.81) {}
                 
        /*!
         * Set values for constructor.
         *
         * \param nt the number of timesteps.
         * \param dt the timestep (seconds).
         * \param g the gravitational acceleration
         */    
        Simulation(int nt, float dt, float g):
            nt(nt),
            dt(dt),
            g(g) {}
            
        /*!
         * Increment the simulation step by one.
         *
         */ 
        void incrSimStep();
        
        /*!
         * Get the current simulation step.
         *
         * \return simStep
         */
        int getSimStep();
};
