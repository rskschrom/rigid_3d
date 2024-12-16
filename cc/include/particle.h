#include <iostream>
#include <iomanip>
#include <fstream>
#include <math.h>
#include <vector>
#include "params.h"
#include "state.h"
#include "quat.h"
#include "geometry.h"

/// Class for particle object
class Particle
{   
    public:
    
        /// initial properties
        std::vector<float> relPoints;
        std::vector<float> comPos = std::vector<float>(3);
        std::vector<float> comVel = std::vector<float>(3);
        std::vector<float> orient = std::vector<float>(4);
        float pointMass;
    
        /*!
         * Default constructor.
         *
         * \param relPoints the relative point mass positions.
         * \param comPos the center of mass of the particle.
         * \param comVel the the velocity of the body.
         * \param orient the quaternion defining the particle orientation.
         * \param pointMass the mass of each point.
         */    

        Particle(std::vector<float> relPoints,
                 std::vector<float> comPos,
                 std::vector<float> comVel,
                 std::vector<float> orient,
                 float pointMass);
                 
        /*!
         * Constructor with everthing set to zero except relative points
         * and point masses.
         *
         * \param relPoints the relative point mass positions.
         * \param pointMass the mass of each point.
         */
         Particle(std::vector<float> relPoints, float pointMass)
         {
             std::vector<float> comPos(3, 0.);
             std::vector<float> comVel(3, 0.);
             std::vector<float> orient = {0.,0.,0.,1.};     
         }
        
        /*!
         * Calculate the total particle mass.
         *
         * \return the total particle mass.
         */
        float totalMass();
        
        /*!
         * Calculate the particle inertia momentum tensor.
         *
         * \return the particle inertia momentum tensor.
         */
        std::vector<float> inertiaMomentTensor();
        
        /*!
         * Transform relative points to absolute points.
         *
         * \return the absolute point positions.
         */
        std::vector<float> transformPoints();
        
        /*!
         * Write out all particle state information.
         *
         * \return.
         */
        void write();

};
