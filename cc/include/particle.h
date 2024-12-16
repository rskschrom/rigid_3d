#include <iostream>
#include <iomanip>
#include <fstream>
#include <math.h>
#include <vector>
#include "io.h"

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
                 float pointMass):
            relPoints(relPoints),
            comPos(comPos),
            comVel(comVel),
            orient(orient),
            pointMass(pointMass) {}
                 
        /*!
         * Constructor with everthing set to zero except relative points
         * and point masses.
         *
         * \param relPoints the relative point mass positions.
         * \param pointMass the mass of each point.
         */
         Particle(std::vector<float> relPoints, float pointMass):
            relPoints(relPoints),
            comPos(3, 0.),
            comVel(3, 0.),
            orient({0.,0.,0.,1.}),
            pointMass(pointMass) {}  
         
        /*!
         * Constructor with relative points read in from file
         *
         * \param relPointsFile a text file containing the relative
         * point positions.
         * \param pointMass the mass of each point.
         */
         Particle(std::string relPointFile, float pointMass):
            relPoints(readPoints(relPointFile)),
            comPos(3, 0.),
            comVel(3, 0.),
            orient({0.,0.,0.,1.}),
            pointMass(pointMass) {}

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
