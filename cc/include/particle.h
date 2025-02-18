#include <iostream>
#include <iomanip>
#include <fstream>
#include <math.h>
#include <vector>
#include <Eigen/Dense>
#include "io.h"

#ifndef PARTICLE_H
#define PARTICLE_H

/// Class for particle object
class Particle
{   

    public:
    
        /// initial properties
        std::vector<float> relPoints;
        float pointMass;

        /*!
         * Default constructor.
         *
         * \param relPoints the relative point mass positions.
         * \param pointMass the mass of each point.
         */    

        Particle(std::vector<float> relPoints, float pointMass):
            relPoints(relPoints),
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
        Eigen::Matrix3f inertiaMomentTensor();

};

#endif
