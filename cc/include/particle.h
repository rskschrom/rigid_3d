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

    private:
        Eigen::Matrix3f matInerm;

    public:
    
        /// initial properties
        std::vector<float> points;
        float pointMass;
        
        std::vector<float> relPoints;

        /*!
         * Default constructor.
         *
         * \param points the point mass positions.
         * \param pointMass the mass of each point.
         */    

        Particle(std::vector<float> points, float pointMass):
            points(points),
            pointMass(pointMass) { initialize(); }
         
        /*!
         * Constructor with points read in from file
         *
         * \param pointsFile a text file containing the relative
         * point positions.
         * \param pointMass the mass of each point.
         */
         Particle(std::string pointsFile, float pointMass):
            points(readPoints(pointsFile)),
            pointMass(pointMass) { initialize(); }
            
        /*!
         * Initializing by reorienting the particle to the principal axes of its inertia momentum tensor.
         *
         */
         void initialize();

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
        
        void setMatInerm(Eigen::Matrix3f imat){ matInerm = imat; }
        
        Eigen::Matrix3f getMatInerm() { return matInerm;}

};

#endif
