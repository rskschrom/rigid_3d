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
        /// tensors (defined in principal axes reference frame of body)
        Eigen::Matrix3f matInerm, matIInerm;

    public:
    
        /// initial properties
        std::vector<float> relPoints;
        std::vector<float> comPos = std::vector<float>(3);
        std::vector<float> comVel = std::vector<float>(3);
        std::vector<float> orient = std::vector<float>(4);
        std::vector<float> omega = std::vector<float>(3);
        float pointMass;

        /*!
         * Default constructor.
         *
         * \param relPoints the relative point mass positions.
         * \param comPos the center of mass of the particle.
         * \param comVel the the velocity of the body.
         * \param orient the quaternion defining the particle orientation.
         * \param orient the quaternion defining the particle angular velocity.
         * \param pointMass the mass of each point.
         */    

        Particle(std::vector<float> relPoints,
                 std::vector<float> comPos,
                 std::vector<float> comVel,
                 std::vector<float> orient,
                 std::vector<float> omega,
                 float pointMass):
            relPoints(relPoints),
            comPos(comPos),
            comVel(comVel),
            orient(orient),
            omega(omega),
            pointMass(pointMass) { initialize(); }
                 
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
            orient({1.,0.,0.,0.}),
            omega(3, 0.),
            pointMass(pointMass) { initialize(); }  
         
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
            orient({1.,0.,0.,0.}),
            omega(3, 0.),
            pointMass(pointMass) { initialize(); }
            
        /*!
         * Constructor with relative points read in from file
         *
         * \param relPointsFile a text file containing the relative
         * point positions.
         * \param pointMass the mass of each point.
         */
         Particle(std::string relPointFile,
                 std::vector<float> comPos,
                 std::vector<float> comVel,
                 std::vector<float> orient,
                 std::vector<float> omega,
                 float pointMass):
            relPoints(readPoints(relPointFile)),
            comPos(comPos),
            comVel(comVel),
            orient(orient),
            omega(omega),
            pointMass(pointMass) { initialize(); }

        /*!
         * Calculate the particle tensors and set base
         * orientation to their principal axes
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
        std::vector<float> inertiaMomentTensor();
        
        /*!
         * Calculate the rotational potential energy of the particle based on buoyancy.
         *
         * \param g the gravitational acceleration
         * \param rhofGrad the linear density gradient of the fluid
         * \param rhob the density of the body
         *
         * \return the potential energy.
         */
        float rotationalPotentialEnergy(float g, float rhofGrad, float rhob);
        
        /*!
         * Calculate the rotational kinetic energy of the particle.
         *
         * \return the kinetic energy.
         */
        float rotationalKineticEnergy();
        
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

        /*!
         * Set the inertia momentum tensor
         *
         * \param imat inertia momentum tensor.
         */
        void setMatInerm(Eigen::Matrix3f imat){ matInerm = imat; }

        /*!
         * Set the inverse inertia momentum tensor
         *
         * \param iimat inverse inertia momentum tensor.
         */
        void setMatIInerm(Eigen::Matrix3f iimat){ matIInerm = iimat; }

        /*!
         * Set particle center of mass position
         *
         * \param pos
         */
        void setComPos(std::vector<float> pos){ comPos = pos; }

        /*!
         * Set particle velocity
         *
         * \param vel
         */
        void setComVel(std::vector<float> vel){ comVel = vel; }

        /*!
         * Set particle orientation
         *
         * \param ori
         */
        void setOrient(std::vector<float> ori){ orient = ori; }

        /*!
         * Set particle angular velocity
         *
         * \param ome
         */
        void setOmega(std::vector<float> ome){ omega = ome; }

        /*!
         * Get the inertia momentum tensor
         *
         * \return matInerm inertia momentum tensor.
         */
        Eigen::Matrix3f getMatInerm() { return matInerm; }

        /*!
         * Get the inverse inertia momentum tensor
         *
         * \return matIInerm inverse inertia momentum tensor.
         */
        Eigen::Matrix3f getMatIInerm(){ return matIInerm; }
};

#endif
