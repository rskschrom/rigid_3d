#include <iostream>
#include <iomanip>
#include <fstream>
#include <math.h>
#include <vector>
#include <Eigen/Dense>
#include "io.h"

#ifndef STATE_H
#define STATE_H

/// Class for state object
class State
{   

    private:
        /// tensors (defined in principal axes reference frame of body)
        Eigen::Matrix3f matInerm, matIInerm;

    public:
    
        /// initial properties
        std::vector<float> comPos = std::vector<float>(3);
        std::vector<float> comVel = std::vector<float>(3);
        std::vector<float> orient = std::vector<float>(4);
        std::vector<float> omega = std::vector<float>(3);

        /*!
         * Default constructor.
         *
         * \param matInerm the inertia momentum tensor of the body
         * \param comPos the center of mass of the particle.
         * \param comVel the the velocity of the body.
         * \param orient the quaternion defining the particle orientation.
         * \param orient the quaternion defining the particle angular velocity.
         */    

        State(Eigen::Matrix3f matInerm,
              std::vector<float> comPos,
              std::vector<float> comVel,
              std::vector<float> orient,
              std::vector<float> omega):
            matInerm(matInerm),
            comPos(comPos),
            comVel(comVel),
            orient(orient),
            omega(omega) { initialize(); }
                 
        /*!
         * Constructor with everthing set to zero except inertia momentum tensor.
         *
         * \param matInerm the inertia momentum tensor of the body
         */
         State(Eigen::Matrix3f matInerm):
               matInerm(matInerm),
               comPos(3, 0.),
               comVel(3, 0.),
               orient({1.,0.,0.,0.}),
               omega(3, 0.) { initialize(); }  
         
        /*!
         * Calculate the particle tensors and set base
         * orientation to their principal axes
         *
         */
        void initialize();
        
        /*!
         * Calculate the rotational potential energy based on buoyancy.
         *
         * \param g the gravitational acceleration
         * \param rhofGrad the linear density gradient of the fluid
         * \param rhob the density of the body
         *
         * \return the potential energy.
         */
        float rotationalPotentialEnergy(float g, float rhofGrad, float rhob);
        
        /*!
         * Calculate the rotational kinetic energy.
         *
         * \return the kinetic energy.
         */
        float rotationalKineticEnergy();
        
        /*!
         * Write out all state information.
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
         * Set center of mass position
         *
         * \param pos
         */
        void setComPos(std::vector<float> pos){ comPos = pos; }

        /*!
         * Set velocity
         *
         * \param vel
         */
        void setComVel(std::vector<float> vel){ comVel = vel; }

        /*!
         * Set orientation
         *
         * \param ori
         */
        void setOrient(std::vector<float> ori){ orient = ori; }

        /*!
         * Set angular velocity
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
