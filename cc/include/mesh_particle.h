#include <iostream>
#include <iomanip>
#include <fstream>
#include <math.h>
#include <vector>
#include <Eigen/Dense>
#include "io.h"

#ifndef MESH_PARTICLE_H
#define MESH_PARTICLE_H

/// Class for triangulated mesh particle object
class MeshParticle
{

    private:
        Eigen::Matrix3f matInerm;

    public:
    
        /// initial properties
        Eigen::MatrixX3f vertices;
        Eigen::MatrixX3i faces;
        float density;
        
        // computed properties of the mesh
        Eigen::VectorXf faceAreas;
        Eigen::MatrixX3f faceNorms;

        /*!
         * Default constructor.
         *
         * \param vertices the vertices of the mesh.
         * \param faces the indices of each face corresponding to a triangular patch.
         * \param density the density of the particle.
         */    

        MeshParticle(Eigen::MatrixX3f vertices, Eigen::MatrixX3i faces, float density):
                     vertices(vertices),
                     faces(faces),
                     density(density) {}

        /*!
         * Initializing by reorienting the particle to the principal axes of its inertia momentum tensor.
         *
         */
         //void initialize();
         
         /*!
         * Calculate the area and normal vector of each face.
         *
         */
         void calculateFaceAreasNorms();

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
        //Eigen::Matrix3f inertiaMomentTensor();
        
        void setMatInerm(Eigen::Matrix3f imat){ matInerm = imat; }
        
        Eigen::Matrix3f getMatInerm() { return matInerm;}

};

#endif
